"""
dbsnp_freq.py — Population Frequency Fetcher (dbSNP freq.vcf.gz, GRAF-pop format)
===================================================================================
Queries a tabix-indexed dbSNP population-frequency VCF for:
  - rsid                  (e.g. rs123456)
  - af_global             (global allele frequency across all SAMN populations)
  - af_max_population     (highest AF across individual populations)
  - af_population_summary (semicolon-separated key=value pairs)
  - rarity_class          ("common" | "low_frequency" | "rare" | "ultra_rare" | "novel")

File format (NCBI GRAF-pop):
    Allele frequencies are stored in FORMAT + per-sample (SAMN*) columns, NOT in INFO.
    FORMAT fields used: AN (total allele count, including REF), AC (alt count per ALT).
    Each sample column corresponds to one NCBI BioSample population group.

    Header line example:
      #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMN10492695 SAMN10492696 ...

    Data row example (FORMAT=AN:AC:HWEP:GR:GV:GA):
      1  925952  rs2286963  G  A  .  .  .  AN:AC:HWEP:GR:GV:GA  2000:8:0:996:8:0  ...

    Global AF  = sum(AC[alt_idx] for all SAMNs with AN>0) / sum(AN for all SAMNs with AN>0)
    Pop AF     = AC[alt_idx] / AN  (per SAMN column)

Rarity thresholds:
  common        AF >= 0.01      (≥1%)
  low_frequency 0.001 <= AF < 0.01
  rare          0.0001 <= AF < 0.001
  ultra_rare    0 < AF < 0.0001
  novel         AF not available in dbSNP

Usage:
    fetcher = DbSNPFreqFetcher("/dbsnp_freq/freq.vcf.gz")
    result  = fetcher.fetch("1", 925952, "G", "A")
    # {"rsid": "rs2286963", "af_global": 0.004, ...}
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ── GRCh38 chromosome → RefSeq accession (NC_) mapping ───────────────────────
# The dbSNP freq.vcf.gz uses NC_ accession IDs for contig names.
# These are the canonical GRCh38 (hg38) accessions.
_GRCH38_NC: Dict[str, str] = {
    "1":  "NC_000001.11", "2":  "NC_000002.12", "3":  "NC_000003.12",
    "4":  "NC_000004.12", "5":  "NC_000005.10", "6":  "NC_000006.12",
    "7":  "NC_000007.14", "8":  "NC_000008.11", "9":  "NC_000009.12",
    "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
    "13": "NC_000013.11", "14": "NC_000014.9",  "15": "NC_000015.10",
    "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10",
    "19": "NC_000019.10", "20": "NC_000020.11", "21": "NC_000021.9",
    "22": "NC_000022.11",
    "X":  "NC_000023.11", "Y":  "NC_000024.10",
    "MT": "NC_012920.1",  "M":  "NC_012920.1",
}

# ── Rarity classification thresholds ─────────────────────────────────────────
_THRESHOLDS = [
    ("common",        0.01),
    ("low_frequency", 0.001),
    ("rare",          0.0001),
]


def _classify_rarity(af: Optional[float]) -> str:
    """Classify a variant by allele frequency."""
    if af is None or af <= 0.0:
        return "novel"
    for label, threshold in _THRESHOLDS:
        if af >= threshold:
            return label
    return "ultra_rare"


def _parse_format_sample(
    format_str: str,
    sample_str: str,
    alt_idx: int,
) -> Tuple[Optional[int], Optional[int]]:
    """
    Parse FORMAT and one sample column, returning (AN, AC_for_alt_idx).

    Parameters
    ----------
    format_str : str
        Colon-separated FORMAT field order, e.g. "AN:AC:HWEP:GR:GV:GA"
    sample_str : str
        Colon-separated sample values, e.g. "2000:8:0:996:8:0"
    alt_idx : int
        0-based index of the ALT allele we want among all ALT alleles.

    Returns
    -------
    (AN, AC) or (None, None) on parse error / missing data.
    """
    fmt_fields = format_str.split(":")
    smp_values = sample_str.split(":")

    if len(smp_values) < len(fmt_fields):
        return None, None

    try:
        an_idx = fmt_fields.index("AN")
        ac_idx = fmt_fields.index("AC")
    except ValueError:
        return None, None

    try:
        an_val = smp_values[an_idx]
        ac_val = smp_values[ac_idx]
    except IndexError:
        return None, None

    if an_val in (".", "") or ac_val in (".", ""):
        return None, None

    try:
        an = int(an_val)
    except ValueError:
        return None, None

    # AC is Number=A: one value per ALT allele, comma-separated for multi-allelic
    ac_parts = ac_val.split(",")
    try:
        ac = int(ac_parts[alt_idx])
    except (IndexError, ValueError):
        # For mono-allelic rows (no comma), use the only value
        try:
            ac = int(ac_parts[0])
        except (IndexError, ValueError):
            return None, None

    return an, ac


class DbSNPFreqFetcher:
    """
    Tabix-backed fetcher for population allele frequencies (GRAF-pop VCF).

    Reads FORMAT + SAMN* sample columns; does NOT use INFO fields for AF.

    Parameters
    ----------
    vcf_path : Path | str
        Path to the tabix-indexed freq.vcf.gz file.
    """

    # Running counters (reset per annotate_variants.py run via module-level logging)
    n_matched: int = 0
    n_af_found: int = 0
    n_novel: int = 0

    def __init__(self, vcf_path) -> None:
        import pysam  # only required inside Apptainer

        self._path = Path(vcf_path)
        if not self._path.exists():
            raise FileNotFoundError(f"dbSNP freq VCF not found: {self._path}")

        tbi = Path(str(self._path) + ".tbi")
        csi = Path(str(self._path) + ".csi")
        if not tbi.exists() and not csi.exists():
            raise FileNotFoundError(
                f"Tabix index not found for {self._path}. "
                f"Expected {self._path}.tbi or {self._path}.csi"
            )

        self._tbx = pysam.TabixFile(str(self._path))
        self._sample_names: List[str] = []  # ordered SAMN* column names
        self._format_fields_seen: Optional[str] = None  # last FORMAT string seen

        self._parse_header()
        logger.info(
            "DbSNPFreqFetcher opened: %s  |  %d population columns: %s…",
            self._path,
            len(self._sample_names),
            " ".join(self._sample_names[:4]),
        )

    def _parse_header(self) -> None:
        """
        Parse the VCF #CHROM header line to extract ordered SAMN* sample names.
        Columns 0–8 are: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
        Columns 9+ are sample/population names (SAMN*).
        """
        try:
            header_lines = list(self._tbx.header)
        except Exception as exc:
            logger.warning("Cannot read VCF header: %s", exc)
            return

        for line in header_lines:
            if isinstance(line, bytes):
                line = line.decode("utf-8", errors="replace")
            if line.startswith("#CHROM"):
                cols = line.lstrip("#").split("\t")
                # Columns after FORMAT (index 8) are sample names
                if len(cols) > 9:
                    self._sample_names = [c.strip() for c in cols[9:]]
                break

        if not self._sample_names:
            logger.warning(
                "dbSNP VCF: no sample columns found in #CHROM header. "
                "AF computation will be unavailable."
            )

    def fetch(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> Optional[Dict]:
        """
        Query the freq.vcf.gz for population frequencies at chrom:pos.

        Match strategy (in order):
          1. exact_ref_alt_match: VCF REF==ref and ALT contains alt
          2. swapped_alleles_match: caller may have passed (alt, ref) — try swap
          3. multiallelic_alt_match: alt matches one ALT in a multi-allelic record

        Parameters
        ----------
        chrom : str
            Bare chromosome (e.g. "1", "X", "MT"). Tries NC_ accession first.
        pos : int
            1-based position.
        ref : str
            Reference allele (OTHER_ALLELE when IS_FLIP=0, EFFECT_ALLELE when IS_FLIP=1).
        alt : str
            Alternate allele (EFFECT_ALLELE when IS_FLIP=0, OTHER_ALLELE when IS_FLIP=1).

        Returns
        -------
        dict or None
        """
        # Build list of contig keys to try, in priority order:
        # 1. NC_ accession (GRCh38) — canonical key in this VCF
        # 2. Bare number ("1", "X")
        # 3. chr-prefixed ("chr1", "chrX")
        bare = chrom[3:] if chrom.startswith("chr") else chrom
        nc_acc = _GRCH38_NC.get(bare.upper()) or _GRCH38_NC.get(bare)
        chroms_to_try: List[str] = []
        if nc_acc:
            chroms_to_try.append(nc_acc)
        chroms_to_try.append(bare)
        chroms_to_try.append("chr" + bare)

        # pysam.fetch uses 0-based half-open coords
        start0 = pos - 1

        ref_u = ref.upper()
        alt_u = alt.upper()

        for chrom_key in chroms_to_try:
            try:
                raw_rows = list(self._tbx.fetch(chrom_key, start0, pos))
            except (ValueError, KeyError):
                continue

            if not raw_rows:
                continue

            # Decode all rows once; collect candidates at this position
            decoded: List[Tuple[str, ...]] = []
            for row in raw_rows:
                if isinstance(row, bytes):
                    row = row.decode("utf-8", errors="replace")
                parts = row.split("\t")
                if len(parts) < 9:
                    continue
                try:
                    if int(parts[1]) != pos:
                        continue
                except ValueError:
                    continue
                decoded.append(tuple(parts))

            if not decoded:
                continue

            # ── Pass 1: exact match (ref == VCF REF, alt in VCF ALTs) ────────
            match_parts: Optional[Tuple] = None
            match_alt_idx: int = 0
            match_type: str = ""

            for parts in decoded:
                _chr, pos_str, id_field, ref_vcf, alt_vcf = parts[:5]
                if ref_vcf.upper() != ref_u:
                    continue
                vcf_alts = [a.strip() for a in alt_vcf.split(",")]
                alt_upper_list = [a.upper() for a in vcf_alts]
                if alt_u in alt_upper_list:
                    match_parts = parts
                    match_alt_idx = alt_upper_list.index(alt_u)
                    match_type = (
                        "multiallelic_alt_match"
                        if len(vcf_alts) > 1
                        else "exact_ref_alt_match"
                    )
                    logger.debug(
                        "dbSNP %s  %s:%d  %s>%s  id=%s  type=%s  alt_idx=%d",
                        match_type, chrom, pos, ref, alt, parts[2], match_type, match_alt_idx,
                    )
                    break

            # ── Pass 2: swapped alleles (caller sent ref/alt reversed) ────────
            if match_parts is None:
                for parts in decoded:
                    _chr, pos_str, id_field, ref_vcf, alt_vcf = parts[:5]
                    if ref_vcf.upper() != alt_u:
                        continue
                    vcf_alts = [a.strip() for a in alt_vcf.split(",")]
                    alt_upper_list = [a.upper() for a in vcf_alts]
                    if ref_u in alt_upper_list:
                        match_parts = parts
                        match_alt_idx = alt_upper_list.index(ref_u)
                        match_type = "swapped_alleles_match"
                        logger.debug(
                            "dbSNP swapped_alleles_match  %s:%d  "
                            "query_ref=%s query_alt=%s  vcf_ref=%s vcf_alts=%s  id=%s",
                            chrom, pos, ref, alt, ref_vcf, alt_vcf, parts[2],
                        )
                        break

            # ── No match at this position ──────────────────────────────────────
            if match_parts is None:
                if decoded:
                    logger.debug(
                        "dbSNP position_match_but_allele_mismatch  %s:%d  "
                        "query=%s>%s  vcf_records=[%s]",
                        chrom, pos, ref, alt,
                        " | ".join(f"{p[3]}>{p[4]}" for p in decoded[:3]),
                    )
                continue

            # ── Extract rsid and compute AFs ───────────────────────────────────
            parts = match_parts
            id_field = parts[2]
            format_str = parts[8]

            rsid = id_field if id_field not in (".", "", "*") else None
            if rsid:
                logger.debug("dbSNP recovered_rsid  %s:%d  rsid=%s", chrom, pos, rsid)

            # Track FORMAT fields seen
            if self._format_fields_seen != format_str:
                self._format_fields_seen = format_str
                logger.debug("dbSNP FORMAT fields: %s", format_str)

            # ── Compute per-population AFs ────────────────────────────────────
            pop_afs: Dict[str, float] = {}
            total_an = 0
            total_ac = 0

            sample_cols = parts[9:]
            n_samples = min(len(self._sample_names), len(sample_cols))

            for i in range(n_samples):
                sample_name = self._sample_names[i]
                sample_val = sample_cols[i]

                an, ac = _parse_format_sample(format_str, sample_val, match_alt_idx)
                if an is None or ac is None or an == 0:
                    continue

                total_an += an
                total_ac += ac

                pop_af = ac / an
                if 0.0 < pop_af <= 1.0:
                    pop_afs[sample_name] = pop_af

            # ── Global AF ─────────────────────────────────────────────────────
            af_global: Optional[float] = None
            if total_an > 0:
                af_global = total_ac / total_an
                if af_global <= 0.0:
                    af_global = None  # monomorphic in all populations → novel

            if af_global is not None:
                logger.debug(
                    "dbSNP recovered_af  %s:%d  rsid=%s  af_global=%.6g  "
                    "n_pops=%d  match_type=%s",
                    chrom, pos, rsid, af_global, len(pop_afs), match_type,
                )

            # ── Max population AF ──────────────────────────────────────────────
            af_max: Optional[float] = None
            if pop_afs:
                af_max = max(pop_afs.values())
            elif af_global is not None:
                af_max = af_global

            # ── Population summary string ──────────────────────────────────────
            pop_summary: Optional[str] = None
            if pop_afs:
                pop_summary = ";".join(
                    f"{name}={v:.6g}" for name, v in pop_afs.items()
                )

            # ── Rarity + counters + return ────────────────────────────────────
            rarity = _classify_rarity(af_global if af_global is not None else af_max)
            self.n_matched += 1
            if af_global is not None:
                self.n_af_found += 1
            else:
                self.n_novel += 1

            return {
                "rsid":                  rsid,
                "af_global":             af_global,
                "af_max_population":     af_max,
                "af_population_summary": pop_summary,
                "rarity_class":          rarity,
            }

        return None

    def log_summary(self) -> None:
        """Log cumulative match/AF/novel counts. Call after annotation loop."""
        logger.info(
            "dbSNP summary — matched: %d | af_computed: %d | novel: %d | FORMAT: %s",
            self.n_matched,
            self.n_af_found,
            self.n_novel,
            self._format_fields_seen or "unknown",
        )

    def close(self) -> None:
        """Close the tabix file handle."""
        try:
            self._tbx.close()
        except Exception:
            pass

    def __del__(self):
        self.close()
