#!/usr/bin/env python3
"""
build_regulatory_beds.py — Convert UCSC table-format regulatory files to
sorted, gzip-compressed BED files for use with RegulatoryBEDIndex.
==========================================================================

Sources found in /mnt/cephfs/hot_nvme/ucsc/:
  encRegTfbsClustered.txt.gz    (196 MB) — ENCODE TF binding site clusters
  wgEncodeRegDnaseClustered.txt.gz (71 MB) — DNase I open-chromatin clusters
  oreganno.txt.gz               (20 MB)  — ORegAnno curated regulatory elements
  v3/encodeCcreCombined.bb      (48 MB)  — ENCODE cCRE V3 (BigBed, GRCh38)

UCSC table format:
  All three txt.gz files have a leading "bin" integer column (UCSC B+-tree key).
  The true BED columns start at column 2 (1-indexed):
    col2=chrom, col3=chromStart, col4=chromEnd, col5=name/score/...

cCRE V3 BigBed format:
  Pure-Python reader (scripts/_bigbed_reader.py) — no external tools needed.
  extra[0]=EH38E accession, extra[7]=encodeLabel (PLS/pELS/dELS/CTCF-only/DNase-H3K4me3)

Output BED format (tab-separated):
  UCSC sources: chrom  start  end  name            (BED4, no chr prefix)
  cCRE V3:      chrom  start  end  ccreType  accession  (BED5, no chr prefix, sorted)

Outputs (written to WORK_DIR/regulatory_derived/):
  encodeCcreCombined_GRCh38.bed.gz — ENCODE cCRE V3 (highest priority)
  encode_tfbs_GRCh38.bed.gz       — ENCODE TFBS cluster names (TF identity)
  encode_dhs_GRCh38.bed.gz        — DNase DHS clusters (name="open_chromatin")
  oreganno_GRCh38.bed.gz          — ORegAnno regulatory elements

Usage:
    python scripts/build_regulatory_beds.py [--all] [--ccre] [--tfbs] [--dhs] [--oreganno]
    python scripts/build_regulatory_beds.py --status

Estimated times:
    cCRE   : ~2 min   (48 MB BigBed → ~8 MB BED.GZ)
    TFBS   : ~3-5 min (196 MB → 47 MB)
    DHS    : ~1 min   (71 MB  → 14 MB)
    ORegAnno: ~30 s   (20 MB  → 4 MB)

Build: GRCh38/hg38  (confirmed by chr-prefix + GRCh38 alt scaffolds in data)
"""
import argparse
import gzip
import logging
import os
import sys
import time
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
logger = logging.getLogger("build_reg_beds")

# ── Paths ─────────────────────────────────────────────────────────────────────

UCSC_DIR  = Path("/mnt/cephfs/hot_nvme/ucsc")
WORK_DIR  = Path("/mnt/cephfs/hot/pgscat/work")
OUT_DIR   = WORK_DIR / "regulatory_derived"

# Source files
CCRE_BB_SRC  = UCSC_DIR / "v3" / "encodeCcreCombined.bb"
TFBS_SRC     = UCSC_DIR / "encRegTfbsClustered.txt.gz"
DHS_SRC      = UCSC_DIR / "wgEncodeRegDnaseClustered.txt.gz"
OREGANNO_SRC = UCSC_DIR / "oreganno.txt.gz"

# Output files
CCRE_OUT     = OUT_DIR / "encodeCcreCombined_GRCh38.bed.gz"
TFBS_OUT     = OUT_DIR / "encode_tfbs_GRCh38.bed.gz"
DHS_OUT      = OUT_DIR / "encode_dhs_GRCh38.bed.gz"
OREGANNO_OUT = OUT_DIR / "oreganno_GRCh38.bed.gz"

# Done sentinel paths
CCRE_DONE     = OUT_DIR / "encodeCcreCombined_GRCh38.done"
TFBS_DONE     = OUT_DIR / "encode_tfbs_GRCh38.done"
DHS_DONE      = OUT_DIR / "encode_dhs_GRCh38.done"
OREGANNO_DONE = OUT_DIR / "oreganno_GRCh38.done"

# ── Chromosome normalisation ──────────────────────────────────────────────────

_CHROM_REMAP = {"chrM": "MT", "M": "MT"}

def _strip_chr(name: str) -> str:
    """chr1 → 1, chrX → X, chrM → MT; keeps alt scaffolds."""
    bare = name[3:] if name.startswith("chr") else name
    return _CHROM_REMAP.get(bare, bare)


def _canonical(chrom: str) -> bool:
    """True for canonical chromosomes only (1-22, X, Y, MT)."""
    bare = _strip_chr(chrom)
    canonical = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
    return bare in canonical


# ── cCRE chromosome sort order ────────────────────────────────────────────────

_CHROM_ORDER: dict = {str(i): i for i in range(1, 23)}
_CHROM_ORDER.update({"X": 23, "Y": 24, "MT": 25})


def _ccre_sort_key(rec: tuple) -> tuple:
    return (_CHROM_ORDER.get(rec[0], 99), rec[1], rec[2])


# ── Converters ───────────────────────────────────────────────────────────────

def convert_ccre_bb(src: Path, out: Path, canonical_only: bool = True) -> int:
    """
    Convert ENCODE cCRE V3 BigBed → sorted BED5 (chrom start end ccreType accession).

    Uses scripts/_bigbed_reader.py — no external tools required.

    BigBed extra fields (0-indexed after chrom/start/end):
      extra[0] = EH38E accession  (e.g. EH38E1310153)
      extra[7] = encodeLabel      (PLS / pELS / dELS / CTCF-only / DNase-H3K4me3)

    Records are collected from the R-tree (unordered), sorted by chrom+start,
    then written as gzip-compressed BED5 without chr prefix.

    Returns: number of records written.
    """
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).parent))
    from _bigbed_reader import iter_records  # pure-Python, no deps

    logger.info("Converting ENCODE cCRE V3 BigBed: %s → %s", src.name, out.name)
    t0 = time.time()

    canonical = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
    records: list = []
    n_in = 0
    n_skipped = 0

    for chrom_raw, start, end, extra in iter_records(str(src)):
        n_in += 1
        bare = _strip_chr(chrom_raw)
        if canonical_only and bare not in canonical:
            n_skipped += 1
            continue
        ccre_type = (extra[7].strip() if len(extra) > 7 else "") or "cCRE"
        accession = extra[0].strip() if len(extra) > 0 else ""
        records.append((bare, start, end, ccre_type, accession))
        if n_in % 200_000 == 0:
            logger.info("  cCRE: read %d records …", n_in)

    logger.info("  cCRE: sorting %d records …", len(records))
    records.sort(key=_ccre_sort_key)

    with gzip.open(out, "wt", encoding="utf-8") as fw:
        for bare, start, end, ccre_type, accession in records:
            fw.write(f"{bare}\t{start}\t{end}\t{ccre_type}\t{accession}\n")

    n_out = len(records)
    elapsed = time.time() - t0
    logger.info(
        "  cCRE done: %d / %d records → %s  (%d skipped, %.1f s)",
        n_out, n_in, out.name, n_skipped, elapsed,
    )
    return n_out


def convert_tfbs(src: Path, out: Path, canonical_only: bool = True) -> int:
    """
    Convert UCSC encRegTfbsClustered table → BED3+name.

    UCSC columns (1-indexed):
      1=bin  2=chrom  3=chromStart  4=chromEnd  5=name(TF)  6=score  ...

    Output: chrom(bare)  start  end  TF_name
    Returns: number of records written.
    """
    logger.info("Converting ENCODE TFBS: %s → %s", src.name, out.name)
    n_in = 0
    n_out = 0
    t0 = time.time()

    with gzip.open(src, "rt", encoding="utf-8", errors="replace") as fh, \
         gzip.open(out, "wt", encoding="utf-8") as fw:
        for line in fh:
            n_in += 1
            if n_in % 2_000_000 == 0:
                logger.info("  TFBS: processed %s M lines …", n_in // 1_000_000)
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            # col[0]=bin, col[1]=chrom, col[2]=start, col[3]=end, col[4]=name
            chrom = parts[1]
            if canonical_only and not _canonical(chrom):
                continue
            bare  = _strip_chr(chrom)
            start = parts[2]
            end   = parts[3]
            name  = parts[4].strip() or "TFBS"
            fw.write(f"{bare}\t{start}\t{end}\t{name}\n")
            n_out += 1

    elapsed = time.time() - t0
    logger.info("  TFBS done: %d / %d records → %s  (%.1f s)",
                n_out, n_in, out.name, elapsed)
    return n_out


def convert_dhs(src: Path, out: Path, canonical_only: bool = True) -> int:
    """
    Convert UCSC wgEncodeRegDnaseClustered table → BED3+name.

    UCSC columns (1-indexed):
      1=bin  2=chrom  3=chromStart  4=chromEnd
      5=count(cell_types)  6=score  ...

    Output: chrom(bare)  start  end  open_chromatin
    Returns: number of records written.
    """
    logger.info("Converting ENCODE DHS: %s → %s", src.name, out.name)
    n_in = 0
    n_out = 0
    t0 = time.time()

    with gzip.open(src, "rt", encoding="utf-8", errors="replace") as fh, \
         gzip.open(out, "wt", encoding="utf-8") as fw:
        for line in fh:
            n_in += 1
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[1]
            if canonical_only and not _canonical(chrom):
                continue
            bare  = _strip_chr(chrom)
            start = parts[2]
            end   = parts[3].strip()
            fw.write(f"{bare}\t{start}\t{end}\topen_chromatin\n")
            n_out += 1

    elapsed = time.time() - t0
    logger.info("  DHS done: %d / %d records → %s  (%.1f s)",
                n_out, n_in, out.name, elapsed)
    return n_out


def convert_oreganno(src: Path, out: Path, canonical_only: bool = True) -> int:
    """
    Convert UCSC oreganno table → BED3+name.

    UCSC columns (1-indexed):
      1=bin  2=chrom  3=chromStart  4=chromEnd  5=name(OREG_ID)  6=strand  7=OREG_ID

    Output: chrom(bare)  start  end  OREG_name
    Returns: number of records written.
    """
    logger.info("Converting ORegAnno: %s → %s", src.name, out.name)
    n_in = 0
    n_out = 0
    t0 = time.time()

    with gzip.open(src, "rt", encoding="utf-8", errors="replace") as fh, \
         gzip.open(out, "wt", encoding="utf-8") as fw:
        for line in fh:
            n_in += 1
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            chrom = parts[1]
            if canonical_only and not _canonical(chrom):
                continue
            bare  = _strip_chr(chrom)
            start = parts[2]
            end   = parts[3]
            name  = parts[4].strip() or "regulatory_region"
            fw.write(f"{bare}\t{start}\t{end}\t{name}\n")
            n_out += 1

    elapsed = time.time() - t0
    logger.info("  ORegAnno done: %d / %d records → %s  (%.1f s)",
                n_out, n_in, out.name, elapsed)
    return n_out


# ── Status ────────────────────────────────────────────────────────────────────

def show_status() -> None:
    logger.info("Regulatory BED derived files — status:")
    for label, out, done, src in [
        ("ENCODE cCRE V3", CCRE_OUT,     CCRE_DONE,     CCRE_BB_SRC),
        ("ENCODE TFBS",    TFBS_OUT,     TFBS_DONE,     TFBS_SRC),
        ("ENCODE DHS",     DHS_OUT,      DHS_DONE,      DHS_SRC),
        ("ORegAnno",       OREGANNO_OUT, OREGANNO_DONE, OREGANNO_SRC),
    ]:
        src_exists = src.exists()
        out_exists = out.exists()
        done_exists = done.exists()
        out_mb = out.stat().st_size / 1e6 if out_exists else 0
        logger.info(
            "  %-15s  source=%s  output=%s  done=%s  %.0f MB",
            label,
            "✅" if src_exists else "❌",
            "✅" if out_exists else "❌",
            "✅" if done_exists else "❌",
            out_mb,
        )


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser(
        description="Convert UCSC regulatory table files to BED.GZ for the PGS annotator.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--all",      action="store_true", help="Build all BED files")
    p.add_argument("--ccre",     action="store_true", help="Build ENCODE cCRE V3 BED (from BigBed)")
    p.add_argument("--tfbs",     action="store_true", help="Build ENCODE TFBS BED")
    p.add_argument("--dhs",      action="store_true", help="Build ENCODE DHS (open chromatin) BED")
    p.add_argument("--oreganno", action="store_true", help="Build ORegAnno BED")
    p.add_argument("--force",    action="store_true", help="Overwrite existing outputs")
    p.add_argument("--status",   action="store_true", help="Show build status and exit")
    p.add_argument("--all-chroms", action="store_true",
                   help="Include alt/patch chromosomes (default: canonical only)")
    args = p.parse_args()

    if args.status:
        show_status()
        sys.exit(0)

    if not any([args.all, args.ccre, args.tfbs, args.dhs, args.oreganno]):
        p.print_help()
        sys.exit(1)

    canonical_only = not args.all_chroms
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    results: dict = {}

    if args.ccre or args.all:
        if not CCRE_BB_SRC.exists():
            logger.error("cCRE V3 BigBed not found: %s", CCRE_BB_SRC)
            results["ccre"] = False
        elif CCRE_DONE.exists() and CCRE_OUT.exists() and not args.force:
            logger.info("ENCODE cCRE V3 already built — skipping (use --force to rebuild)")
            results["ccre"] = True
        else:
            try:
                n = convert_ccre_bb(CCRE_BB_SRC, CCRE_OUT, canonical_only=canonical_only)
                CCRE_DONE.write_text(
                    f"n_records={n}\nsource={CCRE_BB_SRC}\n"
                    f"build_time={time.strftime('%Y-%m-%dT%H:%M:%S')}\n"
                )
                results["ccre"] = True
            except Exception as exc:
                logger.error("cCRE V3 conversion FAILED: %s", exc)
                CCRE_OUT.unlink(missing_ok=True)
                CCRE_DONE.unlink(missing_ok=True)
                results["ccre"] = False

    if args.tfbs or args.all:
        if not TFBS_SRC.exists():
            logger.error("TFBS source not found: %s", TFBS_SRC)
            results["tfbs"] = False
        elif TFBS_DONE.exists() and TFBS_OUT.exists() and not args.force:
            logger.info("ENCODE TFBS already built — skipping (use --force to rebuild)")
            results["tfbs"] = True
        else:
            try:
                n = convert_tfbs(TFBS_SRC, TFBS_OUT, canonical_only=canonical_only)
                TFBS_DONE.write_text(
                    f"n_records={n}\nsource={TFBS_SRC}\n"
                    f"build_time={time.strftime('%Y-%m-%dT%H:%M:%S')}\n"
                )
                results["tfbs"] = True
            except Exception as exc:
                logger.error("TFBS conversion FAILED: %s", exc)
                TFBS_OUT.unlink(missing_ok=True)
                TFBS_DONE.unlink(missing_ok=True)
                results["tfbs"] = False

    if args.dhs or args.all:
        if not DHS_SRC.exists():
            logger.error("DHS source not found: %s", DHS_SRC)
            results["dhs"] = False
        elif DHS_DONE.exists() and DHS_OUT.exists() and not args.force:
            logger.info("ENCODE DHS already built — skipping (use --force to rebuild)")
            results["dhs"] = True
        else:
            try:
                n = convert_dhs(DHS_SRC, DHS_OUT, canonical_only=canonical_only)
                DHS_DONE.write_text(
                    f"n_records={n}\nsource={DHS_SRC}\n"
                    f"build_time={time.strftime('%Y-%m-%dT%H:%M:%S')}\n"
                )
                results["dhs"] = True
            except Exception as exc:
                logger.error("DHS conversion FAILED: %s", exc)
                DHS_OUT.unlink(missing_ok=True)
                DHS_DONE.unlink(missing_ok=True)
                results["dhs"] = False

    if args.oreganno or args.all:
        if not OREGANNO_SRC.exists():
            logger.error("ORegAnno source not found: %s", OREGANNO_SRC)
            results["oreganno"] = False
        elif OREGANNO_DONE.exists() and OREGANNO_OUT.exists() and not args.force:
            logger.info("ORegAnno already built — skipping (use --force to rebuild)")
            results["oreganno"] = True
        else:
            try:
                n = convert_oreganno(OREGANNO_SRC, OREGANNO_OUT, canonical_only=canonical_only)
                OREGANNO_DONE.write_text(
                    f"n_records={n}\nsource={OREGANNO_SRC}\n"
                    f"build_time={time.strftime('%Y-%m-%dT%H:%M:%S')}\n"
                )
                results["oreganno"] = True
            except Exception as exc:
                logger.error("ORegAnno conversion FAILED: %s", exc)
                OREGANNO_OUT.unlink(missing_ok=True)
                OREGANNO_DONE.unlink(missing_ok=True)
                results["oreganno"] = False

    logger.info("=" * 60)
    logger.info("Summary:")
    for name, ok in results.items():
        logger.info("  %-12s %s", name, "OK" if ok else "FAILED")

    if results and not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
