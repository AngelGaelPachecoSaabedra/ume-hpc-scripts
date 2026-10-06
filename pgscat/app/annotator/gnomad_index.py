#!/usr/bin/env python3
"""
gnomad_index.py — Build per-chromosome parquet index from gnomAD 4.1.1 VCF files.
===================================================================================
One-time admin task (~20-40 min/chrom on CephFS hot_nvme; run in background).
After building, GnomadAnnotator uses DuckDB batch queries instead of per-variant
tabix fetches — eliminating CephFS seek latency for annotation runs.

Produces two parquet files per chromosome (exome takes priority at query time):
    {out_dir}/chrN_exome.parquet
    {out_dir}/chrN_genome.parquet
    {out_dir}/chrN.done           ← sentinel (rows, elapsed, build_time)

Source layout expected (same as run_annotation.sh GNOMAD_BASE):
    {gnomad_dir}/exome/gnomad.exomes.v4.1.1.sites.chrN.vcf.bgz
    {gnomad_dir}/genome/gnomad.genomes.v4.1.1.sites.chrN.vcf.bgz

Usage (inside Apptainer or host with pysam + pyarrow):
    # All chromosomes (1-22, X):
    python gnomad_index.py \\
        --gnomad-dir /gnomad/variants \\
        --out-dir    /gnomad/parquet

    # Specific chromosomes only:
    python gnomad_index.py \\
        --gnomad-dir /gnomad/variants \\
        --out-dir    /gnomad/parquet \\
        --chroms 1 2 X

    # Force rebuild:
    python gnomad_index.py ... --force

Typical output path:
    GNOMAD_PARQUET_DIR=/mnt/cephfs/hot_nvme/gnomad_4.1.1/parquet
"""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import pyarrow as pa
import pyarrow.parquet as pq
import pysam

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
logger = logging.getLogger("gnomad_index")

_CHROMS = [str(i) for i in range(1, 23)] + ["X"]

# Parquet schema — float32 for AF (saves ~50% vs float64, precision sufficient)
_SCHEMA = pa.schema([
    pa.field("pos",     pa.int32()),
    pa.field("ref",     pa.string()),
    pa.field("alt",     pa.string()),
    pa.field("af_raw",  pa.float32()),
    pa.field("af_afr",  pa.float32()),
    pa.field("af_amr",  pa.float32()),
    pa.field("af_eas",  pa.float32()),
    pa.field("af_nfe",  pa.float32()),
    pa.field("af_sas",  pa.float32()),
    pa.field("af_mid",  pa.float32()),
    pa.field("an",      pa.int32()),
    pa.field("nhomalt", pa.int32()),
])

# Write in batches of this many rows to cap peak RAM usage
_BATCH_SIZE = 500_000


# ── VCF INFO helpers ──────────────────────────────────────────────────────────

def _parse_info(info_str: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    for part in info_str.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d


def _getf(info: Dict[str, str], key: str, idx: int = 0) -> Optional[float]:
    v = info.get(key)
    if not v or v == ".":
        return None
    parts = v.split(",")
    try:
        s = parts[idx] if idx < len(parts) else parts[0]
        return float(s) if s not in (".", "") else None
    except (ValueError, IndexError):
        return None


def _geti(info: Dict[str, str], key: str, idx: int = 0) -> Optional[int]:
    v = info.get(key)
    if not v or v == ".":
        return None
    parts = v.split(",")
    try:
        s = parts[idx] if idx < len(parts) else parts[0]
        return int(float(s)) if s not in (".", "") else None
    except (ValueError, IndexError):
        return None


# ── Per-source build ──────────────────────────────────────────────────────────

def _build_source(
    chrom: str,
    vcf_path: Path,
    out_par: Path,
    which: str,
) -> int:
    """
    Stream one gnomAD VCF BGZ (exome or genome) and write a parquet file.

    Returns the number of alt-allele rows written.
    Note: htslib "index file is older than data file" warnings are non-fatal.
    """
    mb = vcf_path.stat().st_size / 1e6
    logger.info("chr%s %s: streaming %.0f MB → %s", chrom, which, mb, out_par.name)
    t0 = time.time()

    batch: list = []
    n_rows = 0
    chrom_vcf = f"chr{chrom}"

    tf = pysam.TabixFile(str(vcf_path))
    try:
        with pq.ParquetWriter(str(out_par), _SCHEMA, compression="snappy") as writer:
            for line in tf.fetch(chrom_vcf):
                parts = line.split("\t")
                if len(parts) < 8:
                    continue
                try:
                    pos = int(parts[1])
                except ValueError:
                    continue

                row_ref  = parts[3].upper()
                vcf_alts = [a.upper() for a in parts[4].split(",")]
                info     = _parse_info(parts[7])

                for i, vcf_alt in enumerate(vcf_alts):
                    batch.append({
                        "pos":     pos,
                        "ref":     row_ref,
                        "alt":     vcf_alt,
                        "af_raw":  _getf(info, "AF", i),
                        "af_afr":  _getf(info, "AF_afr", i),
                        "af_amr":  _getf(info, "AF_amr", i),
                        "af_eas":  _getf(info, "AF_eas", i),
                        "af_nfe":  _getf(info, "AF_nfe", i),
                        "af_sas":  _getf(info, "AF_sas", i),
                        "af_mid":  _getf(info, "AF_mid", i),
                        "an":      _geti(info, "AN", 0),
                        "nhomalt": _geti(info, "nhomalt", i),
                    })
                    n_rows += 1

                if len(batch) >= _BATCH_SIZE:
                    writer.write_batch(
                        pa.RecordBatch.from_pylist(batch, schema=_SCHEMA)
                    )
                    batch = []
                    logger.info("chr%s %s: %s k rows …", chrom, which, f"{n_rows // 1000:,}")

            if batch:
                writer.write_batch(pa.RecordBatch.from_pylist(batch, schema=_SCHEMA))

    finally:
        tf.close()

    elapsed = time.time() - t0
    par_mb  = out_par.stat().st_size / 1e6
    logger.info(
        "chr%s %s: done — %s rows, %.0f MB parquet, %.1f s (%.0f MB/s source)",
        chrom, which, f"{n_rows:,}", par_mb, elapsed, mb / elapsed if elapsed > 0 else 0,
    )
    return n_rows


# ── Per-chromosome build ──────────────────────────────────────────────────────

def build_chrom(
    chrom: str,
    gnomad_dir: Path,
    out_dir: Path,
    force: bool = False,
) -> bool:
    done = out_dir / f"chr{chrom}.done"
    exome_par  = out_dir / f"chr{chrom}_exome.parquet"
    genome_par = out_dir / f"chr{chrom}_genome.parquet"

    if done.exists() and not force:
        logger.info("chr%s: already built — skipping (use --force to rebuild)", chrom)
        return True

    if force:
        for f in (done, exome_par, genome_par):
            f.unlink(missing_ok=True)

    t_total = time.time()
    n_total = 0
    sources_built: List[str] = []

    for which, prefix, out_par in [
        ("exome",  "exomes",  exome_par),
        ("genome", "genomes", genome_par),
    ]:
        vcf = gnomad_dir / which / f"gnomad.{prefix}.v4.1.1.sites.chr{chrom}.vcf.bgz"
        if not vcf.exists():
            logger.warning("chr%s: %s VCF not found — %s", chrom, which, vcf)
            continue
        try:
            n = _build_source(chrom, vcf, out_par, which)
            n_total += n
            sources_built.append(which)
        except Exception as exc:
            logger.error("chr%s %s: FAILED — %s", chrom, which, exc)
            out_par.unlink(missing_ok=True)
            return False

    if not sources_built:
        logger.warning("chr%s: no VCF sources found — skipping", chrom)
        return False

    elapsed = time.time() - t_total
    done.write_text(
        f"chrom={chrom}\n"
        f"sources={','.join(sources_built)}\n"
        f"rows={n_total}\n"
        f"elapsed={elapsed:.1f}\n"
        f"build_time={time.strftime('%Y-%m-%dT%H:%M:%S')}\n"
    )
    logger.info(
        "chr%s: build complete — %s rows total, %.1f s",
        chrom, f"{n_total:,}", elapsed,
    )
    return True


# ── CLI ───────────────────────────────────────────────────────────────────────

def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Build gnomAD 4.1.1 per-chromosome parquet index (one-time admin task)."
    )
    p.add_argument(
        "--gnomad-dir", required=True, type=Path, metavar="DIR",
        help="gnomAD variants dir (contains exome/ and genome/ subdirs)",
    )
    p.add_argument(
        "--out-dir", required=True, type=Path, metavar="DIR",
        help="Output dir for parquet files (set as GNOMAD_PARQUET_DIR in run_annotation.sh)",
    )
    p.add_argument(
        "--chroms", nargs="*", default=None,
        help=f"Chromosomes to build (default: all {len(_CHROMS)}: {' '.join(_CHROMS)})",
    )
    p.add_argument("--force", action="store_true", help="Rebuild even if already built")
    args = p.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    chroms = [c.lstrip("chr") for c in (args.chroms or _CHROMS)]

    logger.info("gnomAD parquet index build")
    logger.info("  source : %s", args.gnomad_dir)
    logger.info("  output : %s", args.out_dir)
    logger.info("  chroms : %s", " ".join(chroms))

    ok = True
    for chrom in chroms:
        ok = build_chrom(chrom, args.gnomad_dir, args.out_dir, force=args.force) and ok

    built = sorted(p.stem.split("_")[0][3:] for p in args.out_dir.glob("chr*_exome.parquet"))
    logger.info(
        "Build complete — %d/%d chromosomes ready%s",
        len(built), len(chroms),
        "" if ok else " (some failed — check logs)",
    )
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
