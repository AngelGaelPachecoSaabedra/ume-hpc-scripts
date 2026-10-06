#!/usr/bin/env python3
"""
build_external_indexes.py — Build derived indexes for PGS platform
====================================================================
Builds fast-lookup indexes from raw data sources.

Usage examples
--------------
# Build everything from scratch
python build_external_indexes.py --all

# Build ClinVar only
python build_external_indexes.py --clinvar

# Build MCPS for all chromosomes (incremental — skips already-built ones)
python build_external_indexes.py --mcps

# Build MCPS for specific chromosomes only
python build_external_indexes.py --mcps --mcps-chroms 22 21 X

# Force rebuild even if done sentinel exists
python build_external_indexes.py --mcps --mcps-chroms 1 --force

# Resume an interrupted build (skips done chromosomes)
python build_external_indexes.py --mcps --resume

# Show build status
python build_external_indexes.py --mcps --status

Outputs
-------
ClinVar: /mnt/cephfs/hot/pgscat/work/clinvar_idx/clinvar_grch38.parquet  (~100MB)
MCPS:    /mnt/cephfs/hot/pgscat/work/mcps_idx/parquet/chrN.parquet       (~20-120MB each)
         /mnt/cephfs/hot/pgscat/work/mcps_idx/done/chrN.done             (sentinels)

Estimated times (v2 — DuckDB native CSV reader)
------------------------------------------------
ClinVar:       1-3 minutes
MCPS per chr:  1-5 seconds  (was 5+ hours/chr with old Python approach)
MCPS all:      ~2-5 minutes total
"""
import argparse
import logging
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
logger = logging.getLogger("build_indexes")

# Add src/ to path so service modules are importable
REPO_ROOT = Path(__file__).parent.parent
SRC_DIR   = REPO_ROOT / "src"
sys.path.insert(0, str(SRC_DIR))

WORK_DIR = Path("/mnt/cephfs/hot/pgscat/work")


# ── ClinVar ───────────────────────────────────────────────────────────────────

def build_clinvar() -> bool:
    from services.clinvar_fetcher import ClinVarFetcher
    logger.info("=" * 60)
    logger.info("Building ClinVar GRCh38 index …")
    cf = ClinVarFetcher(work_dir=WORK_DIR)
    if not cf.is_available():
        logger.error("ClinVar source not found: %s", cf._source)
        return False
    ok = cf._build_index()
    if ok:
        logger.info("ClinVar index ready: %s", cf._parquet)
    else:
        logger.error("ClinVar index build FAILED")
    return ok


# ── MCPS ─────────────────────────────────────────────────────────────────────

def mcps_status() -> None:
    from services.mcps_fetcher import McpsFetcher
    mf = McpsFetcher(work_dir=WORK_DIR)
    s = mf.status()
    logger.info("=" * 60)
    logger.info("MCPS index status")
    logger.info("  Parquet dir : %s", s["parquet_dir"])
    logger.info("  Built       : %d / %d chromosomes", s["n_built"], s["n_total"])
    logger.info("  Complete    : %s", s["complete"])
    logger.info("  Chromosomes : %s", ", ".join(s["built_chroms"]) or "(none)")
    if s["built_chroms"]:
        done_dir = WORK_DIR / "mcps_idx" / "done"
        for c in s["built_chroms"]:
            par = WORK_DIR / "mcps_idx" / "parquet" / f"chr{c}.parquet"
            sz  = par.stat().st_size / 1e6 if par.exists() else 0
            logger.info("    chr%-3s  %.0fMB", c, sz)


def build_mcps(
    chroms=None,
    force: bool = False,
    resume: bool = False,
) -> bool:
    from services.mcps_fetcher import McpsFetcher, _CHROMS
    logger.info("=" * 60)
    logger.info("Building MCPS AF Parquet index (DuckDB native reader) …")
    mf = McpsFetcher(work_dir=WORK_DIR)

    if not mf.is_available():
        logger.error("MCPS source not found: %s", mf._source_dir)
        return False

    # --resume is just "don't force" + default to all chroms
    if resume and chroms is None:
        chroms = _CHROMS
        force = False

    ok = mf.build_index(chroms=chroms, force=force)

    if ok:
        s = mf.status()
        logger.info(
            "MCPS: %d chromosomes ready in %s",
            s["n_built"], s["parquet_dir"],
        )
    else:
        logger.error("MCPS index build encountered errors (check log above)")
    return ok


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build derived lookup indexes for PGS platform",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    # What to build
    parser.add_argument("--clinvar", action="store_true", help="Build ClinVar Parquet index")
    parser.add_argument("--mcps",    action="store_true", help="Build MCPS Parquet index")
    parser.add_argument("--all",     action="store_true", help="Build all indexes")

    # MCPS options
    mcps_grp = parser.add_argument_group("MCPS options")
    mcps_grp.add_argument(
        "--mcps-chroms", nargs="+", default=None, metavar="CHROM",
        help="Build only these chromosomes (e.g. 1 22 X). Accepts with or without 'chr' prefix.",
    )
    mcps_grp.add_argument(
        "--force",  action="store_true",
        help="Rebuild chromosomes even if done sentinel exists",
    )
    mcps_grp.add_argument(
        "--resume", action="store_true",
        help="Resume interrupted build: process all chromosomes, skip already-done ones",
    )
    mcps_grp.add_argument(
        "--status", action="store_true",
        help="Show MCPS build status and exit",
    )

    args = parser.parse_args()

    # --status implies --mcps
    if args.status:
        mcps_status()
        sys.exit(0)

    if not any([args.clinvar, args.mcps, args.all]):
        parser.print_help()
        sys.exit(1)

    results: dict = {}

    if args.clinvar or args.all:
        results["clinvar"] = build_clinvar()

    if args.mcps or args.all:
        results["mcps"] = build_mcps(
            chroms=args.mcps_chroms,
            force=args.force,
            resume=args.resume,
        )

    logger.info("=" * 60)
    logger.info("Summary:")
    for name, ok in results.items():
        logger.info("  %-12s %s", name, "OK" if ok else "FAILED")

    if not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
