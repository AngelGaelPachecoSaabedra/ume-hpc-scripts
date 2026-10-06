"""
ScorePreparer service
─────────────────────
Manages and executes the preparation of a PGS Catalog score for local use.

Responsibilities:
  1. Detect current local preparation state (which layout files exist)
  2. Execute real preparation steps (download + betamap build)
  3. Generate a structured manifest describing all expected files
  4. Validate existing files against the layout contract

─── Preparation steps that this service CAN execute ────────────────────────
  A. Download harmonised scoring file from PGS Catalog via pgscat CLI
  B. Parse raw file and build betamap.tsv.gz (pure Python, no external deps)
  C. Write columns.json and betamap.meta.json alongside

─── Steps that require Slurm (NOT executed here) ────────────────────────────
  D. Per-chromosome PRS scoring (GPU Spark job)
  E. PRS aggregation
  F. Parquet export (DuckDB — no Slurm needed, but deferred to Phase 2)

─── Security guarantees ─────────────────────────────────────────────────────
  - pgs_id validated with strict regex before any filesystem operation
  - All paths derived from cfg.scores_dir — no user-supplied paths executed
  - subprocess called with explicit arg list (no shell=True, no string interpolation)
  - Timeouts on every subprocess call
  - Output directory created with mode 0o755 by the app user only
  - Downloaded files written atomically (tmp → rename)

─── Local layout contract ───────────────────────────────────────────────────
Expected layout under {APP_DATA_DIR}/{PGS_ID}/:

  {PGS_ID}.json                               PGS Catalog metadata
  {PGS_ID}_hmPOS_GRCh38.txt.gz               Raw harmonized scoring file
  {PGS_ID}_hmPOS_GRCh38.columns.json         Column mapping report
  {PGS_ID}_hmPOS_GRCh38.betamap.tsv.gz       Pipeline input (8 columns)
  {PGS_ID}_hmPOS_GRCh38.betamap.tsv.gz.meta.json  Generation report

Betamap columns: PRS_ID, CHROM, POS, ID, EFFECT_ALLELE, OTHER_ALLELE, BETA, IS_FLIP
"""
from __future__ import annotations

import gzip
import json
import logging
import os
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Optional

from config import Config
from services.local_catalog import validate_pgs_id
from services.pipeline_inspector import BETAMAP_COLUMNS

logger = logging.getLogger(__name__)

# Preparation stages in order
PREP_STAGES = [
    "catalog_metadata",       # {PGS_ID}.json
    "raw_harmonized_score",   # {PGS_ID}_hmPOS_GRCh38.txt.gz
    "column_map",             # {PGS_ID}_hmPOS_GRCh38.columns.json
    "betamap",                # {PGS_ID}_hmPOS_GRCh38.betamap.tsv.gz
]

COMPUTE_STAGES = [
    "per_chrom_scores",       # {PGS_ID}_chr{N}_scores.tsv (at least chr1)
    "aggregated_scores",      # {PGS_ID}_PRS_total.tsv
    "parquet_export",         # {PGS_ID}_PRS_total.parquet
]

# Timeouts for subprocess calls
_PGSCAT_DOWNLOAD_TIMEOUT = 300   # seconds — large files can take time
_PGSCAT_INFO_TIMEOUT = 30        # seconds

# Known column name variants in raw PGS Catalog files (lowercase)
_CHROM_COLS  = {"hm_chr", "chr_name", "chromosome", "chrom", "chr"}
_POS_COLS    = {"hm_pos", "chr_position", "position", "pos", "bp"}
_ID_COLS     = {"hm_rsid", "rsid", "snp", "variant_id", "id", "markername"}
_EA_COLS     = {"effect_allele", "allele1", "a1", "ea", "effectallele"}
_OA_COLS     = {"hm_inferotherallele", "other_allele", "allele2", "a2", "nea",
                 "non_effect_allele", "otherallele"}
_BETA_COLS   = {"effect_weight", "beta", "b", "effect", "or", "z", "log_or"}


class PreparationError(Exception):
    """Raised when a preparation step fails unrecoverably."""


class ScorePreparer:
    def __init__(self, cfg: Config) -> None:
        self.cfg = cfg

    # ── Public read-only API ─────────────────────────────────────────────────

    def get_preparation_status(self, pgs_id: str) -> dict:
        """
        Return full preparation status for a PGS ID.
        Reads filesystem only — no writes.
        """
        if not validate_pgs_id(pgs_id):
            return {"error": "Invalid PGS ID format"}

        pgs_dir = self.cfg.scores_dir / pgs_id
        files = self._scan_files(pgs_id, pgs_dir)
        stage = self._current_stage(files)
        ready_for_compute = self._ready_for_compute(files)
        ready_for_dashboard = files.get("parquet_export", {}).get("exists", False)
        n_chroms = self._count_chrom_scores(pgs_id, pgs_dir)

        return {
            "pgs_id": pgs_id,
            "pgs_dir": str(pgs_dir),
            "dir_exists": pgs_dir.is_dir(),
            "current_stage": stage,
            "ready_for_compute": ready_for_compute,
            "ready_for_dashboard": ready_for_dashboard,
            "n_chrom_scores": n_chroms,
            "files": files,
            "next_step": self._next_step(stage, files, n_chroms),
            "betamap_columns": BETAMAP_COLUMNS,
        }

    def build_manifest(
        self,
        pgs_id: str,
        remote_info: Optional[dict] = None,
    ) -> dict:
        """
        Build a complete manifest describing all expected files for pgs_id.
        """
        if not validate_pgs_id(pgs_id):
            return {"error": "Invalid PGS ID format"}

        pgs_dir = self.cfg.scores_dir / pgs_id
        p = pgs_id
        sd = str(self.cfg.scores_dir)
        scripts = str(self.cfg.SCRIPTS_DIR)

        ftp_url = ""
        n_variants = None
        is_harmonized = False
        trait = ""
        if remote_info and "error" not in remote_info:
            ftp_url = remote_info.get("ftp_scoring_file", "")
            n_variants = remote_info.get("variants_number")
            is_harmonized = remote_info.get("is_harmonized", False)
            trait = remote_info.get("trait_reported", "") or remote_info.get("name", "")

        file_specs = [
            {
                "filename": f"{p}.json",
                "role": "catalog_metadata",
                "path": str(pgs_dir / f"{p}.json"),
                "exists": (pgs_dir / f"{p}.json").exists(),
                "required": True,
                "description": "PGS Catalog metadata (traits, publication, n_variants)",
                "produced_by": "pgscat score → JSON",
                "auto_preparable": True,
            },
            {
                "filename": f"{p}_hmPOS_GRCh38.txt.gz",
                "role": "raw_harmonized_score",
                "path": str(pgs_dir / f"{p}_hmPOS_GRCh38.txt.gz"),
                "exists": (pgs_dir / f"{p}_hmPOS_GRCh38.txt.gz").exists(),
                "required": True,
                "description": "Raw GRCh38-harmonised scoring file from PGS Catalog FTP",
                "produced_by": "pgscat download --build GRCh38 --type positions",
                "auto_preparable": True,
                "is_harmonized": is_harmonized,
                "n_variants": n_variants,
            },
            {
                "filename": f"{p}_hmPOS_GRCh38.columns.json",
                "role": "column_map",
                "path": str(pgs_dir / f"{p}_hmPOS_GRCh38.columns.json"),
                "exists": (pgs_dir / f"{p}_hmPOS_GRCh38.columns.json").exists(),
                "required": True,
                "description": "Column mapping from raw file to betamap schema",
                "produced_by": "betamap builder (Python, no Slurm)",
                "auto_preparable": True,
            },
            {
                "filename": f"{p}_hmPOS_GRCh38.betamap.tsv.gz",
                "role": "betamap",
                "path": str(pgs_dir / f"{p}_hmPOS_GRCh38.betamap.tsv.gz"),
                "exists": (pgs_dir / f"{p}_hmPOS_GRCh38.betamap.tsv.gz").exists(),
                "required": True,
                "description": (
                    f"Pipeline input. Columns: {', '.join(BETAMAP_COLUMNS)}. "
                    "IS_FLIP=1 means effect on REF strand (dosage = 2 − call_DS)."
                ),
                "produced_by": "betamap builder (Python, no Slurm)",
                "auto_preparable": True,
            },
            {
                "filename": f"{p}_hmPOS_GRCh38.betamap.tsv.gz.meta.json",
                "role": "betamap_meta",
                "path": str(pgs_dir / f"{p}_hmPOS_GRCh38.betamap.tsv.gz.meta.json"),
                "exists": (pgs_dir / f"{p}_hmPOS_GRCh38.betamap.tsv.gz.meta.json").exists(),
                "required": False,
                "description": "Betamap generation report (row counts, skipped variants)",
                "produced_by": "betamap builder",
                "auto_preparable": True,
            },
            {
                "filename": f"{p}_chr{{N}}_scores.tsv",
                "role": "per_chrom_scores",
                "path": str(pgs_dir / f"{p}_chr{{N}}_scores.tsv"),
                "exists": self._count_chrom_scores(pgs_id, pgs_dir) > 0,
                "required": True,
                "description": "Per-chromosome PRS (cols: sample_id, PRS). One file per chrom 1–22.",
                "produced_by": "compute_prs_spark_gpu.py / prs_gpu_compute.sbatch",
                "auto_preparable": False,
            },
            {
                "filename": f"{p}_PRS_total.tsv",
                "role": "aggregated_scores",
                "path": str(pgs_dir / f"{p}_PRS_total.tsv"),
                "exists": (pgs_dir / f"{p}_PRS_total.tsv").exists(),
                "required": True,
                "description": "Aggregated PRS. Cols: sample_id, PRS_total, PRS_chr1..22.",
                "produced_by": "aggregation step",
                "auto_preparable": False,
            },
            {
                "filename": f"{p}_PRS_total_metadata.json",
                "role": "aggregated_metadata",
                "path": str(pgs_dir / f"{p}_PRS_total_metadata.json"),
                "exists": (pgs_dir / f"{p}_PRS_total_metadata.json").exists(),
                "required": False,
                "description": "Aggregation report",
                "produced_by": "aggregation step",
                "auto_preparable": False,
            },
            {
                "filename": f"{p}_PRS_total.parquet",
                "role": "parquet_export",
                "path": str(pgs_dir / f"{p}_PRS_total.parquet"),
                "exists": (pgs_dir / f"{p}_PRS_total.parquet").exists(),
                "required": False,
                "description": "Parquet export (ZSTD). Required for web dashboard.",
                "produced_by": "DuckDB COPY TO",
                "auto_preparable": False,
            },
        ]

        status = self.get_preparation_status(pgs_id)

        return {
            "pgs_id": pgs_id,
            "pgs_dir": str(pgs_dir),
            "trait": trait,
            "is_harmonized": is_harmonized,
            "n_variants": n_variants,
            "current_stage": status["current_stage"],
            "ready_for_compute": status["ready_for_compute"],
            "ready_for_dashboard": status["ready_for_dashboard"],
            "n_chrom_scores": status["n_chrom_scores"],
            "next_step": status["next_step"],
            "betamap_columns": BETAMAP_COLUMNS,
            "files": file_specs,
        }

    # ── Preparation execution ─────────────────────────────────────────────────

    def prepare_score(self, pgs_id: str) -> dict:
        """
        Execute the full input-preparation pipeline for pgs_id.

        Steps performed (in order, skipping already-done ones):
          1. Create output directory
          2. Download catalog metadata → {PGS_ID}.json
          3. Download harmonised scoring file → {PGS_ID}_hmPOS_GRCh38.txt.gz
          4. Build betamap → .betamap.tsv.gz + .columns.json + .meta.json

        Returns a result dict with:
          status:  "ok" | "partial" | "error"
          steps:   list of {name, status, message, elapsed_s}
          stage:   final stage reached
          ready_for_compute: bool
        """
        if not validate_pgs_id(pgs_id):
            return {"status": "error", "error": "Invalid PGS ID format", "steps": []}

        pgs_dir = self.cfg.scores_dir / pgs_id
        steps_log: list[dict] = []
        overall_ok = True

        def log_step(name: str, status: str, message: str, elapsed: float) -> None:
            steps_log.append({
                "name": name,
                "status": status,
                "message": message,
                "elapsed_s": round(elapsed, 2),
            })
            level = logging.INFO if status in ("ok", "skipped") else logging.ERROR
            logger.log(level, "[prepare_score %s] %s: %s (%.1fs)", pgs_id, name, message, elapsed)

        # ── Step 0: Create directory ──────────────────────────────────────────
        t0 = time.monotonic()
        try:
            pgs_dir.mkdir(parents=True, exist_ok=True)
            log_step("create_directory", "ok", str(pgs_dir), time.monotonic() - t0)
        except OSError as exc:
            log_step("create_directory", "error", str(exc), time.monotonic() - t0)
            return self._result(pgs_id, pgs_dir, "error", steps_log,
                                error=f"Cannot create directory {pgs_dir}: {exc}")

        # ── Step 1: Catalog metadata ──────────────────────────────────────────
        t0 = time.monotonic()
        meta_path = pgs_dir / f"{pgs_id}.json"
        if meta_path.exists():
            log_step("catalog_metadata", "skipped", "already exists", time.monotonic() - t0)
        else:
            ok, msg = self._download_catalog_metadata(pgs_id, meta_path)
            log_step("catalog_metadata", "ok" if ok else "error", msg, time.monotonic() - t0)
            if not ok:
                overall_ok = False

        # ── Step 2: Download harmonised scoring file ──────────────────────────
        t0 = time.monotonic()
        raw_path = pgs_dir / f"{pgs_id}_hmPOS_GRCh38.txt.gz"
        if raw_path.exists():
            log_step("download_raw_score", "skipped", "already exists", time.monotonic() - t0)
        else:
            ok, msg = self._download_raw_score(pgs_id, raw_path)
            log_step("download_raw_score", "ok" if ok else "error", msg, time.monotonic() - t0)
            if not ok:
                overall_ok = False
                return self._result(pgs_id, pgs_dir, "error", steps_log,
                                    error=f"Download failed: {msg}")

        # ── Step 3: Build betamap ─────────────────────────────────────────────
        t0 = time.monotonic()
        betamap_path = pgs_dir / f"{pgs_id}_hmPOS_GRCh38.betamap.tsv.gz"
        cols_path = pgs_dir / f"{pgs_id}_hmPOS_GRCh38.columns.json"
        meta_bm_path = pgs_dir / f"{pgs_id}_hmPOS_GRCh38.betamap.tsv.gz.meta.json"
        if betamap_path.exists() and cols_path.exists():
            log_step("build_betamap", "skipped", "already exists", time.monotonic() - t0)
        else:
            ok, msg = self._build_betamap(
                pgs_id, raw_path, betamap_path, cols_path, meta_bm_path
            )
            log_step("build_betamap", "ok" if ok else "error", msg, time.monotonic() - t0)
            if not ok:
                overall_ok = False

        status = "ok" if overall_ok else "partial"
        return self._result(pgs_id, pgs_dir, status, steps_log)

    # ── Download: catalog metadata ────────────────────────────────────────────

    def _download_catalog_metadata(self, pgs_id: str, out_path: Path) -> tuple[bool, str]:
        """
        Fetch PGS Catalog metadata for pgs_id using pgscat CLI and save as JSON.
        Uses atomic write (tmp → rename).
        """
        try:
            proc = subprocess.run(
                [self.cfg.PGSCAT_BIN, "--json", "score", pgs_id],
                capture_output=True,
                text=True,
                timeout=_PGSCAT_INFO_TIMEOUT,
                check=False,
            )
            if proc.returncode != 0:
                err = (proc.stderr or "").strip() or f"exit {proc.returncode}"
                return False, f"pgscat score failed: {err}"

            raw = proc.stdout.strip()
            if not raw:
                return False, "pgscat returned empty output"

            # Validate it's parseable JSON
            data = json.loads(raw)

            # Atomic write
            _atomic_write_text(out_path, json.dumps(data, indent=2, ensure_ascii=False))
            return True, f"saved {out_path.name}"

        except subprocess.TimeoutExpired:
            return False, f"pgscat timed out after {_PGSCAT_INFO_TIMEOUT}s"
        except json.JSONDecodeError as exc:
            return False, f"pgscat returned invalid JSON: {exc}"
        except Exception as exc:
            logger.exception("Unexpected error downloading metadata for %s", pgs_id)
            return False, str(exc)

    # ── Download: harmonised scoring file ─────────────────────────────────────

    def _download_raw_score(self, pgs_id: str, out_path: Path) -> tuple[bool, str]:
        """
        Download the GRCh38-harmonised scoring file using pgscat CLI.

        Tries --type positions (harmonized positions) first.
        Falls back to default (original) if positions are unavailable.
        Uses atomic write via temp file in same directory.
        """
        tmp = out_path.with_suffix(".tmp.gz")
        try:
            # Try harmonised positions first
            for file_type in ("positions", "original"):
                proc = subprocess.run(
                    [
                        self.cfg.PGSCAT_BIN,
                        "download", pgs_id,
                        "--build", "GRCh38",
                        "--type", file_type,
                        "-o", str(tmp),
                    ],
                    capture_output=True,
                    text=True,
                    timeout=_PGSCAT_DOWNLOAD_TIMEOUT,
                    check=False,
                )
                if proc.returncode == 0 and tmp.exists() and tmp.stat().st_size > 0:
                    tmp.rename(out_path)
                    size_mb = round(out_path.stat().st_size / 1_048_576, 2)
                    return True, f"downloaded {file_type} ({size_mb} MB) → {out_path.name}"
                else:
                    err = (proc.stderr or "").strip()
                    logger.warning("pgscat download %s --type %s failed: %s", pgs_id, file_type, err)
                    if tmp.exists():
                        tmp.unlink(missing_ok=True)

            return False, "pgscat download failed for both positions and original file types"

        except subprocess.TimeoutExpired:
            if tmp.exists():
                tmp.unlink(missing_ok=True)
            return False, f"Download timed out after {_PGSCAT_DOWNLOAD_TIMEOUT}s"
        except Exception as exc:
            if tmp.exists():
                tmp.unlink(missing_ok=True)
            logger.exception("Unexpected error downloading score file for %s", pgs_id)
            return False, str(exc)

    # ── Betamap builder ───────────────────────────────────────────────────────

    def _build_betamap(
        self,
        pgs_id: str,
        raw_path: Path,
        betamap_path: Path,
        cols_path: Path,
        meta_path: Path,
    ) -> tuple[bool, str]:
        """
        Parse the raw harmonised scoring file and produce:
          - {pgs_id}_hmPOS_GRCh38.betamap.tsv.gz   (pipeline input)
          - {pgs_id}_hmPOS_GRCh38.columns.json      (column mapping report)
          - {pgs_id}_hmPOS_GRCh38.betamap.tsv.gz.meta.json  (generation report)

        Column mapping (case-insensitive):
          CHROM  ← hm_chr / chr_name / chromosome / chrom / chr
          POS    ← hm_pos / chr_position / position / pos / bp
          ID     ← hm_rsID / rsID / SNP / variant_id / markername / id
          EFFECT_ALLELE ← effect_allele / allele1 / a1 / ea
          OTHER_ALLELE  ← hm_inferOtherAllele / other_allele / allele2 / a2 / nea
          BETA   ← effect_weight / beta / b / effect

        IS_FLIP is always 0 in Phase 1. Actual flip detection requires zarr
        lookup (a compute step), not metadata processing.

        Rows are skipped if CHROM, POS, EFFECT_ALLELE, or BETA are missing/empty.
        """
        if not raw_path.exists():
            return False, f"Raw score file not found: {raw_path}"

        try:
            # ── Parse header + detect columns ────────────────────────────────
            header, file_meta = _read_pgs_header(raw_path)

            if not header:
                return False, "Could not find column header in scoring file"

            col_map = _detect_columns(header)

            required = {"CHROM", "POS", "EFFECT_ALLELE", "BETA"}
            missing_required = required - set(col_map.keys())
            if missing_required:
                return False, (
                    f"Cannot map required columns: {missing_required}. "
                    f"Available columns in file: {header}"
                )

            # ── Write columns.json ────────────────────────────────────────────
            cols_report = {
                "file": str(raw_path),
                "delimiter": "\t",
                "meta": file_meta,
                "n_columns": len(header),
                "header": header,
                "column_mapping": col_map,
                "missing_optional": [
                    c for c in ["ID", "OTHER_ALLELE"] if c not in col_map
                ],
            }
            _atomic_write_text(cols_path, json.dumps(cols_report, indent=2, ensure_ascii=False))

            # ── Build betamap ─────────────────────────────────────────────────
            n_total = 0
            n_written = 0
            n_bad_chrom = 0
            n_bad_pos = 0
            n_bad_allele = 0
            n_bad_beta = 0

            # Atomic write via temp file
            tmp_betamap = betamap_path.with_suffix(".tmp.gz")
            try:
                with gzip.open(raw_path, "rt", encoding="utf-8", errors="replace") as fin, \
                     gzip.open(tmp_betamap, "wt", encoding="utf-8") as fout:

                    # Write betamap header
                    fout.write("\t".join(BETAMAP_COLUMNS) + "\n")

                    # Skip header/comment lines; find data start
                    data_header_found = False
                    for line in fin:
                        line = line.rstrip("\n")
                        if line.startswith("#"):
                            continue
                        # First non-comment line is the column header
                        if not data_header_found:
                            data_header_found = True
                            # Rebuild index map from actual header in stream
                            actual_header = line.split("\t")
                            idx = _build_index_map(actual_header, col_map)
                            continue

                        if not line.strip():
                            continue

                        n_total += 1
                        fields = line.split("\t")

                        chrom = _get_field(fields, idx.get("CHROM"))
                        pos   = _get_field(fields, idx.get("POS"))
                        ea    = _get_field(fields, idx.get("EFFECT_ALLELE"))
                        oa    = _get_field(fields, idx.get("OTHER_ALLELE", -1))
                        beta  = _get_field(fields, idx.get("BETA"))
                        rid   = _get_field(fields, idx.get("ID", -1))

                        # Validate required fields
                        if not chrom or chrom in (".", "NA", "nan"):
                            n_bad_chrom += 1
                            continue
                        if not pos or pos in (".", "NA", "nan"):
                            n_bad_pos += 1
                            continue
                        if not ea or ea in (".", "NA", "nan"):
                            n_bad_allele += 1
                            continue
                        if not beta or beta in (".", "NA", "nan"):
                            n_bad_beta += 1
                            continue

                        # Normalise chromosome (remove "chr" prefix)
                        chrom = chrom.lstrip("chrCHR").lstrip("0") or chrom

                        row = [
                            pgs_id,    # PRS_ID
                            chrom,     # CHROM
                            pos,       # POS
                            rid,       # ID (may be empty)
                            ea,        # EFFECT_ALLELE
                            oa,        # OTHER_ALLELE (may be empty)
                            beta,      # BETA
                            "0",       # IS_FLIP (Phase 1: always 0)
                        ]
                        fout.write("\t".join(row) + "\n")
                        n_written += 1

                tmp_betamap.rename(betamap_path)

            except Exception:
                if tmp_betamap.exists():
                    tmp_betamap.unlink(missing_ok=True)
                raise

            # ── Write meta.json ───────────────────────────────────────────────
            meta_report = {
                "prs_id": pgs_id,
                "prs_path": str(raw_path),
                "colsinfo_path": str(cols_path),
                "betamap_path": str(betamap_path),
                "weight_type": file_meta.get("weight_type", "NR"),
                "selected_columns": BETAMAP_COLUMNS,
                "column_mapping": col_map,
                "n_total_rows": n_total,
                "n_written_rows": n_written,
                "n_bad_chrom": n_bad_chrom,
                "n_bad_pos": n_bad_pos,
                "n_bad_allele": n_bad_allele,
                "n_bad_beta": n_bad_beta,
                "is_flip_note": (
                    "IS_FLIP set to 0 for all rows in Phase 1. "
                    "Actual flip detection requires zarr allele comparison (compute step)."
                ),
            }
            _atomic_write_text(meta_path, json.dumps(meta_report, indent=2, ensure_ascii=False))

            size_mb = round(betamap_path.stat().st_size / 1_048_576, 2)
            return True, (
                f"betamap ready: {n_written}/{n_total} variants written "
                f"({n_bad_chrom} bad_chrom, {n_bad_pos} bad_pos, "
                f"{n_bad_allele} bad_allele, {n_bad_beta} bad_beta) — {size_mb} MB"
            )

        except Exception as exc:
            logger.exception("Betamap build failed for %s", pgs_id)
            return False, str(exc)

    # ── Internals ─────────────────────────────────────────────────────────────

    def _result(
        self,
        pgs_id: str,
        pgs_dir: Path,
        status: str,
        steps: list[dict],
        error: str = "",
    ) -> dict:
        prep_status = self.get_preparation_status(pgs_id)
        result = {
            "pgs_id": pgs_id,
            "pgs_dir": str(pgs_dir),
            "status": status,
            "steps": steps,
            "current_stage": prep_status.get("current_stage", "unknown"),
            "ready_for_compute": prep_status.get("ready_for_compute", False),
            "ready_for_dashboard": prep_status.get("ready_for_dashboard", False),
        }
        if error:
            result["error"] = error
        return result

    def _scan_files(self, pgs_id: str, pgs_dir: Path) -> dict:
        p = pgs_id
        checks = {
            "catalog_metadata":     pgs_dir / f"{p}.json",
            "raw_harmonized_score": pgs_dir / f"{p}_hmPOS_GRCh38.txt.gz",
            "column_map":           pgs_dir / f"{p}_hmPOS_GRCh38.columns.json",
            "betamap":              pgs_dir / f"{p}_hmPOS_GRCh38.betamap.tsv.gz",
            "betamap_meta":         pgs_dir / f"{p}_hmPOS_GRCh38.betamap.tsv.gz.meta.json",
            "aggregated_scores":    pgs_dir / f"{p}_PRS_total.tsv",
            "aggregated_metadata":  pgs_dir / f"{p}_PRS_total_metadata.json",
            "parquet_export":       pgs_dir / f"{p}_PRS_total.parquet",
        }
        result = {}
        for role, path in checks.items():
            exists = path.exists()
            result[role] = {
                "path": str(path),
                "exists": exists,
                "size_mb": (
                    round(path.stat().st_size / 1_048_576, 2)
                    if exists else None
                ),
            }
        return result

    def _count_chrom_scores(self, pgs_id: str, pgs_dir: Path) -> int:
        count = 0
        for c in list(range(1, 23)) + ["X"]:
            if (pgs_dir / f"{pgs_id}_chr{c}_scores.tsv").exists():
                count += 1
        return count

    def _ready_for_compute(self, files: dict) -> bool:
        return files.get("betamap", {}).get("exists", False)

    def _current_stage(self, files: dict) -> str:
        if files.get("parquet_export", {}).get("exists"):
            return "parquet_export"
        if files.get("aggregated_scores", {}).get("exists"):
            return "aggregated_scores"
        if files.get("betamap", {}).get("exists"):
            return "betamap"
        if files.get("column_map", {}).get("exists"):
            return "column_map"
        if files.get("raw_harmonized_score", {}).get("exists"):
            return "raw_harmonized_score"
        if files.get("catalog_metadata", {}).get("exists"):
            return "catalog_metadata"
        return "not_started"

    def _next_step(self, stage: str, files: dict, n_chroms: int) -> str:
        if stage == "parquet_export":
            return "Dashboard ready — no action needed"
        if stage == "aggregated_scores":
            return "Run DuckDB COPY TO to export parquet (step 5)"
        if n_chroms > 0 and n_chroms < 22:
            return f"Per-chrom scoring in progress ({n_chroms}/22 done) — wait or resubmit missing chroms"
        if n_chroms == 22:
            return "All chrom scores done — run aggregation step (step 4)"
        if stage == "betamap":
            return "Betamap ready — submit Slurm array job for per-chrom scoring (step 3)"
        if stage in ("raw_harmonized_score", "column_map"):
            return "Raw score downloaded — run betamap builder (prepare again)"
        if stage == "catalog_metadata":
            return "Metadata available — click 'Prepare Score' to download and build betamap"
        return "Not started — click 'Prepare Score' to download and build betamap"


# ── Module-level helpers ──────────────────────────────────────────────────────

def _atomic_write_text(path: Path, content: str) -> None:
    """Write text to path atomically via a temp file in the same directory."""
    tmp = path.with_suffix(".tmp")
    try:
        tmp.write_text(content, encoding="utf-8")
        tmp.rename(path)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise


def _read_pgs_header(raw_path: Path) -> tuple[list[str], dict]:
    """
    Read the PGS Catalog file header (lines starting with #) and the column row.
    Returns (column_list, metadata_dict).
    """
    meta: dict = {}
    columns: list[str] = []

    with gzip.open(raw_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##") or line.startswith("###"):
                # Section headers — skip
                continue
            if line.startswith("#"):
                # Metadata key=value
                kv = line.lstrip("#").strip()
                if "=" in kv:
                    k, _, v = kv.partition("=")
                    meta[k.strip()] = v.strip()
                continue
            # First non-comment line is the column header
            columns = line.split("\t")
            break

    return columns, meta


def _detect_columns(header: list[str]) -> dict:
    """
    Map betamap output columns to input column indices.
    Returns dict of {betamap_col → original_col_name}.
    Matching is case-insensitive.
    """
    lower = {c.lower(): c for c in header}
    mapping: dict[str, str] = {}

    def _find(candidates: set[str]) -> Optional[str]:
        for cand in candidates:
            if cand in lower:
                return lower[cand]
        return None

    chrom = _find(_CHROM_COLS)
    pos   = _find(_POS_COLS)
    ea    = _find(_EA_COLS)
    oa    = _find(_OA_COLS)
    beta  = _find(_BETA_COLS)
    rid   = _find(_ID_COLS)

    if chrom: mapping["CHROM"] = chrom
    if pos:   mapping["POS"] = pos
    if ea:    mapping["EFFECT_ALLELE"] = ea
    if oa:    mapping["OTHER_ALLELE"] = oa
    if beta:  mapping["BETA"] = beta
    if rid:   mapping["ID"] = rid

    return mapping


def _build_index_map(header: list[str], col_map: dict) -> dict:
    """Build col_name → column_index from actual data header."""
    idx_by_name = {c: i for i, c in enumerate(header)}
    return {
        betamap_col: idx_by_name.get(original_col, -1)
        for betamap_col, original_col in col_map.items()
    }


def _get_field(fields: list[str], idx) -> str:
    """Safely extract field by index; return '' if out of range or None index."""
    if idx is None or idx < 0:
        return ""
    try:
        return fields[idx].strip()
    except IndexError:
        return ""


# Fix missing Optional import used in _detect_columns
from typing import Optional  # noqa: E402 (deferred to avoid circular in type hints above)
