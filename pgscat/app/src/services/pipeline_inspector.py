"""
PipelineInspector service
──────────────────────────
Generates structured, read-only processing plans for three distinct flows:

  local_existing         PGS already computed locally → describe what exists
  catalog_remote_prepare PGS in PGS Catalog but not local → full prep + compute plan
  custom_prs             User-defined score → format + compute plan

IMPORTANT: This service NEVER executes any job or modifies any file.
It only returns informational plans and file-layout descriptions.

─── Local file layout contract ──────────────────────────────────────────────
Expected layout under {APP_DATA_DIR}/{PGS_ID}/:

  {PGS_ID}.json                               PGS Catalog metadata (downloaded once)
  {PGS_ID}_hmPOS_GRCh38.txt.gz               Raw harmonized scoring file from PGS Catalog
  {PGS_ID}_hmPOS_GRCh38.columns.json         Column mapping (auto-detected by betamap builder)
  {PGS_ID}_hmPOS_GRCh38.betamap.tsv.gz       Pipeline input — betamap with columns:
                                               PRS_ID, CHROM, POS, ID,
                                               EFFECT_ALLELE, OTHER_ALLELE, BETA, IS_FLIP
  {PGS_ID}_hmPOS_GRCh38.betamap.tsv.gz.meta.json  Betamap generation report
  {PGS_ID}_chr{N}_scores.tsv                 Per-chrom PRS (cols: sample_id, PRS)
  {PGS_ID}_chr{N}_metadata.json              Per-chrom run metadata
  {PGS_ID}_PRS_total.tsv                     Aggregated (sample_id, PRS_total, PRS_chr1..22)
  {PGS_ID}_PRS_total_metadata.json           Aggregation report
  {PGS_ID}_PRS_total.parquet                 Parquet export (enables web dashboard)

─── Scripts source ──────────────────────────────────────────────────────────
Pipeline scripts directory: configured via APP_SCRIPTS_DIR (default: /pipeline_scripts)
Key scripts:
  compute_prs.py             CPU/DuckDB additive PRS scoring
  compute_prs_spark_gpu.py   Spark variant matching + CuPy GPU scoring
  prs_gpu_compute.sbatch     Slurm array job (22 chromosomes) via GPU
  run_prs_spark_gpu.sbatch   Slurm single-node Spark+GPU job
  sbrc_01_download_format.sbatch  Download + format scoring file
"""
import json
import logging
from pathlib import Path
from typing import Optional

from config import Config
from services.local_catalog import validate_pgs_id

logger = logging.getLogger(__name__)

# Strict allowlist for script content reads (prevents path traversal)
_ALLOWED_SCRIPTS = frozenset([
    "compute_prs.py",
    "compute_prs_spark_gpu.py",
    "prs_gpu_compute.sbatch",
    "run_prs_spark_gpu.sbatch",
    "prs_end2end.sbatch",
    "prs_end2end_v3.sbatch",
    "sbrc_01_download_format.sbatch",
])

# ── Zarr / index constants ────────────────────────────────────────────────────
_ZARR_BASE = "{GENOTYPES_DIR}/chr{N}.zarr"
_ZARR_INDEX = "{GENOTYPES_DIR}/zarr_index.duckdb"

# ── Betamap column contract ───────────────────────────────────────────────────
BETAMAP_COLUMNS = ["PRS_ID", "CHROM", "POS", "ID", "EFFECT_ALLELE", "OTHER_ALLELE", "BETA", "IS_FLIP"]


class PipelineInspector:
    def __init__(self, cfg: Config) -> None:
        self.cfg = cfg

    # ── Public API ───────────────────────────────────────────────────────────

    def get_pipeline_plan(
        self,
        pgs_id: str,
        flow_mode: str = "auto",
        local_exists: bool = False,
        remote_info: Optional[dict] = None,
    ) -> dict:
        """
        Return a structured processing plan for pgs_id.

        flow_mode:
          "auto"                  — detect from local_exists / remote_info
          "local_existing"        — PGS computed locally
          "catalog_remote_prepare" — PGS in catalog but not local
          "custom_prs"            — user-defined score (pgs_id may be arbitrary)

        No jobs are submitted or files modified.
        """
        if not validate_pgs_id(pgs_id):
            return {"error": "Invalid PGS ID format", "pgs_id": pgs_id}

        # Auto-detect flow
        if flow_mode == "auto":
            if local_exists:
                flow_mode = "local_existing"
            elif remote_info:
                flow_mode = "catalog_remote_prepare"
            else:
                flow_mode = "catalog_remote_prepare"

        base = {
            "pgs_id": pgs_id,
            "flow_mode": flow_mode,
            "pgs_scores_dir": str(self.cfg.scores_dir / pgs_id),
            "zarr_base": _ZARR_BASE,
            "duckdb_index": _ZARR_INDEX,
            "scripts_dir": str(self.cfg.SCRIPTS_DIR),
            "scripts_available": self.list_scripts(),
            "local_file_layout": self._local_layout(pgs_id),
            "note": (
                "This plan is informational only. "
                "No jobs have been submitted and no files have been modified."
            ),
        }

        if flow_mode == "local_existing":
            base["steps"] = self._steps_local(pgs_id)
            base["description"] = "PGS is already computed locally. Pipeline steps describe what was run."
        elif flow_mode == "catalog_remote_prepare":
            base["steps"] = self._steps_remote_prepare(pgs_id, remote_info or {})
            base["description"] = (
                "PGS exists in the PGS Catalog but not locally. "
                "Steps describe download, format preparation, and compute plan."
            )
        elif flow_mode == "custom_prs":
            base["steps"] = self._steps_custom_prs(pgs_id)
            base["description"] = (
                "Custom PRS: bring your own scoring file and run the same pipeline."
            )
        else:
            base["error"] = f"Unknown flow_mode: {flow_mode!r}"

        return base

    def get_local_layout_status(self, pgs_id: str) -> dict:
        """
        Check which local layout files exist for pgs_id.
        Returns dict with file → {path, exists, required}.
        """
        if not validate_pgs_id(pgs_id):
            return {"error": "Invalid PGS ID"}
        d = self.cfg.scores_dir / pgs_id
        layout = {}
        for item in self._local_layout(pgs_id)["files"]:
            p = d / item["filename"]
            layout[item["role"]] = {
                "path": str(p),
                "filename": item["filename"],
                "exists": p.exists(),
                "required": item["required"],
                "description": item["description"],
            }
        return layout

    def list_scripts(self) -> list[dict]:
        """List known pipeline scripts with filesystem presence check."""
        result = []
        for name in sorted(_ALLOWED_SCRIPTS):
            path = self.cfg.SCRIPTS_DIR / name
            exists = path.exists()
            result.append({
                "name": name,
                "path": str(path),
                "exists": exists,
                "size_kb": round(path.stat().st_size / 1024, 1) if exists else None,
            })
        return result

    def read_script(self, script_name: str) -> Optional[str]:
        """
        Return the text content of a known pipeline script for display.
        Returns None if script is not in the allowlist or doesn't exist.
        """
        if script_name not in _ALLOWED_SCRIPTS:
            logger.warning("Attempted to read non-allowed script: %r", script_name)
            return None
        path = self.cfg.SCRIPTS_DIR / script_name
        if not path.exists():
            return None
        try:
            return path.read_text(encoding="utf-8", errors="replace")
        except Exception as exc:
            logger.error("Cannot read script %s: %s", path, exc)
            return None

    # ── Local layout contract ─────────────────────────────────────────────────

    def _local_layout(self, pgs_id: str) -> dict:
        """
        Return the complete expected file layout for a locally prepared PGS.
        This is the definitive file layout contract for a locally prepared PGS.
        See docs/ARCHITECTURE.md for the full column and file specification.
        """
        p = pgs_id
        return {
            "base_dir": f"{{SCORES_DIR}}/{p}/",
            "betamap_columns": BETAMAP_COLUMNS,
            "files": [
                {
                    "filename": f"{p}.json",
                    "role": "catalog_metadata",
                    "description": "PGS Catalog metadata (JSON from API or pgscat download)",
                    "required": True,
                    "produced_by": "pgscat / manual download",
                },
                {
                    "filename": f"{p}_hmPOS_GRCh38.txt.gz",
                    "role": "raw_harmonized_score",
                    "description": "Raw GRCh38-harmonised scoring file from PGS Catalog FTP",
                    "required": True,
                    "produced_by": "pgscat download / sbrc_01_download_format.sbatch",
                },
                {
                    "filename": f"{p}_hmPOS_GRCh38.columns.json",
                    "role": "column_map",
                    "description": "Auto-detected column mapping from raw file to betamap schema",
                    "required": True,
                    "produced_by": "betamap builder (format step)",
                },
                {
                    "filename": f"{p}_hmPOS_GRCh38.betamap.tsv.gz",
                    "role": "betamap",
                    "description": (
                        f"Pipeline input — tab-separated, columns: {', '.join(BETAMAP_COLUMNS)}. "
                        "IS_FLIP=1 means effect is on REF allele; dosage = 2 − call_DS."
                    ),
                    "required": True,
                    "produced_by": "betamap builder (format step)",
                },
                {
                    "filename": f"{p}_hmPOS_GRCh38.betamap.tsv.gz.meta.json",
                    "role": "betamap_meta",
                    "description": "Betamap generation report (n_total_rows, n_written_rows, n_skipped_*)",
                    "required": False,
                    "produced_by": "betamap builder (format step)",
                },
                {
                    "filename": f"{p}_chr{{N}}_scores.tsv",
                    "role": "per_chrom_scores",
                    "description": "Per-chromosome PRS results (cols: sample_id, PRS). One file per chrom.",
                    "required": True,
                    "produced_by": "compute_prs_spark_gpu.py / prs_gpu_compute.sbatch",
                },
                {
                    "filename": f"{p}_chr{{N}}_metadata.json",
                    "role": "per_chrom_metadata",
                    "description": "Per-chromosome run metadata (n_samples, n_variants_matched, timing, etc.)",
                    "required": False,
                    "produced_by": "compute_prs_spark_gpu.py",
                },
                {
                    "filename": f"{p}_PRS_total.tsv",
                    "role": "aggregated_scores",
                    "description": "Aggregated PRS (cols: sample_id, PRS_total, PRS_chr1..PRS_chr22). 140k+ rows.",
                    "required": True,
                    "produced_by": "aggregation step (post per-chrom scoring)",
                },
                {
                    "filename": f"{p}_PRS_total_metadata.json",
                    "role": "aggregated_metadata",
                    "description": "Aggregation report (n_samples, n_chr_scores, per_chromosome stats)",
                    "required": False,
                    "produced_by": "aggregation step",
                },
                {
                    "filename": f"{p}_PRS_total.parquet",
                    "role": "parquet_export",
                    "description": "Parquet export of PRS_total.tsv (ZSTD). Required for web dashboard.",
                    "required": False,
                    "produced_by": "DuckDB COPY TO (post-aggregation step)",
                },
            ],
        }

    # ── Step builders ─────────────────────────────────────────────────────────

    def _steps_local(self, pgs_id: str) -> list[dict]:
        """Plan for a PGS that already exists locally."""
        sd = str(self.cfg.scores_dir)
        p = pgs_id
        return [
            {
                "step": 1,
                "status": "done",
                "name": "Score downloaded & betamap prepared",
                "description": (
                    f"Files already present: {p}_hmPOS_GRCh38.txt.gz, "
                    f"{p}_hmPOS_GRCh38.betamap.tsv.gz"
                ),
                "output_pattern": f"{p}_hmPOS_GRCh38.betamap.tsv.gz",
                "tool": "sbrc_01_download_format.sbatch",
                "notes": f"Betamap columns: {', '.join(BETAMAP_COLUMNS)}",
            },
            {
                "step": 2,
                "status": "done",
                "name": "Per-chromosome PRS scoring",
                "description": (
                    "Spark + CuPy GPU scoring completed for all chromosomes. "
                    "Zarr dosage data matched to betamap by position + allele."
                ),
                "output_pattern": f"{p}_chr{{N}}_scores.tsv  +  {p}_chr{{N}}_metadata.json",
                "tool": "compute_prs_spark_gpu.py / prs_gpu_compute.sbatch",
                "notes": "IS_FLIP variants: dosage = 2 − call_DS",
            },
            {
                "step": 3,
                "status": "done",
                "name": "PRS aggregation",
                "description": f"All per-chromosome scores summed into {p}_PRS_total.tsv.",
                "output_pattern": f"{p}_PRS_total.tsv  +  {p}_PRS_total_metadata.json",
                "tool": "aggregation step",
                "notes": "PRS_total[s] = Σ_chr PRS_chr[s]",
            },
            {
                "step": 4,
                "status": "done" ,
                "name": "Parquet export",
                "description": f"TSV converted to {p}_PRS_total.parquet for web dashboard.",
                "output_pattern": f"{p}_PRS_total.parquet",
                "tool": "DuckDB COPY TO",
                "command_example": (
                    f'duckdb -c "COPY (SELECT * FROM read_csv_auto(\'{sd}/{p}/{p}_PRS_total.tsv\', '
                    f"delim=chr(9), header=true)) "
                    f"TO '{sd}/{p}/{p}_PRS_total.parquet' "
                    f'(FORMAT PARQUET, COMPRESSION ZSTD)"'
                ),
                "notes": "ZSTD compression gives ~3× size reduction.",
            },
        ]

    def _steps_remote_prepare(self, pgs_id: str, remote_info: dict) -> list[dict]:
        """Full plan: download from PGS Catalog → prepare → score → aggregate → parquet."""
        sd = str(self.cfg.scores_dir)
        p = pgs_id
        scripts = str(self.cfg.SCRIPTS_DIR)
        ftp_url = remote_info.get("ftp_scoring_file", "")
        is_harmonized = remote_info.get("is_harmonized", False)
        n_variants = remote_info.get("variants_number")

        harmonized_note = (
            "Score is GRCh38-harmonised — use HmPOS columns (hm_chr, hm_pos, hm_rsID)."
            if is_harmonized else
            "Score may NOT be GRCh38-harmonised. Verify build before running pipeline. "
            "Consider using pgscat to get the harmonized version."
        )

        return [
            {
                "step": 1,
                "status": "pending",
                "name": "Download scoring file from PGS Catalog",
                "description": (
                    f"Download the GRCh38-harmonised scoring file for {p} from PGS Catalog. "
                    f"Expected output: {p}_hmPOS_GRCh38.txt.gz"
                ),
                "output_pattern": f"{p}_hmPOS_GRCh38.txt.gz",
                "tool": "pgscat CLI",
                "command_example": (
                    f"mkdir -p {sd}/{p}\n"
                    f"pgscat download {p} --build GRCh38 --out {sd}/{p}/{p}_hmPOS_GRCh38.txt.gz"
                    + (f"\n# Direct FTP URL: {ftp_url}" if ftp_url else "")
                ),
                "submit_example": (
                    f"# sbatch wrapper (optional):\n"
                    f"sbatch --export=PGS_ID={p},OUT_DIR={sd}/{p} "
                    f"{scripts}/sbrc_01_download_format.sbatch"
                ),
                "notes": (
                    f"{harmonized_note} "
                    f"n_variants in catalog: {n_variants if n_variants else 'unknown'}"
                ),
                "automation": "automatic",
            },
            {
                "step": 2,
                "status": "pending",
                "name": "Build betamap (format score for pipeline)",
                "description": (
                    f"Map columns from raw scoring file to betamap schema: "
                    f"{', '.join(BETAMAP_COLUMNS)}. "
                    f"Auto-detect CHROM/POS/EFFECT_ALLELE/OTHER_ALLELE/BETA columns. "
                    f"Compute IS_FLIP flag (1 = effect on REF allele). "
                    f"Output: {p}_hmPOS_GRCh38.betamap.tsv.gz + .columns.json + .meta.json"
                ),
                "output_pattern": (
                    f"{p}_hmPOS_GRCh38.betamap.tsv.gz  +  "
                    f"{p}_hmPOS_GRCh38.columns.json  +  "
                    f"{p}_hmPOS_GRCh38.betamap.tsv.gz.meta.json"
                ),
                "tool": "betamap builder (sbrc_01_download_format.sbatch or custom script)",
                "command_example": (
                    f"# Included in sbrc_01_download_format.sbatch after download step\n"
                    f"# Input:  {sd}/{p}/{p}_hmPOS_GRCh38.txt.gz\n"
                    f"# Output: {sd}/{p}/{p}_hmPOS_GRCh38.betamap.tsv.gz"
                ),
                "notes": (
                    "IS_FLIP=1 when effect allele matches REF in zarr (dosage = 2 − call_DS). "
                    "Rows with unresolvable alleles are excluded (recorded in .meta.json)."
                ),
                "automation": "automatic",
            },
            {
                "step": 3,
                "status": "pending",
                "name": "Per-chromosome PRS scoring (Spark + CuPy GPU)",
                "description": (
                    f"For each chromosome 1–22: match betamap variants to zarr dosage data "
                    f"by position + allele (strand-flip aware). "
                    f"Compute additive PRS: Σ_i(BETA_i × dosage_i) per sample (~140k samples). "
                    f"Zarr: {_ZARR_BASE}  |  Index: {_ZARR_INDEX}"
                ),
                "output_pattern": (
                    f"{p}_chr{{N}}_scores.tsv  +  {p}_chr{{N}}_metadata.json"
                ),
                "tool": "compute_prs_spark_gpu.py",
                "slurm_script": "prs_gpu_compute.sbatch",
                "slurm_params": {
                    "partition": "gpu",
                    "account": "researchers",
                    "qos": "vip",
                    "nodes": 1,
                    "cpus_per_task": 16,
                    "mem": "120G",
                    "gres": "gpu:1",
                    "time": "06:00:00",
                },
                "submit_example": (
                    f"# Array job — all 22 autosomes in parallel:\n"
                    f"sbatch --array=1-22 --export=PGS_ID={p},SCORES_DIR={sd} "
                    f"{scripts}/prs_gpu_compute.sbatch\n\n"
                    f"# Or single chromosome:\n"
                    f"sbatch --export=PGS_ID={p},CHROM=1,SCORES_DIR={sd} "
                    f"{scripts}/run_prs_spark_gpu.sbatch"
                ),
                "notes": (
                    "gpu_batch_matches=4000 avoids OOM on RTX 3060 (12 GB VRAM). "
                    "CPU fallback: compute_prs.py with DuckDB index. "
                    "Chromosome X excluded from standard array (add manually if needed)."
                ),
                "automation": "informational — requires Slurm / manual submit",
            },
            {
                "step": 4,
                "status": "pending",
                "name": "Aggregate total PRS",
                "description": (
                    f"Sum all per-chromosome PRS values per sample. "
                    f"Output columns: sample_id, PRS_total, PRS_chr1, …, PRS_chr22."
                ),
                "output_pattern": (
                    f"{p}_PRS_total.tsv  +  {p}_PRS_total_metadata.json"
                ),
                "tool": "aggregation step (post per-chrom scoring)",
                "command_example": (
                    f"# After all per-chrom scores are complete:\n"
                    f"python3 {scripts}/compute_prs.py --aggregate "
                    f"--pgs-id {p} --scores-dir {sd}/{p}"
                ),
                "notes": "PRS_total[s] = Σ_chr PRS_chr[s] across all 22 autosomes.",
                "automation": "informational — depends on step 3 completion",
            },
            {
                "step": 5,
                "status": "pending",
                "name": "Export to Parquet (enables web dashboard)",
                "description": (
                    f"Convert {p}_PRS_total.tsv to {p}_PRS_total.parquet using DuckDB. "
                    "Once this file exists the web dashboard becomes available."
                ),
                "output_pattern": f"{p}_PRS_total.parquet",
                "tool": "DuckDB COPY TO",
                "command_example": (
                    f'duckdb -c "COPY (SELECT * FROM read_csv_auto(\'{sd}/{p}/{p}_PRS_total.tsv\', '
                    f"delim=chr(9), header=true)) "
                    f"TO '{sd}/{p}/{p}_PRS_total.parquet' "
                    f'(FORMAT PARQUET, COMPRESSION ZSTD)"'
                ),
                "notes": "ZSTD compression gives ~3× size reduction. DuckDB reads both formats.",
                "automation": "automatic — no Slurm needed",
            },
        ]

    def _steps_custom_prs(self, pgs_id: str) -> list[dict]:
        """Plan for a user-defined custom PRS score."""
        sd = str(self.cfg.scores_dir)
        p = pgs_id
        scripts = str(self.cfg.SCRIPTS_DIR)
        return [
            {
                "step": 1,
                "status": "pending",
                "name": "Prepare scoring file",
                "description": (
                    "Provide your scoring file in GRCh38 coordinates. "
                    "Required columns: CHROM, POS, EFFECT_ALLELE, OTHER_ALLELE, BETA. "
                    "Optionally: rsID (ID column), allele frequency."
                ),
                "output_pattern": f"{p}_hmPOS_GRCh38.txt.gz",
                "tool": "User-provided",
                "notes": (
                    "If your file is GRCh37/hg19, lift over to GRCh38 first. "
                    "Scripts available: sbrc_fix3_liftover_hg19_to_hg38.sbatch"
                ),
                "automation": "manual — user provides file",
            },
            {
                "step": 2,
                "status": "pending",
                "name": "Build betamap",
                "description": (
                    f"Map your columns to betamap schema: {', '.join(BETAMAP_COLUMNS)}. "
                    "Compute IS_FLIP flag. Output must match the pipeline input contract."
                ),
                "output_pattern": (
                    f"{p}_hmPOS_GRCh38.betamap.tsv.gz  +  "
                    f"{p}_hmPOS_GRCh38.columns.json  +  "
                    f"{p}_hmPOS_GRCh38.betamap.tsv.gz.meta.json"
                ),
                "tool": "betamap builder",
                "command_example": (
                    f"# Adapt sbrc_01_download_format.sbatch to your input file path\n"
                    f"# Input:  {sd}/{p}/<your_score_file>.txt.gz\n"
                    f"# Output: {sd}/{p}/{p}_hmPOS_GRCh38.betamap.tsv.gz"
                ),
                "notes": "IS_FLIP=1 when effect allele = REF in zarr (dosage = 2 − call_DS).",
                "automation": "semi-automatic — requires column mapping review",
            },
            {
                "step": 3,
                "status": "pending",
                "name": "Per-chromosome PRS scoring",
                "description": (
                    "Same as catalog flow. Run per-chrom scoring with your betamap. "
                    f"Zarr: {_ZARR_BASE}"
                ),
                "output_pattern": f"{p}_chr{{N}}_scores.tsv  +  {p}_chr{{N}}_metadata.json",
                "tool": "compute_prs_spark_gpu.py",
                "submit_example": (
                    f"sbatch --array=1-22 --export=PGS_ID={p},SCORES_DIR={sd} "
                    f"{scripts}/prs_gpu_compute.sbatch"
                ),
                "notes": "CPU fallback: compute_prs.py with DuckDB index.",
                "automation": "informational — requires Slurm / manual submit",
            },
            {
                "step": 4,
                "status": "pending",
                "name": "Aggregate + Parquet export",
                "description": "Same as catalog flow. Sum per-chrom scores, export to parquet.",
                "output_pattern": f"{p}_PRS_total.tsv  +  {p}_PRS_total.parquet",
                "tool": "aggregation step + DuckDB",
                "automation": "automatic — no Slurm needed",
            },
        ]
