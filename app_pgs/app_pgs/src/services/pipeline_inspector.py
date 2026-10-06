"""
PipelineInspector service
──────────────────────────
Provides a read-only view of the existing HPC pipeline scripts and
generates a structured processing plan for any PGS ID.

IMPORTANT: This service NEVER executes any job.
It only reads script files and returns informational plans.

Pipeline derived from actual inspection of:
  compute_prs.py             – CPU/DuckDB additive PRS scoring
  compute_prs_spark_gpu.py   – Spark variant matching + CuPy GPU scoring
  prs_gpu_compute.sbatch     – Slurm array job (22 chromosomes) via GPU
  run_prs_spark_gpu.sbatch   – Slurm single-node Spark+GPU job
  prs_end2end.sbatch / v3    – Clumping + scoring for GWAS-derived PRS

Key pipeline facts (verified from code):
  Input weights:   {PGS_ID}_hmPOS_GRCh38.betamap.tsv.gz
  Columns:         PRS_ID, CHROM, POS, ID, EFFECT_ALLELE, OTHER_ALLELE, BETA, IS_FLIP
  Genotype data:   /mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files/chr{N}.zarr
  Per-chrom out:   {PGS_ID}_chr{N}_scores.tsv  (cols: sample_id, PRS)
  Aggregated out:  {PGS_ID}_PRS_total.tsv  (cols: sample_id, PRS_total, PRS_chr1..chr22)
  Parquet export:  {PGS_ID}_PRS_total.parquet  (recommended new step, not yet in scripts)
"""
import json
import logging
from pathlib import Path
from typing import Optional

from config import Config
from services.local_catalog import validate_pgs_id

logger = logging.getLogger(__name__)

# Scripts whose content may be read via read_script() – strict allowlist
_ALLOWED_SCRIPTS = frozenset([
    "prs_gpu_compute.sbatch",
    "run_prs_spark_gpu.sbatch",
    "compute_prs.py",
    "compute_prs_spark_gpu.py",
    "prs_end2end.sbatch",
    "prs_end2end_v3.sbatch",
])

# Static plan template (substituted with actual pgs_id at call time)
_PLAN_TEMPLATE: list[dict] = [
    {
        "step": 1,
        "name": "Download & prepare betamap",
        "description": (
            "Download the GRCh38-harmonised scoring file from PGS Catalog "
            "using `pgscat download {pgs_id} --build GRCh38`. "
            "Preprocess to betamap TSV: columns PRS_ID, CHROM, POS, ID, "
            "EFFECT_ALLELE, OTHER_ALLELE, BETA, IS_FLIP."
        ),
        "output_pattern": "{pgs_id}_hmPOS_GRCh38.betamap.tsv.gz",
        "tool": "pgscat CLI + custom preprocessing",
        "slurm_script": None,
        "notes": "IS_FLIP=1 means effect is on REF allele; dosage = 2 − call_DS.",
    },
    {
        "step": 2,
        "name": "Per-chromosome PRS scoring (Spark + CuPy GPU)",
        "description": (
            "For each chromosome 1–22 (+X): match betamap variants to zarr "
            "dosage data by position+allele (strand-flip aware). "
            "Compute additive PRS: Σ_i(BETA_i × dosage_i) per sample. "
            "Recommended: Spark broadcast join for matching + CuPy for GPU dot products."
        ),
        "output_pattern": "{pgs_id}_chr{N}_scores.tsv  +  {pgs_id}_chr{N}_metadata.json",
        "tool": "compute_prs_spark_gpu.py",
        "slurm_script": "run_prs_spark_gpu.sbatch",
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
            "sbatch --export=PGS_ID={pgs_id},CHROM=1 run_prs_spark_gpu.sbatch\n"
            "# or array job for all chromosomes:\n"
            "sbatch --array=1-22 --export=PGS_ID={pgs_id} prs_gpu_compute.sbatch"
        ),
        "notes": (
            "gpu_batch_matches=4000 avoids OOM on RTX 3060 (12 GB VRAM). "
            "CPU fallback: compute_prs.py with DuckDB index."
        ),
    },
    {
        "step": 3,
        "name": "Aggregate total PRS",
        "description": (
            "Sum all per-chromosome PRS values per sample. "
            "Output has 24 columns: sample_id, PRS_total, PRS_chr1 … PRS_chr22."
        ),
        "output_pattern": "{pgs_id}_PRS_total.tsv  +  {pgs_id}_PRS_total_metadata.json",
        "tool": "aggregation step (post-chromosome scoring)",
        "slurm_script": None,
        "notes": "PRS_total[s] = Σ_chr PRS_chr[s] across all 22 autosomes.",
    },
    {
        "step": 4,
        "name": "Export to Parquet (recommended — enables web dashboard)",
        "description": (
            "Convert {pgs_id}_PRS_total.tsv to {pgs_id}_PRS_total.parquet "
            "using DuckDB. This step is NOT yet in the pipeline scripts and "
            "must be added after step 3."
        ),
        "output_pattern": "{pgs_id}_PRS_total.parquet",
        "tool": "DuckDB COPY TO",
        "slurm_script": None,
        "command_example": (
            "duckdb -c \""
            "COPY (SELECT * FROM read_csv_auto('{scores_dir}/{pgs_id}/{pgs_id}_PRS_total.tsv', "
            "delim=chr(9), header=true)) "
            "TO '{scores_dir}/{pgs_id}/{pgs_id}_PRS_total.parquet' "
            "(FORMAT PARQUET, COMPRESSION ZSTD)\""
        ),
        "notes": (
            "ZSTD compression gives ~3× size reduction over uncompressed TSV. "
            "DuckDB reads both TSV and parquet; parquet preferred for speed."
        ),
    },
]


class PipelineInspector:
    def __init__(self, cfg: Config) -> None:
        self.cfg = cfg

    # ── Public API ───────────────────────────────────────────────────────────

    def get_pipeline_plan(self, pgs_id: str) -> dict:
        """
        Return a structured processing plan for pgs_id.
        No jobs are submitted or files modified.
        """
        if not validate_pgs_id(pgs_id):
            return {"error": "Invalid PGS ID format", "pgs_id": pgs_id}

        scores_dir = str(self.cfg.scores_dir)
        pgs_scores_dir = str(self.cfg.scores_dir / pgs_id)

        def sub(s: str) -> str:
            return (
                s.replace("{pgs_id}", pgs_id)
                 .replace("{N}", "<chrom>")
                 .replace("{scores_dir}", scores_dir)
            )

        steps = []
        for tpl in _PLAN_TEMPLATE:
            step = json.loads(json.dumps(tpl))  # deep copy
            for key in ("description", "output_pattern", "submit_example", "command_example"):
                if key in step:
                    step[key] = sub(step[key])
            if step.get("submit_example"):
                step["submit_example"] = sub(step["submit_example"])
            steps.append(step)

        return {
            "pgs_id": pgs_id,
            "pgs_scores_dir": pgs_scores_dir,
            "zarr_base": "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files/chr{N}.zarr",
            "duckdb_index": "/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zarr_index.duckdb",
            "steps": steps,
            "scripts_dir": str(self.cfg.SCRIPTS_DIR),
            "scripts_available": self.list_scripts(),
            "note": (
                "This plan is informational only. "
                "No jobs have been submitted. "
                "Verify all paths and resource requirements before running."
            ),
        }

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
        Returns None if the script is not in the allowlist or doesn't exist.
        Security: allowlist prevents any path traversal.
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
