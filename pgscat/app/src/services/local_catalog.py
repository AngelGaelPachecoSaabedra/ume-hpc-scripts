"""
LocalCatalog service
────────────────────
Scans the scores directory (read-only mount) and exposes structured
information about locally available PGS datasets.

Directory layout assumed:
  $APP_DATA_DIR/
    PGS000004/
      PGS000004.json                          ← PGS Catalog metadata
      PGS000004_hmPOS_GRCh38.betamap.tsv.gz  ← input weights
      PGS000004_chr1_scores.tsv              ← per-chrom scores
      PGS000004_chr1_metadata.json
      ...
      PGS000004_PRS_total.tsv                ← aggregated scores (24 cols)
      PGS000004_PRS_total_metadata.json      ← aggregation metadata
      PGS000004_PRS_total.parquet            ← preferred for web queries
"""
import json
import re
import logging
from pathlib import Path
from typing import Optional

from config import Config

logger = logging.getLogger(__name__)

# PGS IDs: PGS followed by exactly 6 digits (e.g. PGS000001 ... PGS999999)
_PGS_ID_RE = re.compile(r"^PGS\d{6}$")

CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X"]


def validate_pgs_id(pgs_id: str) -> bool:
    """Return True if pgs_id has the canonical PGS ID format."""
    return bool(_PGS_ID_RE.match(pgs_id))


class LocalCatalog:
    def __init__(self, cfg: Config) -> None:
        self.cfg = cfg

    # ── Public API ───────────────────────────────────────────────────────────

    def list_pgs(self) -> list[dict]:
        """
        Scan scores_dir and return a list of quick-info dicts for every
        valid PGS subdirectory found.  Does not raise on missing/empty dir.
        """
        scores_dir = self.cfg.scores_dir
        if not scores_dir.exists():
            logger.warning("scores_dir does not exist: %s", scores_dir)
            return []

        results = []
        for entry in sorted(scores_dir.iterdir()):
            if entry.is_dir() and validate_pgs_id(entry.name):
                try:
                    results.append(self._quick_info(entry))
                except Exception as exc:
                    logger.error("Error reading %s: %s", entry.name, exc)
        return results

    def get_pgs_info(self, pgs_id: str) -> Optional[dict]:
        """
        Return full info dict for a single PGS ID, or None if not found.
        Includes metadata, file listing, total_metadata, and chromosome list.
        """
        if not validate_pgs_id(pgs_id):
            return None
        pgs_dir = self._pgs_dir(pgs_id)
        if not pgs_dir.is_dir():
            return None

        info = self._quick_info(pgs_dir)
        info["total_metadata"] = self._load_json(pgs_dir / f"{pgs_id}_PRS_total_metadata.json")
        info["pgs_dir"] = str(pgs_dir)
        info["files"] = self._list_files(pgs_dir)
        return info

    def exists_locally(self, pgs_id: str) -> bool:
        """True if the PGS directory exists."""
        return validate_pgs_id(pgs_id) and self._pgs_dir(pgs_id).is_dir()

    def has_results(self, pgs_id: str) -> bool:
        """True if at least one results file (parquet or total TSV) exists."""
        if not validate_pgs_id(pgs_id):
            return False
        pgs_dir = self._pgs_dir(pgs_id)
        return (
            (pgs_dir / f"{pgs_id}_PRS_total.parquet").exists()
            or (pgs_dir / f"{pgs_id}_PRS_total.tsv").exists()
        )

    # ── Internal helpers ─────────────────────────────────────────────────────

    def _pgs_dir(self, pgs_id: str) -> Path:
        return self.cfg.scores_dir / pgs_id

    def _quick_info(self, pgs_dir: Path) -> dict:
        pgs_id = pgs_dir.name
        meta = self._load_json(pgs_dir / f"{pgs_id}.json")
        chromosomes = self._detect_chromosomes(pgs_dir, pgs_id)
        has_parquet = (pgs_dir / f"{pgs_id}_PRS_total.parquet").exists()
        has_tsv = (pgs_dir / f"{pgs_id}_PRS_total.tsv").exists()
        return {
            "pgs_id": pgs_id,
            # trait_reported is canonical; fall back to name, then Unknown
            "trait_name": meta.get("trait_reported") or meta.get("name") or "Unknown",
            "trait_efo": meta.get("trait_efo", []),
            "n_variants": meta.get("variants_number"),
            "has_parquet": has_parquet,
            "has_tsv": has_tsv,
            "has_results": has_parquet or has_tsv,
            "chromosomes": chromosomes,
            "n_chromosomes": len(chromosomes),
            "meta_available": bool(meta),
        }

    def _detect_chromosomes(self, pgs_dir: Path, pgs_id: str) -> list[str]:
        found = []
        for c in CHROMOSOMES:
            if (pgs_dir / f"{pgs_id}_chr{c}_scores.tsv").exists():
                found.append(c)
        return found

    def _list_files(self, pgs_dir: Path) -> list[dict]:
        files = []
        for f in sorted(pgs_dir.iterdir()):
            if f.is_file():
                try:
                    size_mb = round(f.stat().st_size / 1_048_576, 2)
                except OSError:
                    size_mb = None
                files.append({
                    "name": f.name,
                    "size_mb": size_mb,
                    "suffix": f.suffix,
                })
        return files

    def _load_json(self, path: Path) -> dict:
        try:
            if path.exists():
                with path.open() as fh:
                    return json.load(fh)
        except Exception as exc:
            logger.debug("Could not load %s: %s", path, exc)
        return {}
