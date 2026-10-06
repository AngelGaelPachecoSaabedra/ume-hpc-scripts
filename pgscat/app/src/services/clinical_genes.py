"""
clinical_genes.py
=================
Runtime service for clinical genes lookup.

Loads resources/clinical_genes/clinical_genes_v1.json (built by clinical_genes_builder.py)
and exposes fast in-memory lookup methods used by GeneBrowserService and
VariantAnnotationService.

Usage:
    from services.clinical_genes import get_clinical_genes_service

    svc = get_clinical_genes_service()
    svc.is_clinical_gene("CFTR")          # True
    svc.get_gene_info("CFTR")             # {"gene_symbol": "CFTR", "sources": [...], ...}
    svc.get_clinical_sources("CFTR")      # ["LocalExcel_GenCC_AR", "GenCC"]
    svc.get_confidence("CFTR")            # "high"
    svc.status()                          # {"available": True, "n_genes": 612, ...}
"""
from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

logger = logging.getLogger(__name__)

# ── Path resolution ────────────────────────────────────────────────────────────
# Priority:
#   1. APP_CLINICAL_GENES_JSON env var (explicit override)
#   2. Default path relative to this file   (/app/resources/…)
#   3. Inferred from APP_ANNOTATIONS_DIR  (container w/ CephFS mount but no COPY)
_HERE      = Path(__file__).parent          # src/services/
_APP_ROOT  = _HERE.parent.parent            # app/
_DEFAULT   = _APP_ROOT / "resources" / "clinical_genes" / "clinical_genes_v1.json"

_env_path = os.environ.get("APP_CLINICAL_GENES_JSON", "").strip()
if _env_path:
    _JSON_PATH = Path(_env_path)
elif _DEFAULT.exists():
    _JSON_PATH = _DEFAULT
else:
    # Infer resources path relative to APP_ANNOTATIONS_DIR
    _ann = os.environ.get("APP_ANNOTATIONS_DIR", "").strip()
    _candidate = (
        Path(_ann).parent / "app" / "resources" / "clinical_genes" / "clinical_genes_v1.json"
        if _ann else _DEFAULT
    )
    _JSON_PATH = _candidate if _candidate.exists() else _DEFAULT


class ClinicalGenesService:
    """
    In-memory clinical genes lookup service.

    Loads from clinical_genes_v1.json on first instantiation.
    Returns empty/safe values if the JSON is not yet built.
    """

    def __init__(self, json_path: Path = _JSON_PATH) -> None:
        self._json_path = json_path
        self._genes: Dict[str, Dict[str, Any]] = {}   # symbol → entry
        self._available = False
        self._n_genes = 0
        self._load()

    def _load(self) -> None:
        if not self._json_path.exists():
            logger.info(
                "Clinical genes JSON not found at %s. "
                "Run scripts/download_clinical_genes.sh then src/services/clinical_genes_builder.py",
                self._json_path,
            )
            return

        try:
            data = json.loads(self._json_path.read_text(encoding="utf-8"))
            genes_list: List[Dict[str, Any]] = []

            if isinstance(data, dict):
                genes_list = data.get("genes", [])
            elif isinstance(data, list):
                genes_list = data

            for entry in genes_list:
                sym = entry.get("gene_symbol", "").strip().upper()
                if sym:
                    self._genes[sym] = entry

            self._n_genes   = len(self._genes)
            self._available = self._n_genes > 0
            logger.info("ClinicalGenesService: loaded %d genes from %s", self._n_genes, self._json_path)

        except Exception as exc:
            logger.warning("ClinicalGenesService: cannot load %s: %s", self._json_path, exc)

    # ── Public API ──────────────────────────────────────────────────────────────

    def is_available(self) -> bool:
        """True if the clinical genes dataset has been loaded."""
        return self._available

    @property
    def n_genes(self) -> int:
        return self._n_genes

    def is_clinical_gene(self, gene_symbol: str) -> bool:
        """Return True if gene_symbol is in the clinical genes dataset."""
        if not gene_symbol:
            return False
        return gene_symbol.strip().upper() in self._genes

    def get_gene_info(self, gene_symbol: str) -> Optional[Dict[str, Any]]:
        """
        Return the full clinical gene entry for gene_symbol, or None.

        Entry format:
          {
            "gene_symbol": "CFTR",
            "sources": ["LocalExcel_GenCC_AR", "GenCC"],
            "confidence": "high",
            "moi": "Autosomal recessive",
            "evidence": "definitive",
            "disease": "cystic fibrosis",
            "hgnc_id": "HGNC:1884",
            "in_clinvar": true
          }
        """
        if not gene_symbol:
            return None
        return self._genes.get(gene_symbol.strip().upper())

    def get_clinical_sources(self, gene_symbol: str) -> List[str]:
        """Return list of source labels for gene_symbol, or []."""
        entry = self.get_gene_info(gene_symbol)
        if not entry:
            return []
        sources = entry.get("sources", [])
        return list(sources) if isinstance(sources, list) else [str(sources)]

    def get_confidence(self, gene_symbol: str) -> Optional[str]:
        """Return confidence tier ('high', 'medium', 'standard') or None."""
        entry = self.get_gene_info(gene_symbol)
        return entry.get("confidence") if entry else None

    def get_all_symbols(self) -> Set[str]:
        """Return the full set of clinical gene symbols."""
        return set(self._genes.keys())

    def status(self) -> Dict[str, Any]:
        """Return a status dict suitable for API/UI display."""
        return {
            "available":   self._available,
            "n_genes":     self._n_genes,
            "json_path":   str(self._json_path),
            "json_exists": self._json_path.exists(),
            "instructions": (
                "To build the clinical genes dataset:\n"
                "  1. bash scripts/download_clinical_genes.sh\n"
                "  2. python3 src/services/clinical_genes_builder.py"
            ) if not self._available else None,
        }


# ── Module-level singleton ─────────────────────────────────────────────────────
_service_instance: Optional[ClinicalGenesService] = None


def get_clinical_genes_service() -> ClinicalGenesService:
    """Return the module-level singleton ClinicalGenesService."""
    global _service_instance
    if _service_instance is None:
        _service_instance = ClinicalGenesService()
    return _service_instance
