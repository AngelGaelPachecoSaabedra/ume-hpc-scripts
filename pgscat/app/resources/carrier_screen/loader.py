"""
carrier_screen/loader.py
========================
Loader for the Infinium GDA Carrier Screen gene list.

BACKGROUND
----------
The Illumina Infinium Global Diversity Array (GDA) with Carrier Screening content
covers ~602 genes associated with recessive conditions, as specified by ACMG
carrier screening recommendations (Carrier Screening in the Age of Genomic Medicine,
ACMG 2021 and carrier screening guidelines).

There are two authoritative routes to obtain the gene list:

  Route A — Illumina Manifest CSV
    The official manifest for the GDA array (Infinium_GlobalDiversity_Array-8
    v1.0 Kit) is available from Illumina. The manifest file contains:
      - Variant identifiers (SNP, indel)
      - Gene symbols
      - Chromosomal coordinates
    The gene list can be extracted by selecting unique ILMN_Strand-bearing
    rows and collecting the unique values in the 'Name' / gene-symbol column.

    File expected at:
      resources/carrier_screen/gda_carrier_manifest.csv
    or as specified in MANIFEST_CSV_PATH below.

  Route B — ACMG Carrier Screening Recommendations
    The ACMG 2023 secondary findings list and the 2021 Carrier Screening
    guidelines cover an overlapping but not identical set of genes.
    These can be used to validate or supplement the manifest-derived list.
    See: https://www.acmg.net/ACMG/Medical-Genetics-Practice-Resources/Practice-Guidelines.aspx

IMPORTANT: The gene list is NOT fabricated here. This loader only reads from
authoritative files. If no file is available, it returns an empty structure
with documentation about what is needed.

USAGE
-----
    from resources.carrier_screen.loader import CarrierScreenLoader
    loader = CarrierScreenLoader()
    genes  = loader.get_genes()   # list of {"gene_symbol": str, "n_variants": int, ...}
    loader.is_available()         # True if a gene list file has been loaded
"""
from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

# ── Paths ──────────────────────────────────────────────────────────────────────
_HERE = Path(__file__).parent

# Primary: Illumina manifest CSV (not included — must be obtained from Illumina)
MANIFEST_CSV_PATH = _HERE / "gda_carrier_manifest.csv"

# Secondary: curated gene list JSON (can be built from manifest or ACMG guidelines)
GENE_LIST_JSON_PATH = _HERE / "carrier_genes.json"

# Schema reference
SCHEMA_PATH = _HERE / "schema.json"


class CarrierScreenLoader:
    """
    Loads and exposes the carrier screen gene list.

    Loading priority:
      1. carrier_genes.json  (prebuilt from manifest or curated manually)
      2. gda_carrier_manifest.csv  (raw Illumina manifest — parsed on-the-fly)
      3. Empty list with documentation if neither is available
    """

    def __init__(self) -> None:
        self._genes: Optional[List[Dict[str, Any]]] = None
        self._source: Optional[str] = None
        self._load()

    def _load(self) -> None:
        """Load gene list from the highest-priority available source."""
        # Priority 1: prebuilt JSON
        if GENE_LIST_JSON_PATH.exists():
            try:
                data = json.loads(GENE_LIST_JSON_PATH.read_text(encoding="utf-8"))
                if isinstance(data, dict):
                    self._genes = data.get("genes", [])
                elif isinstance(data, list):
                    self._genes = data
                self._source = "carrier_genes.json"
                logger.info(
                    "Carrier screen: loaded %d genes from %s",
                    len(self._genes), GENE_LIST_JSON_PATH
                )
                return
            except Exception as exc:
                logger.warning("Cannot load carrier_genes.json: %s", exc)

        # Priority 2: raw Illumina manifest CSV
        if MANIFEST_CSV_PATH.exists():
            try:
                self._genes = self._parse_manifest(MANIFEST_CSV_PATH)
                self._source = "gda_carrier_manifest.csv"
                logger.info(
                    "Carrier screen: parsed %d unique genes from manifest %s",
                    len(self._genes), MANIFEST_CSV_PATH
                )
                return
            except Exception as exc:
                logger.warning("Cannot parse GDA manifest CSV: %s", exc)

        # Priority 3: not available
        logger.info(
            "Carrier screen gene list not available. "
            "See resources/carrier_screen/README.md for instructions."
        )
        self._genes = []
        self._source = None

    @staticmethod
    def _parse_manifest(csv_path: Path) -> List[Dict[str, Any]]:
        """
        Parse the Illumina GDA manifest CSV and extract unique gene symbols.

        The Illumina manifest CSV has a [Assay] section and a [Controls] section.
        Data rows start after the header row that contains 'IlmnID' or 'Name'.

        Expected columns (may vary by manifest version):
          IlmnID, Name, IlmnStrand, SNP, AddressA_ID, AlleleA_ProbeSeq,
          AddressB_ID, AlleleB_ProbeSeq, GenomeBuild, Chr, MapInfo,
          Ploidy, Species, Source, SourceVersion, SourceStrand, SourceSeq,
          TopGenomicSeq, BeadSetID, Exp_Clusters, Intensity_Only,
          RefStrand, ILMN_Gene, ...
        """
        gene_counts: Dict[str, int] = {}
        in_data = False
        col_gene: Optional[int] = None
        col_chr: Optional[int] = None

        with open(csv_path, encoding="utf-8", errors="replace") as fh:
            reader = csv.reader(fh)
            for row in reader:
                if not row:
                    continue

                # Detect header row: contains "IlmnID" or "Name" and "ILMN_Gene" / "Chr"
                if not in_data:
                    row_lower = [c.strip().lower() for c in row]
                    if "ilmnid" in row_lower or "name" in row_lower:
                        in_data = True
                        # Find gene column
                        for candidate in ("ilmn_gene", "gene_name", "gene", "symbol"):
                            if candidate in row_lower:
                                col_gene = row_lower.index(candidate)
                                break
                        # Find chr column
                        for candidate in ("chr", "chromosome", "chrom"):
                            if candidate in row_lower:
                                col_chr = row_lower.index(candidate)
                                break
                        # If we still haven't found a gene column, look for anything with 'gene'
                        if col_gene is None:
                            for i, c in enumerate(row_lower):
                                if "gene" in c:
                                    col_gene = i
                                    break
                    continue

                # Skip section markers
                if row[0].startswith("["):
                    in_data = False
                    continue

                # Extract gene
                if col_gene is not None and len(row) > col_gene:
                    gene = row[col_gene].strip()
                    if gene and gene not in ("", ".", "NA", "nan", "-"):
                        gene_counts[gene.upper()] = gene_counts.get(gene.upper(), 0) + 1

        return [
            {"gene_symbol": g, "n_variants_in_manifest": n}
            for g, n in sorted(gene_counts.items())
        ]

    # ── Public API ──────────────────────────────────────────────────────────────

    def is_available(self) -> bool:
        """Return True if a gene list has been successfully loaded."""
        return bool(self._genes)

    def get_genes(self) -> List[Dict[str, Any]]:
        """
        Return the list of carrier screen genes.

        Each entry is a dict with at least:
          {"gene_symbol": str}
        and optionally:
          {"n_variants_in_manifest": int, "conditions": [...], "acmg_tier": str}
        """
        return list(self._genes or [])

    def get_gene_symbols(self) -> List[str]:
        """Return just the gene symbols as a sorted list."""
        return sorted(g["gene_symbol"] for g in self.get_genes())

    def is_carrier_gene(self, gene_symbol: str) -> bool:
        """Return True if gene_symbol is in the carrier screen panel."""
        return gene_symbol.upper() in {g["gene_symbol"].upper() for g in self.get_genes()}

    @property
    def source(self) -> Optional[str]:
        """Name of the source file that was loaded."""
        return self._source

    @property
    def n_genes(self) -> int:
        """Number of unique genes in the carrier screen panel."""
        return len(self._genes or [])

    def status(self) -> Dict[str, Any]:
        """Return a status dict for API/UI display."""
        return {
            "available":  self.is_available(),
            "n_genes":    self.n_genes,
            "source":     self._source,
            "expected_sources": [
                str(GENE_LIST_JSON_PATH),
                str(MANIFEST_CSV_PATH),
            ],
            "instructions": (
                "To populate the carrier screen gene list, either:\n"
                "  1. Place the Illumina GDA manifest CSV at:\n"
                f"     {MANIFEST_CSV_PATH}\n"
                "     (obtain from Illumina support or product page)\n"
                "  2. Or create carrier_genes.json from the manifest or ACMG guidelines.\n"
                "     See resources/carrier_screen/README.md for the expected JSON schema."
            ),
        }

    def save_genes_json(self, genes: List[Dict[str, Any]], source_note: str = "") -> None:
        """
        Persist a list of gene dicts to carrier_genes.json.
        Called after parsing the manifest to cache the result.

        Parameters
        ----------
        genes : list of dict
            Each dict must have at least {"gene_symbol": str}
        source_note : str
            Free-text note about how the list was generated.
        """
        payload = {
            "source_note": source_note,
            "n_genes": len(genes),
            "genes": genes,
        }
        GENE_LIST_JSON_PATH.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False),
            encoding="utf-8"
        )
        logger.info("Saved %d carrier screen genes to %s", len(genes), GENE_LIST_JSON_PATH)
        self._genes = genes
        self._source = "carrier_genes.json"


# ── Module-level singleton ─────────────────────────────────────────────────────
# Instantiated lazily on first import to avoid blocking startup.
_loader_instance: Optional[CarrierScreenLoader] = None


def get_loader() -> CarrierScreenLoader:
    """Return the module-level singleton loader."""
    global _loader_instance
    if _loader_instance is None:
        _loader_instance = CarrierScreenLoader()
    return _loader_instance
