# Carrier Screen Gene List — Resources

## Overview

This directory holds the gene list for the **Infinium GDA Carrier Screen**
(Illumina Infinium Global Diversity Array with Carrier Screening content, ~602 genes).

The gene list is used by the PGS/PRS Platform to:
- Tag variants in carrier screen genes in the gene browser
- Filter/highlight carrier screen genes in the variant annotation view
- Support a future dedicated carrier screen view at `/carrier-screen/<gene>`

---

## Status

> **Gene list NOT yet populated.**
>
> The list of ~602 genes must be obtained from an authoritative source (see below).
> This directory contains only the loading infrastructure, schema, and documentation.
> No gene list is fabricated or invented.

---

## How to Populate the Gene List

### Option A — Illumina GDA Manifest CSV (Recommended)

1. Obtain the manifest file from Illumina:
   - Product: **Infinium Global Diversity Array-8 Kit** (or equivalent)
   - File: `Infinium_GlobalDiversity_Array-8_v1.0_<build>.csv`
   - Source: Illumina support portal or product download page
   - The manifest contains all probe positions, gene symbols, SNP IDs, etc.

2. Place the file at:
   ```
   resources/carrier_screen/gda_carrier_manifest.csv
   ```

3. Run the extraction script to build `carrier_genes.json`:
   ```python
   from resources.carrier_screen.loader import CarrierScreenLoader
   loader = CarrierScreenLoader()
   # If manifest is present, loader will auto-parse it.
   # Then save as prebuilt JSON for faster future loads:
   loader.save_genes_json(loader.get_genes(), source_note="GDA manifest v1.0, 2026-04-13")
   ```
   This produces `resources/carrier_screen/carrier_genes.json`.

4. Verify:
   ```python
   from resources.carrier_screen.loader import get_loader
   loader = get_loader()
   print(loader.n_genes, "genes loaded from", loader.source)
   ```

### Option B — ACMG Carrier Screening Recommendations

The ACMG publishes carrier screening guidelines that overlap significantly with
the GDA panel. To use these:

1. Reference: *ACMG SF v3.2 (2023)* and *"Carrier Screening in the Age of Genomic Medicine" (ACMG, 2021)*
2. Manually curate or download the ACMG carrier screening gene list
3. Create `carrier_genes.json` following the schema in `schema.json`
4. Include `acmg_tier` and `inheritance` fields for each gene

---

## File Layout

```
resources/carrier_screen/
  README.md                 ← This file
  loader.py                 ← Python loader (auto-parses manifest OR JSON)
  schema.json               ← JSON schema for carrier_genes.json
  carrier_genes.json        ← Prebuilt gene list (NOT included; generate from manifest)
  gda_carrier_manifest.csv  ← Raw Illumina manifest (NOT included; obtain from Illumina)
```

---

## carrier_genes.json Format

See `schema.json` for the full schema. Minimum valid format:

```json
{
  "source_note": "Extracted from Illumina GDA manifest v1.0 on YYYY-MM-DD",
  "n_genes": 602,
  "genes": [
    {
      "gene_symbol": "CFTR",
      "n_variants_in_manifest": 42,
      "conditions": ["Cystic fibrosis"],
      "inheritance": "AR",
      "acmg_tier": "Tier 1",
      "chromosome": "7"
    }
  ]
}
```

---

## Integration with the Platform

Once `carrier_genes.json` is populated:

- **Gene browser** (`/gene/<symbol>`): carrier screen genes will display a badge
- **API endpoint** (`/api/carrier-screen/genes`): returns the full gene list as JSON
- **Variant annotation table**: variants in carrier screen genes can be filtered
- **Future view** (`/carrier-screen`): dedicated carrier screen interpretation page

### Planned routes (future implementation)

```
GET  /carrier-screen           — carrier screen overview
GET  /carrier-screen/<gene>    — gene detail (redirects to /gene/<gene> with CS badge)
GET  /api/carrier-screen/genes — list of all carrier screen genes
GET  /api/carrier-screen/status — loader status JSON
```

---

## Gene Count Note

The ~602 gene count comes from the Infinium GDA array content specification.
The exact number depends on the manifest version and may vary slightly between
array versions. Do not use a hardcoded 602 — use the count from the manifest.

---

## References

- Illumina GDA product page: https://www.illumina.com/products/by-type/microarray-kits/infinium-global-diversity-array.html
- ACMG Carrier Screening Guidelines: https://www.acmg.net
- ACMG SF v3.2: https://www.acmg.net/ACMG/Medical-Genetics-Practice-Resources/Practice-Guidelines.aspx
