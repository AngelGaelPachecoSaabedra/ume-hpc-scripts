# Variant Annotator — Module Documentation (v1.3)

## Overview

The Variant Annotator is an independent HPC pipeline module integrated into the PGS/PRS Platform.
It annotates each variant in a betamap file using:

- **GENCODE GFF3** (required) — region classification, splice-site distance, all-gene overlaps
- **Reference FASTA** (optional, `--fasta`) — FASTA-backed coding consequences (missense,
  synonymous, stop_gained, stop_lost, start_lost, frameshift, inframe indels) and
  splice donor/acceptor classification via GT/AG inspection
- **dbNSFP5** (optional, `--dbnsfp`) — CADD_phred, REVEL_score, SIFT_pred,
  Polyphen2_HDIV_pred, clinvar_clnsig via tabix queries
- **dbSNP population frequency VCF** (optional, `--dbsnp`) — rsid, allele frequencies,
  per-population AFs, rarity classification (v1.3 new)
- **Regulatory BED files** (optional, `--regulatory-bed`) — ENCODE cCREs, Ensembl Regulatory,
  promoters, enhancers, open chromatin

**Key design principle:** the annotation pipeline runs _outside_ the web container, on HPC nodes
via Apptainer. The Bottle web app is read-only with respect to annotation.

---

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    PGS/PRS Platform (Bottle)                │
│                                                             │
│  /variants/<pgs_id>              HTML view (variants.tpl)   │
│  /api/variants/<pgs_id>          paginated variant JSON     │
│  /api/variants/<pgs_id>/summary  annotation summary JSON    │
│  /api/variants/<pgs_id>/ideogram chromosome-level data      │
│                                                             │
│  VariantAnnotationService  ←── reads filesystem only        │
└───────────────────────┬─────────────────────────────────────┘
                        │ reads
                        ▼
          /annotations/{PGS_ID}/
              {PGS_ID}_variants_annotated.tsv.gz
              {PGS_ID}_variants_annotated.parquet
              {PGS_ID}_annotation_summary.json
              annotation.log
                        ▲ writes
                        │
┌───────────────────────┴─────────────────────────────────────┐
│                 HPC Node (Apptainer)                        │
│                                                             │
│  run_annotation.sh PGS000001                                │
│    └─ apptainer exec variant_annotator.sif                  │
│         └─ annotate_variants.py                             │
│              ├─ gff3_parser.py   (GFF3 → IntervalTree,      │
│              │                    transcript_cds_map)        │
│              ├─ classify.py      (region_class, consequence, │
│              │                    splice distance/type,      │
│              │                    codon/AA/dbNSFP fields)    │
│              ├─ fasta_engine.py  (codon extraction, splice   │
│              │                    GT/AG, REF/ALT FASTA)      │
│              ├─ dbnsfp.py        (CADD/REVEL/SIFT tabix)     │
│              ├─ sequence_utils.py (codon table, translate)   │
│              └─ regulatory.py    (BED → IntervalTree)       │
└─────────────────────────────────────────────────────────────┘

External data (bind-mounted at runtime, not in container):
  /path/to/hg38.fa            ← FASTA
  /path/to/dbNSFP5.0a_grch38.gz ← dbNSFP5
  /path/to/gencode.annotation.gff3
```

---

## File Layout

```
app/
├── annotator/
│   ├── __init__.py
│   ├── gff3_parser.py          GFF3 stream parser + interval tree index
│   │                           Stores begin/end in GFF3Feature;
│   │                           transcript_cds_map for codon lookups;
│   │                           nearest_exon_boundary_distance(),
│   │                           get_cds_offset().
│   ├── classify.py             Variant classification logic
│   │                           Region/consequence; splice distance + type;
│   │                           codon/AA consequence (v1.2);
│   │                           dbNSFP5 score fields; is_lof/missense/synonymous.
│   ├── sequence_utils.py       Codon table, translation, reverse complement,
│   │                           snp_codon_consequence, indel_consequence,
│   │                           is_lof_consequence.
│   ├── fasta_engine.py         FASTA-backed consequence engine.
│   │                           get_codon_consequence(), get_splice_type(),
│   │                           multi-exon codon extraction.
│   ├── dbnsfp.py               dbNSFP5 tabix fetcher (pysam).
│   │                           CADD, REVEL, SIFT, PolyPhen2, ClinVar.
│   ├── dbsnp_freq.py           dbSNP population frequency VCF fetcher (v1.3).
│   │                           rsid, af_global, af_max_population,
│   │                           af_population_summary, rarity_class.
│   ├── regulatory.py           Regulatory BED interval index.
│   ├── annotate_variants.py    Main pipeline entry point (CLI v1.3)
│   └── run_annotation.sh       Apptainer wrapper script (v1.3)
├── apptainer/
│   └── variant_annotator.def   Apptainer definition file (v1.3.0)
├── docs/
│   ├── VARIANT_ANNOTATION.md   This document
│   └── ARCHITECTURE.md         Platform architecture
├── resources/
│   └── carrier_screen/         Carrier screen gene list infrastructure
│       ├── README.md            Instructions for populating the gene list
│       ├── loader.py            Python loader (manifest CSV or carrier_genes.json)
│       └── schema.json          JSON schema for carrier_genes.json
└── src/
    ├── services/
    │   ├── variant_annotation.py   Web-layer service (read-only, v1.3 dbSNP cols)
    │   └── gene_browser.py         Gene-centric browser service (Mode E)
    └── views/
        ├── variants.tpl            HTML template for annotation view
        ├── gene.tpl                Gene browser (Plotly stacked tracks, no MTR)
        └── genes.tpl               Gene list / search UI
```

---

## Building the Container

```bash
apptainer build \
    /path/to/variant_annotator.sif \
    apptainer/variant_annotator.def
```

**Container contents (v1.2.0):**
- Python 3.12
- pandas 2.2, numpy 1.26, pyarrow 15 (conda-forge)
- intervaltree 3.1.0
- pyfaidx 0.8.1.2 — random-access FASTA reader
- pysam 0.22.1 — tabix interface for dbNSFP5
- polars 0.20.18, tqdm 4.66

Source code is bind-mounted at `/app/annotator` — no container rebuild needed for code changes.

---

## Running Annotation

### Full v1.3 (FASTA + dbNSFP5 + dbSNP — recommended)

```bash
DATA_DIR=/path/to/data \
ANNOTATIONS_DIR=/annotations \
GFF3_PATH=/path/to/gencode.annotation.gff3 \
SIF_PATH=/path/to/variant_annotator.sif \
FASTA_PATH=/path/to/hg38.fa \
DBNSFP_PATH=/path/to/dbNSFP5.0a_grch38.gz \
DBSNP_FREQ_PATH=/path/to/dbsnp_freq.vcf.gz \
/path/to/annotator/run_annotation.sh PGS000001
```

### Minimal (GFF3 only)

```bash
DATA_DIR=/path/to/data \
ANNOTATIONS_DIR=/annotations \
GFF3_PATH=/path/to/gencode.annotation.gff3 \
SIF_PATH=/path/to/variant_annotator.sif \
/path/to/annotator/run_annotation.sh PGS000001
```

Without `--fasta`, CDS variants are reported as `coding_sequence_variant`. Without
`--dbnsfp`, functional score columns are empty.

### With regulatory BED files

```bash
REGULATORY_BED="/path/to/encode_cCRE.bed.gz \
               /path/to/ensembl_reg.bed.gz" \
... \
/path/to/annotator/run_annotation.sh PGS000001
```

### Via Slurm

```bash
sbatch \
    --partition=highmem \
    --mem=64G \
    --cpus-per-task=4 \
    --job-name=annotate_PGS000001 \
    --output=/annotations/PGS000001/slurm-%j.log \
    /path/to/annotator/run_annotation.sh PGS000001
```

---

## Inputs

### Betamap file

| Column        | Type    | Description                              |
|---------------|---------|------------------------------------------|
| PRS_ID        | string  | PGS identifier (e.g. PGS000001)          |
| CHROM         | string  | Chromosome (bare: 1..22, X, Y, MT)       |
| POS           | int64   | Position, 1-based GRCh38                 |
| ID            | string  | rsID or variant identifier               |
| EFFECT_ALLELE | string  | Effect allele                            |
| OTHER_ALLELE  | string  | Reference/other allele                   |
| BETA          | float64 | Effect size                              |
| IS_FLIP       | int8    | 1 = dosage flipped; 0 = standard         |

### GFF3 reference (required)

- **File:** `/path/to/gencode.annotation.gff3`
- **Release:** GENCODE v49 basic
- **Assembly:** GRCh38

### FASTA reference (optional, `--fasta` / `FASTA_PATH`)

- **File:** `/path/to/hg38.fa`
- **Purpose:** Determines REF vs ALT by FASTA lookup; extracts reference codon; applies ALT substitution in coding-strand orientation; translates and classifies SNP and indel coding consequences; classifies splice donor/acceptor by reading GT/AG dinucleotides.
- **FASTA chromosome convention:** hg38.fa uses `chr` prefix (chr1…chrX, chrM). Bare chromosome names in betamap (1, X, MT) are mapped transparently.

### dbNSFP5 (optional, `--dbnsfp` / `DBNSFP_PATH`)

- **File:** `/path/to/dbNSFP5.0a_grch38.gz` (tabix-indexed)
- **Purpose:** Per-SNP functional scores. Only SNPs are queried (indels are not in dbNSFP5).
- **Columns fetched:** `CADD_phred`, `REVEL_score`, `SIFT_pred`, `Polyphen2_HDIV_pred`, `clinvar_clnsig`

### Regulatory BED files (optional, `--regulatory-bed`)

- Zero or more BED files (plain or `.gz`)
- Standard 3-column BED minimum; 4th column used as `element_type`
- Typical sources: ENCODE cCREs, Ensembl Regulatory Build, tissue-specific promoters

---

## Outputs

All outputs are written to `{ANNOTATIONS_DIR}/{PGS_ID}/`:

### `{PGS_ID}_variants_annotated.tsv.gz`

Tab-separated, gzip-compressed. All betamap columns plus annotation columns.

#### v1.0 annotation columns

| Column                 | Type    | Description                                           |
|------------------------|---------|-------------------------------------------------------|
| gene_name              | string  | HGNC gene symbol of best-hit gene                     |
| gene_id                | string  | Ensembl gene ID without version (ENSG...)             |
| gene_type              | string  | protein_coding, lncRNA, pseudogene, etc.              |
| transcript_id          | string  | Ensembl transcript ID (ENST...) of best hit           |
| transcript_type        | string  | Transcript biotype                                    |
| feature_type           | string  | Raw GFF3 feature driving classification               |
| region_class           | string  | See Region Classes below                              |
| consequence            | string  | SO-inspired consequence term                          |
| is_coding              | bool    | True if CDS overlap                                   |
| is_regulatory          | bool    | True if regulatory BED overlap or TSS heuristic       |
| is_intergenic          | bool    | True if no gene overlap                               |
| n_overlapping_genes    | int     | Number of distinct genes overlapping position         |
| strand                 | string  | Strand of the annotating feature (+/-/.)              |
| distance_nearest_gene  | int/NA  | bp to nearest gene (only for intergenic)              |

#### v1.1 annotation columns

| Column                      | Type    | Description                                      |
|-----------------------------|---------|--------------------------------------------------|
| distance_to_splice_site     | int/NA  | bp to nearest exon boundary (all variants)       |
| consequence_priority        | int     | Numeric impact rank — lower = higher impact      |
| all_overlapping_genes       | string  | Comma-separated gene names of all overlapping genes |
| all_overlapping_transcripts | string  | Comma-separated transcript IDs (best per gene)   |
| all_region_classes          | string  | Comma-separated region_class per gene            |
| regulatory_source           | string  | "heuristic", "bed:\<source\>", or ""             |

#### v1.2 annotation columns (require `--fasta` / `--dbnsfp`)

| Column              | Type      | Description                                           |
|---------------------|-----------|-------------------------------------------------------|
| codon_ref           | string    | Reference codon in coding-strand orientation, e.g. "ATG" |
| codon_alt           | string    | Alternate codon, e.g. "TTG"                          |
| aa_ref              | string    | 1-letter reference amino acid                        |
| aa_alt              | string    | 1-letter alternate amino acid                        |
| aa_ref_3            | string    | 3-letter reference amino acid, e.g. "Met"            |
| aa_alt_3            | string    | 3-letter alternate amino acid, e.g. "Leu"            |
| splice_type         | string    | "donor" \| "acceptor" \| "" (FASTA GT/AG inspection) |
| cadd_phred          | float/NA  | CADD phred-scaled pathogenicity score                |
| revel_score         | float/NA  | REVEL ensemble missense score (0–1)                  |
| sift_pred           | string    | "D" (deleterious) \| "T" (tolerated) \| ""           |
| polyphen2_pred      | string    | "D" (damaging) \| "P" (possibly) \| "B" (benign) \| "" |
| clinvar_clnsig      | string    | ClinVar clinical significance (first non-missing)    |
| is_missense         | bool      | True if consequence == missense_variant              |
| is_synonymous       | bool      | True if consequence == synonymous_variant            |
| is_lof              | bool      | True for stop_gained/lost, start_lost, frameshift, splice_donor/acceptor |

#### v1.3 annotation columns (require `--dbsnp` / `DBSNP_FREQ_PATH`)

Complementary to dbNSFP5 — both can be active simultaneously.

| Column                   | Type      | Description                                          |
|--------------------------|-----------|------------------------------------------------------|
| rsid                     | string/NA | dbSNP rsID (e.g. "rs123456"), or None                |
| af_global                | float/NA  | Global allele frequency (best estimate from VCF)     |
| af_max_population        | float/NA  | Maximum AF across any reported population            |
| af_population_summary    | string/NA | Comma-separated "POP=AF" pairs (e.g. "AFR=0.12,EUR=0.05") |
| rarity_class             | string    | "common" (≥1%) \| "low_frequency" (0.1–1%) \| "rare" (0.01–0.1%) \| "ultra_rare" (<0.01%) \| "novel" (not in dbSNP) |

**Rarity thresholds:**
- `common`: AF ≥ 0.01
- `low_frequency`: 0.001 ≤ AF < 0.01
- `rare`: 0.0001 ≤ AF < 0.001
- `ultra_rare`: 0 < AF < 0.0001
- `novel`: not found in dbSNP population frequency VCF

**dbSNP source file:** `/path/to/dbsnp_freq.vcf.gz`
(tabix-indexed; index at `freq.vcf.gz.tbi`)

### `{PGS_ID}_variants_annotated.parquet`

Same schema, Snappy-compressed Parquet. Preferred for web queries.

### `{PGS_ID}_annotation_summary.json`

```json
{
  "pgs_id": "PGS000001",
  "annotated_at": "2026-04-10T12:00:00+00:00",
  "annotation_tool": "variant_annotator/1.2",
  "schema_version": "1.2",
  "gff3_reference": "/ref/gencode.v49.basic.annotation.gff3",
  "fasta_reference": "/fasta/hg38.fa",
  "dbnsfp_reference": "/dbnsfp/dbNSFP5.0a_grch38.gz",
  "regulatory_beds": [],
  "betamap_input": "/data/PGS000001/...",
  "elapsed_seconds": 185.3,
  "stats": {
    "total_variants":  12345,
    "n_coding":        234,
    "n_regulatory":    89,
    "n_intergenic":    4567,
    "n_genic":         7778,
    "n_splice_site":   3,
    "n_splice_donor":  7,
    "n_splice_acceptor": 5,
    "n_splice_region": 45,
    "n_missense":      89,
    "n_synonymous":    62,
    "n_lof":           19,
    "n_stop_gained":   4,
    "n_frameshift":    3,
    "pct_coding":      1.90,
    "pct_intergenic":  37.0,
    "pct_missense":    0.72,
    "pct_lof":         0.15,
    ...
  },
  "resolved_limitations": [...],
  "remaining_limitations": [...]
}
```

---

## Region Classes

| region_class    | GFF3 feature driving it            | Consequence                        | is_coding |
|-----------------|------------------------------------|------------------------------------|-----------|
| coding          | CDS                                | missense_variant \| synonymous_variant \| stop_gained \| stop_lost \| start_lost \| frameshift_variant \| inframe_insertion \| inframe_deletion \| coding_sequence_variant (no FASTA) | true |
| UTR_5prime      | five_prime_UTR                     | 5_prime_UTR_variant                | false     |
| UTR_3prime      | three_prime_UTR                    | 3_prime_UTR_variant                | false     |
| non_coding_exon | exon (non-CDS)                     | non_coding_exon_variant            | false     |
| intronic        | transcript or gene                 | intron_variant                     | false     |
| intronic        | intronic, distance ≤ 2 bp          | splice_donor_variant \| splice_acceptor_variant (FASTA) \| splice_site_variant (GFF3-only fallback) | false |
| splice_region   | intronic, distance 3–8 bp          | splice_region_variant              | false     |
| near_gene       | none (within 2 kb, no BED)         | upstream_gene_variant              | false     |
| regulatory      | none (overlaps BED element)        | regulatory_region_variant          | false     |
| intergenic      | none                               | intergenic_variant                 | false     |

---

## Consequence Priority

The `consequence_priority` column allows sorting variants by biological impact (lower = higher impact):

| consequence                 | priority |
|-----------------------------|----------|
| stop_gained                 | 50       |
| frameshift_variant          | 60       |
| stop_lost                   | 70       |
| start_lost                  | 80       |
| splice_donor_variant        | 100      |
| splice_acceptor_variant     | 110      |
| splice_site_variant         | 120      |
| inframe_insertion           | 150      |
| inframe_deletion            | 160      |
| splice_region_variant       | 200      |
| missense_variant            | 250      |
| synonymous_variant          | 280      |
| coding_sequence_variant     | 300      |
| 5_prime_UTR_variant         | 400      |
| 3_prime_UTR_variant         | 500      |
| non_coding_exon_variant     | 600      |
| intron_variant              | 700      |
| upstream_gene_variant       | 800      |
| regulatory_region_variant   | 850      |
| intergenic_variant          | 900      |

---

## Coding Consequence Logic (v1.2, requires `--fasta`)

For CDS variants, `FASTAEngine` computes the consequence as follows:

1. **REF/ALT determination** — the FASTA is queried at the variant position to identify which of `EFFECT_ALLELE` / `OTHER_ALLELE` matches the reference.
2. **Indel classification** — if `len(REF) != len(ALT)`: classify as `frameshift_variant` (len diff not divisible by 3), `inframe_insertion`, or `inframe_deletion`.
3. **SNP codon extraction** — `GFF3Index.get_cds_offset()` returns the 0-based CDS position; `FASTAEngine._get_codon_seq()` fetches the 3-nt reference codon from the FASTA, handling multi-exon codons.
4. **Strand orientation** — on minus-strand genes, `EFFECT_ALLELE` (always on + strand) is complemented before codon substitution.
5. **Translation and classification** — `snp_codon_consequence()` classifies: `start_lost`, `stop_lost`, `stop_gained`, `synonymous_variant`, `missense_variant`.

---

## Splice Site Classification (v1.2, requires `--fasta`)

For intronic variants within ≤ 2 bp of an exon boundary:

- `FASTAEngine.get_splice_type()` finds the nearest exon boundary using the IntervalTree.
- Reads the canonical dinucleotide at the boundary from the FASTA:
  - **+ strand, exon end boundary:** reads `FASTA[iv.end : iv.end+2]` → `GT` = donor
  - **+ strand, exon start boundary:** reads `FASTA[iv.begin-2 : iv.begin]` → `AG` = acceptor
  - **− strand:** same logic with reverse complement applied
- Returns `"donor"` → `splice_donor_variant` or `"acceptor"` → `splice_acceptor_variant`.
- Falls back to `splice_site_variant` if FASTA is unavailable or the dinucleotide is non-canonical.

---

## Multi-Gene Overlap Model

When a variant overlaps multiple genes:

1. All overlapping genes are enumerated from the GFF3 interval tree.
2. For each gene, the most specific feature (CDS > UTR > exon > transcript > gene) is selected.
3. The gene with the highest specificity (protein_coding preferred at ties) becomes the **best hit**.
4. All overlapping genes appear in `all_overlapping_genes`, corresponding transcripts in `all_overlapping_transcripts`, and corresponding region classes in `all_region_classes`.
5. `n_overlapping_genes` counts all distinct gene IDs.

---

## Splice Site Distance

`distance_to_splice_site` is populated for every variant (not only intronic).

- The GFF3 exon interval tree is queried in an expanding window (200 → 2 000 → 20 000 bp).
- Distance = `min(|pos − exon_start|, |pos − exon_end|)` across candidate exons.
- Thresholds: ≤ 2 bp → canonical splice site; 3–8 bp → splice region; > 8 bp → unchanged.

---

## Regulatory Annotation

### With `--regulatory-bed` (element-level)

Variants overlapping any BED element receive `is_regulatory = True`, `regulatory_source = "bed:<name>"`.

### Without `--regulatory-bed` (heuristic fallback)

Variants within ±2000 bp of a gene boundary receive `is_regulatory = True`, `regulatory_source = "heuristic"`.

---

## Web Endpoints

| Endpoint                                   | Description                              |
|--------------------------------------------|------------------------------------------|
| `GET /variants/<pgs_id>`                   | HTML annotation view                     |
| `GET /api/variants/<pgs_id>`               | Paginated variant table JSON             |
| `GET /api/variants/<pgs_id>/summary`       | Full annotation summary JSON             |
| `GET /api/variants/<pgs_id>/ideogram`      | Chromosome-level data for ideogram       |
| `GET /api/variants/<pgs_id>/run-info`      | Apptainer run command for this PGS_ID    |

### Paginated variant API query parameters

| Parameter    | Default | Description                                                   |
|--------------|---------|---------------------------------------------------------------|
| page         | 1       | Page number                                                   |
| page_size    | 500     | Rows per page (max 5000)                                      |
| chrom        | —       | Filter by chromosome (e.g. `1`, `X`, `MT`)                   |
| region_class | —       | Filter by region class. `splice_region` matches both `splice_site_variant` and `splice_region_variant`. |
| gene_name    | —       | Match best-hit gene or any gene in `all_overlapping_genes`    |
| only_coding  | false   | If `1`/`true`, return only coding variants                    |

---

## Configuration

```
APP_ANNOTATIONS_DIR=/annotations
```

Default: `/annotations`

---

## Resolved vs Remaining Limitations

### Resolved in v1.1

| Limitation (v1.0)                         | Resolution                                   |
|-------------------------------------------|----------------------------------------------|
| Splice-site variants not detected         | `distance_to_splice_site` computed for all variants; intronic variants near exon boundaries reclassified. |
| Only best overlapping gene reported       | `all_overlapping_genes`, `all_overlapping_transcripts`, `all_region_classes` columns. |
| Regulatory is TSS-only heuristic          | `--regulatory-bed` adds element-level annotation from any BED file. |
| No consequence priority ordering          | `consequence_priority` numeric column added. |

### Resolved in v1.2

| Limitation (v1.1)                                 | Resolution                                   |
|---------------------------------------------------|----------------------------------------------|
| Missense/synonymous/stop-gained not computed      | `FASTAEngine.get_codon_consequence()` computes full codon-level consequence when `--fasta` is provided. |
| Splice donor vs acceptor not distinguished        | `FASTAEngine.get_splice_type()` reads GT/AG from FASTA; classifies `splice_donor_variant` / `splice_acceptor_variant`. |
| No functional scores                              | `DbNSFP5Fetcher` queries CADD, REVEL, SIFT, PolyPhen2, ClinVar from dbNSFP5 tabix when `--dbnsfp` is provided. |
| REF allele not verified                           | `FASTAEngine.determine_ref_alt()` confirms REF by FASTA lookup. |

### Optional / configurable

| Capability                                | Requirement                                  |
|-------------------------------------------|----------------------------------------------|
| Coding consequences (missense, etc.)      | `FASTA_PATH=/path/to/hg38.fa` |
| dbNSFP5 functional scores                 | `DBNSFP_PATH=/path/to/dbNSFP5.0a_grch38.gz` |
| Element-level regulatory annotation      | `REGULATORY_BED=/path/to/encode.bed.gz ...` |

---

## Troubleshooting

**GFF3 index is empty:**
Check that the GFF3 file path is correct inside the container (`/ref/...`) and that the `--bind` mount includes the GFF3 directory.

**FASTA not found / codon consequences absent:**
Verify `FASTA_PATH` is set, the file exists, and a `.fai` index is present in the same directory.

**dbNSFP5 not found:**
Verify `DBNSFP_PATH` is set, the `.gz` file exists, and the `.gz.tbi` tabix index is present alongside it.

**`codon_ref` / `aa_ref` columns empty for CDS variants:**
The CDS variant's transcript_id may not be in `transcript_cds_map`. Check that the GFF3 includes CDS features (GENCODE basic always does).

**distance_to_splice_site is always None:**
The GFF3 has no `exon` features for that chromosome. Verify the GFF3 file includes exon lines.

**Web shows "Not annotated" after pipeline completes:**
Verify `APP_ANNOTATIONS_DIR` on the web container matches the `ANNOTATIONS_DIR` used during the run.

**Old annotation files missing v1.2 columns:**
The service gracefully handles absent columns — they show as null/empty in the table.
Re-run the pipeline with the new container and `FASTA_PATH` + `DBNSFP_PATH` to get full v1.2 output.

**High memory usage:**
The GFF3 index for all 25 chromosomes uses ~6–10 GB. Request at least 64 GB on the HPC node when FASTA + dbNSFP5 are active.
