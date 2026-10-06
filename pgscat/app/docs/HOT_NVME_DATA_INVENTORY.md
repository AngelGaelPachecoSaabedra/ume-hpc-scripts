# HOT_NVME Data Inventory
**Generated:** 2026-04-27  
**Source root:** `/mnt/cephfs/hot_nvme/`

---

## Classification Key

| Label | Meaning |
|-------|---------|
| `SAFE_PUBLIC_DERIVED` | Publicly available data (gnomAD, ClinVar, dbSNP, GENCODE) — safe to serve |
| `INTERNAL_AGGREGATED_ONLY` | MCPS aggregate statistics — no individual-level data; internal use only |
| `SENSITIVE_DO_NOT_EXPOSE` | Individual-level data (CRAMs, GVCFs, pVCFs, phenotypes) — never expose |
| `UNKNOWN_NEEDS_REVIEW` | Not fully reviewed |

---

## gnomAD 4.1.1 — `SAFE_PUBLIC_DERIVED`

**Path:** `/mnt/cephfs/hot_nvme/gnomad_4.1.1/`  
**Status:** ✅ INTEGRATED (v1.4, live-enriched)  
**License:** CC0 (public domain)

| Dataset | Path | Chromosomes | Index |
|---------|------|-------------|-------|
| Exome sites | `variants/exome/gnomad.exomes.v4.1.1.sites.chrN.vcf.bgz` | 1-22, X, Y | `.tbi` ✅ |
| Genome sites | `variants/genome/gnomad.genomes.v4.1.1.sites.chrN.vcf.bgz` | 1-22, X, Y | `.tbi` ✅ |
| SVs | `variants/sv/gnomad.v4.1.sv.sites.vcf.gz` | all | `.tbi` ✅ |

**Integrated fields:**
- `af_gnomad` — global AF (exome preferred, genome fallback)
- `af_gnomad_afr`, `af_gnomad_amr`, `af_gnomad_eas`, `af_gnomad_nfe`, `af_gnomad_sas`, `af_gnomad_mid`
- `an_gnomad`, `nhomalt_gnomad`, `rarity_class_gnomad`

**Integration points:**
- `src/services/gnomad_fetcher.py` — runtime pysam tabix fetcher
- `annotator/gnomad_freq.py` — annotator pipeline module
- `src/services/variant_annotation.py` — live enrichment of variant tables
- `src/services/gene_browser.py` — AF track uses gnomAD preferentially
- `src/services/variant_ranking.py` — rarity score uses gnomAD AF

**Notes:**
- AN > 1M in exome — highest quality AF source available
- AMR population = Admixed American = best gnomAD proxy for Mexican/Latin American
- SV dataset not yet integrated (future work)

---

## ClinVar GRCh38 — `SAFE_PUBLIC_DERIVED`

**Path:** `/mnt/cephfs/hot_nvme/clinvar/GRCh38/clinvar.vcf.gz`  
**Updated:** 2026-04-15  
**Status:** ✅ INTEGRATED (derived Parquet index built)  
**License:** Public domain (NCBI)

**Files:**
- `clinvar.vcf.gz` (BGZF) — 4,404,292 variants, no .tbi (built derived Parquet)
- Derived: `/mnt/cephfs/hot/pgscat/work/clinvar_idx/clinvar_grch38.parquet` (94MB, Snappy)

**Integrated fields:**
- `clnsig` / `clinvar_clnsig` — germline pathogenicity classification
- `clndn` / `clinvar_clndn` — ClinVar preferred disease name
- `clnrevstat` / `clinvar_clnrevstat` — review status (evidence level)
- `alleleid` / `clinvar_alleleid` — ClinVar Allele ID
- `clnvc` / `clinvar_clnvc` — variant type (SNV, Deletion, etc.)

**Integration points:**
- `src/services/clinvar_fetcher.py` — DuckDB-backed runtime fetcher
- `annotator/clinvar_anno.py` — annotator in-memory index

**Notes:**
- Chromosome naming: bare numbers (1, 2, ..., X, Y — no chr prefix)
- Build time: ~31 seconds
- Replaces/enriches the dbNSFP5-derived `clinvar_clnsig` with authoritative VCF source
- To rebuild: `python scripts/build_external_indexes.py --clinvar`

---

## MCPS Variant Browser AFs — `INTERNAL_AGGREGATED_ONLY`

**Path:** `/mnt/cephfs/hot_nvme/mcps/mcps-variant-browser-afs/`  
**Status:** ✅ INTEGRATED (DuckDB index building in background ~30 min)  
**Privacy:** Aggregated population frequencies only — no individual data

**Files:** `chr{1-22,X}.freq.tsv.gz` — 24 files, ~2.2GB total compressed  
~10.9M variants on chr1 alone (est. ~80-100M total across all chromosomes)  
Derived: `/mnt/cephfs/hot/pgscat/work/mcps_idx/mcps_af.duckdb`

**Format columns:**
```
ID  CPRA  SOURCE  STON  GHOM  AHOM
AN_RAW  AC_RAW  AF_RAW
AN_AFR  AC_AFR  AF_AFR
AN_EUR  AC_EUR  AF_EUR
AN_MEX  AC_MEX  AF_MEX
```

**Integrated fields:**
- `af_mcps` (= AF_RAW) — global MCPS allele frequency
- `af_mcps_afr` — African ancestry AF
- `af_mcps_eur` — European ancestry AF
- `af_mcps_mex` — Mexican ancestry AF (**key population for this cohort**)

**Integration points:**
- `src/services/mcps_fetcher.py` — DuckDB-backed runtime fetcher
- `src/services/variant_annotation.py` — live enrichment
- `src/services/gene_browser.py` — AF hover labels

**Notes:**
- Population `MEX` = ~185K alleles/site → statistically robust Mexican AF
- Source column: `WESOnly` | `WGSOnly` | `WES+WGS`
- To rebuild: `python scripts/build_external_indexes.py --mcps`
- Build time: ~15-30 minutes for all chromosomes

---

## dbSNP Population Frequencies — `SAFE_PUBLIC_DERIVED`

**Path:** `/mnt/cephfs/hot_nvme/dbsnp/population_frequency/`  
**Status:** ✅ Already integrated in annotator (v1.3 `dbsnp_freq.py`)  

**Files:** `freq.vcf.gz` + `.tbi` — NCBI GRAF-pop format  
`population_map.expanded.tsv`, `contigs.txt`

**Notes:**
- Used by `annotator/dbsnp_freq.py` during pipeline annotation
- Produces: rsid, af_global, af_max_population, af_population_summary, rarity_class
- gnomAD 4.1.1 now preferred for AF quality — dbSNP retained as fallback

---

## dbNSFP 5.0a — `SAFE_PUBLIC_DERIVED`

**Path:** `/mnt/cephfs/hot_nvme/dbNSFP5/`  
**Status:** ✅ Already integrated in annotator (v1.2 `dbnsfp.py`)  

**Files:**
- `dbNSFP5.0a_grch38.gz` + `.tbi` — combined lookup file
- `dbNSFP5.0a_variant.chrN.gz` — per-chromosome variant files

**Provides:** CADD phred, REVEL, SIFT, PolyPhen2, clinvar_clnsig (from dbNSFP)

---

---

## Regulatory Annotation — `SAFE_PUBLIC_DERIVED`

**Status:** ✅ INTEGRATED (v1.5 — real BED annotation replaces ±2 kb heuristic)  
**Warning:** UCSC data pre-dates ENCODE phase 3 cCREs. Data from 2015-2022.  
**Build:** GRCh38/hg38 confirmed (chr-prefix + GRCh38 alt scaffolds present in data)

### ENCODE TFBS Clusters

**Source:** `/mnt/cephfs/hot_nvme/ucsc/encRegTfbsClustered.txt.gz`  
**Size:** 196 MB (UCSC table format, bin column leading)  
**Derived:** `/mnt/cephfs/hot/pgscat/work/regulatory_derived/encode_tfbs_GRCh38.bed.gz`  
**Format:** `bin chrom chromStart chromEnd TF_name score nCellTypes ...`  
**Content:** ENCODE project transcription factor binding site clusters across cell types.  
**Utility:** ★★★★★ — Gold standard for TFBS. TF name (col 4) parsed to infer element type:
- `CTCF`, `CTCFL`, `RAD21`, `SMC3` → **insulator** (TAD boundaries)
- `H3K4me3`, `H3K9ac` → **promoter**
- `H3K27ac`, `H3K4me1`, `EP300` → **enhancer**
- Other named TFs → **TFBS**

**Recommendation:** ✅ USE — primary regulatory annotation source

### ENCODE DNase I Hypersensitivity Site (DHS) Clusters

**Source:** `/mnt/cephfs/hot_nvme/ucsc/wgEncodeRegDnaseClustered.txt.gz`  
**Size:** 71 MB (UCSC table format)  
**Derived:** `/mnt/cephfs/hot/pgscat/work/regulatory_derived/encode_dhs_GRCh38.bed.gz`  
**Format:** `bin chrom chromStart chromEnd nCellTypes maxScore ...`  
**Content:** DNase I hypersensitivity clusters across cell types — marks open chromatin / active regulatory regions.  
**Utility:** ★★★★☆ — Complements TFBS; identifies accessible chromatin without TF identity.  
**Element type:** `open_chromatin`

**Recommendation:** ✅ USE — secondary regulatory source (open chromatin)

### ORegAnno Curated Regulatory Elements

**Source:** `/mnt/cephfs/hot_nvme/ucsc/oreganno.txt.gz`  
**Size:** 20 MB (UCSC table format)  
**Derived:** `/mnt/cephfs/hot/pgscat/work/regulatory_derived/oreganno_GRCh38.bed.gz`  
**Format:** `bin chrom chromStart chromEnd OREG_ID strand OREG_ID`  
**Content:** ORegAnno database of curated regulatory regions (promoters, enhancers, etc.).  
**Utility:** ★★★☆☆ — Literature-curated; smaller than ENCODE but manually annotated.  
**Element type:** `regulatory_region` (IDs like OREG1234567)

**Recommendation:** ✅ USE — tertiary source for curated elements

### Individual ENCODE TF ChIP-seq Experiments (1,256 files)

**Source:** `/mnt/cephfs/hot_nvme/ucsc/encTfChipPkENCFF*.txt.gz`  
**Size:** ~400 KB–2 MB each  
**Content:** Individual ENCODE ChIP-seq peak calls per TF per cell type (ENCFF accession IDs).  
**Build:** **Unknown — may be hg19 or hg38** (UCSC mixes builds for these; verify before use).  
**Utility:** ★★★☆☆ — More granular than clusters, but requires build verification per file.  
**Not integrated** — Use `encRegTfbsClustered.txt.gz` (clustered version) instead.

### Genomic Build Confirmation

**Build confirmation method:** The presence of `chr10_KN196480v1_fix` (a GRCh38-exclusive alt scaffold)  
in `encRegTfbsClustered.txt.gz` and `wgEncodeRegDnaseClustered.txt.gz` **confirms GRCh38**.  
GRCh37/hg19 does not contain `KN*` alt scaffolds.

⚠️ The individual `encTfChipPkENCFF*.txt.gz` files (1,256 ChIP-seq experiments) have  
**unverified build** — check file content before using. The clustered files above are GRCh38.

### Regulatory Data NOT Found in hot_nvme

The following state-of-the-art sources were **not found** and would further improve annotation:
- **ENCODE cCRE V3 (2023)** — `encodeCcreCombined.txt.gz` or `ENCFF833FGR.bed.gz` — current gold standard with 5 element types (PLS/pELS/dELS/CTCF-only/DNase-H3K4me3). **Recommend downloading.**
- **Ensembl Regulatory Build** — `homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz` — ensemble of many cell types
- **ATAC-seq files** — not found in hot_nvme
- **ChromHMM chromatin state segmentation** — not found in hot_nvme

### Integration Points (v1.5)

| File | Source | Build | Status |
|------|--------|-------|--------|
| `encode_tfbs_GRCh38.bed.gz` | `ucsc/encRegTfbsClustered.txt.gz` | GRCh38 | Build with script |
| `encode_dhs_GRCh38.bed.gz` | `ucsc/wgEncodeRegDnaseClustered.txt.gz` | GRCh38 | Build with script |
| `oreganno_GRCh38.bed.gz` | `ucsc/oreganno.txt.gz` | GRCh38 | Build with script |

**New annotation columns (v1.5):**
- `regulatory_element_type` — promoter | enhancer | insulator | TFBS | open_chromatin | regulatory_region
- `regulatory_element_id` — TF name (CTCF, GATA1, …) or OREG ID
- `distance_to_regulatory` — 0 if overlapping; bp to nearest element if within 100 kb; None if no BED

**Build command:**
```bash
cd /mnt/cephfs/hot/pgscat/app
python scripts/build_regulatory_beds.py --all
# Estimated: ~6-7 min total (196+71+20 MB input)
# Output: /mnt/cephfs/hot/pgscat/work/regulatory_derived/
```

**Runtime integration:**
```bash
export REGULATORY_BED="/mnt/cephfs/hot/pgscat/work/regulatory_derived/encode_tfbs_GRCh38.bed.gz,\
/mnt/cephfs/hot/pgscat/work/regulatory_derived/encode_dhs_GRCh38.bed.gz,\
/mnt/cephfs/hot/pgscat/work/regulatory_derived/oreganno_GRCh38.bed.gz"
```
`run_annotation.sh` auto-detects these files when built — no manual export needed.

**Ranking bonus:** Non-coding variants in regulatory elements now receive a consequence score boost:
- promoter → 0.80, insulator → 0.70, enhancer → 0.65, open_chromatin → 0.40, TFBS → 0.35

---

## MCPS Local Ancestry — `SENSITIVE_DO_NOT_EXPOSE`

**Path:** `/mnt/cephfs/hot_nvme/mcps/local-ancestry/rfmix/3-way/`  
**Status:** ⚠️ NOT integrated — individual-level ancestry data  

Contains RFMix local ancestry calls (3-way: European/African/Native American).  
**Individual-level data — NEVER expose in any API or UI.**  
Only aggregate proportions (e.g., cohort-level admixture fractions) would be safe.

---

## MCPS Phenotypes — `SENSITIVE_DO_NOT_EXPOSE`

**Path:** `/mnt/cephfs/hot_nvme/mcps/phenotypes/`  
**Status:** ⚠️ NOT integrated — individual phenotype data  

Files: BASELINE.csv, BASE_LAB.csv, BASE_DRUGS.csv, BASE_NMR.csv, DEATHS.csv, RGN_LINK.csv, data-dictionaries/  
**Individual-level phenotype data — NEVER expose.**

---

## MCPS Whole Exome Sequencing — `SENSITIVE_DO_NOT_EXPOSE`

**Path:** `/mnt/cephfs/hot_nvme/mcps/whole_exome_sequencing/`  
**Status:** ⚠️ NOT integrated  

Contains individual CRAMs, GVCFs, pVCFs.  
Only the manifest/QC summary is safe to reference:  
- `2-11-21 - MCPS Exome Sequencing Update (F145K).xlsx` — QC summary (145K samples)

---

## MCPS Whole Genome Sequencing — `SENSITIVE_DO_NOT_EXPOSE`

**Path:** `/mnt/cephfs/hot_nvme/mcps/whole_genome_sequencing/`  
**Status:** ⚠️ NOT integrated  

Contains individual CRAM/GVCF/pVCF files. Freeze Two manifest available.  
**Individual-level genomic data — NEVER expose.**

**Safe metadata only:**
- `dataStructure.README.MCPS_WGS_Freeze_Two.pVCF.txt`
- `MCPS_WGS_Freeze_Two_pVCF.manifest.md5_summary.tsv`

---

## MCPS PGS — `INTERNAL_AGGREGATED_ONLY`

**Path:** `/mnt/cephfs/hot_nvme/mcps/mcps-pgs/`  
**Status:** `UNKNOWN_NEEDS_REVIEW`  
Requires review before integration.

---

## MCPS Imputed TOPMed — `SENSITIVE_DO_NOT_EXPOSE`

**Path:** `/mnt/cephfs/hot_nvme/mcps/imputed-topmed/`  
**Status:** ⚠️ NOT integrated — individual-level imputed genotype data

---

## MCPS Genotyping — `SENSITIVE_DO_NOT_EXPOSE`

**Path:** `/mnt/cephfs/hot_nvme/mcps/genotyping/`  
**Status:** ⚠️ NOT integrated — individual-level array genotype data

---

## Ancestry Reference — `SAFE_PUBLIC_DERIVED`

**Path:** `/mnt/cephfs/hot_nvme/ancestry/`  
**Status:** Not integrated (contains 1KGP + MCPS exome reference)  

Files: `ALL_chr.MCPS-exome_1KGP_noSAS.GRCh38.vcf.gz` — 1000 Genomes + MCPS exome reference panel.  
Could be used for ancestry PCA or population stratification analyses.

---

## GENCODE GRCh38.14 — `SAFE_PUBLIC_DERIVED`

**Path:** `/mnt/cephfs/hot_nvme/gencode/GRCh38.14/`  
**Status:** ✅ Already used by annotator pipeline (GFF3 annotation)

---

## hg38 Reference — `SAFE_PUBLIC_DERIVED`

**Path:** `/mnt/cephfs/hot_nvme/hg38/`  
**Status:** ✅ Already used by annotator pipeline (FASTA engine)  

---

## Integration Summary

| Source | Integration | Fields Added | Index Type |
|--------|------------|--------------|-----------|
| gnomAD 4.1.1 Exome | ✅ Live (pysam) | af_gnomad + 8 pop AFs | tabix BGZ |
| gnomAD 4.1.1 Genome | ✅ Fallback | same | tabix BGZ |
| ClinVar GRCh38 | ✅ Live (DuckDB) | clnsig, clndn, clnrevstat | Parquet 94MB |
| MCPS AF | ✅ Live (DuckDB) | af_mcps + MEX/AFR/EUR | DuckDB ~1.5GB |
| dbSNP pop freq | ✅ Annotator | rsid, af_global, rarity_class | tabix BGZ |
| dbNSFP5.0a | ✅ Annotator | CADD, REVEL, SIFT, PolyPhen2 | tabix |
| GENCODE GRCh38.14 | ✅ Annotator | gene/transcript annotation | GFF3 index |
| hg38 FASTA | ✅ Annotator | codon/aa changes | FASTA |
| ENCODE TFBS (UCSC) | ✅ Annotator v1.5 | regulatory_element_type, _id, distance | BED.GZ (derived) |
| ENCODE DHS (UCSC) | ✅ Annotator v1.5 | open_chromatin annotation | BED.GZ (derived) |
| ORegAnno (UCSC) | ✅ Annotator v1.5 | curated regulatory regions | BED.GZ (derived) |

---

## Privacy Warnings

1. **NEVER** expose individual-level MCPS data (CRAMs, GVCFs, phenotypes, local ancestry).
2. MCPS aggregate AFs (mcps-variant-browser-afs) are safe — these are population-level statistics.
3. The `af_mcps`, `af_mcps_mex` fields served by the API contain no individual-level information.
4. All gnomAD and ClinVar data is public-domain.

---

## Commands to Rebuild Indexes

```bash
cd /mnt/cephfs/hot/pgscat/app

# ClinVar Parquet (~31 seconds)
python scripts/build_external_indexes.py --clinvar

# MCPS DuckDB (~15-30 minutes)
python scripts/build_external_indexes.py --mcps

# Both
python scripts/build_external_indexes.py --all

# Via API (admin endpoint, synchronous):
curl -X POST http://localhost:8080/api/admin/build-indexes \
     -H 'Content-Type: application/json' \
     -d '{"sources": ["clinvar"]}'
```

---

## Is It Necessary to Re-run the Annotator?

**For gnomAD + ClinVar + MCPS:** NO — these are live-enriched at query time from the new fetchers. Existing annotation parquets are automatically enriched when served via `/api/variants/<pgs_id>` and `/api/gene/<gene>/tracks`.

**When to re-run the annotator:**
- When you want gnomAD/ClinVar/MCPS data baked INTO the parquets (faster reads, no live fetch overhead)
- For this, add `--gnomad-dir`, `--clinvar`, `--mcps-dir` flags to `annotate_variants.py` (see `annotator/gnomad_freq.py` and `annotator/clinvar_anno.py`)

**Annotator re-run command (HPC/Apptainer):**
```bash
apptainer exec variant_annotator.sif python annotate_variants.py \
    --betamap /data/PGS000001/PGS000001_hmPOS_GRCh38.betamap.tsv.gz \
    --gff3    /ref/gencode.v49.basic.annotation.gff3 \
    --outdir  /annotations/PGS000001 \
    --pgs-id  PGS000001 \
    --dbsnp   /mnt/cephfs/hot_nvme/dbsnp/population_frequency/freq.vcf.gz \
    --dbnsfp  /mnt/cephfs/hot_nvme/dbNSFP5/dbNSFP5.0a_grch38.gz \
    # Future flags (not yet wired to CLI):
    # --gnomad-dir /mnt/cephfs/hot_nvme/gnomad_4.1.1/variants
    # --clinvar    /mnt/cephfs/hot_nvme/clinvar/GRCh38/clinvar.vcf.gz
```
