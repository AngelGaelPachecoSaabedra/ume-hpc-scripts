#!/usr/bin/env python3
"""
diagnose_pgs.py — Quick column/data diagnostic for an annotated PGS parquet.

Usage:
    python3 scripts/diagnose_pgs.py PGS000004
    python3 scripts/diagnose_pgs.py PGS000004 --annotations-dir /annotations
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq


def main():
    ap = argparse.ArgumentParser(description="Diagnose annotated PGS parquet columns.")
    ap.add_argument("pgs_id", help="PGS ID, e.g. PGS000004")
    ap.add_argument(
        "--annotations-dir",
        default="/mnt/cephfs/hot/pgscat/annotations",
        help="Root annotations directory (default: /mnt/cephfs/hot/pgscat/annotations)",
    )
    args = ap.parse_args()

    pgs_id = args.pgs_id.strip().upper()
    ann_dir = Path(args.annotations_dir)
    pq_path = ann_dir / pgs_id / f"{pgs_id}_variants_annotated.parquet"

    if not pq_path.exists():
        print(f"[ERROR] Parquet not found: {pq_path}", file=sys.stderr)
        sys.exit(1)

    print(f"\n{'='*60}")
    print(f"  Diagnostic: {pgs_id}")
    print(f"  File: {pq_path}")
    print(f"{'='*60}\n")

    df = pq.read_table(pq_path).to_pandas()
    print(f"Shape: {df.shape[0]} rows × {df.shape[1]} columns\n")

    GROUPS = {
        "Consequence / Region": ["consequence", "region_class", "is_coding", "is_regulatory",
                                  "is_lof", "is_missense", "is_synonymous"],
        "CADD / Pathogenicity": ["cadd_phred", "revel_score", "sift_pred", "polyphen2_pred"],
        "ClinVar":              ["clinvar_clnsig", "clinvar_clndn", "clinvar_clnrevstat", "clinvar_alleleid"],
        "dbSNP / Global AF":   ["rsid", "af_global", "af_max_population", "af_population_summary",
                                  "rarity_class"],
        "gnomAD AF (live)":    ["af_gnomad", "af_gnomad_afr", "af_gnomad_amr", "af_gnomad_eas",
                                  "af_gnomad_nfe", "af_gnomad_sas", "rarity_class_gnomad"],
        "MCPS AF (live)":      ["af_mcps", "af_mcps_afr", "af_mcps_eur", "af_mcps_mex",
                                  "an_mcps", "ac_mcps"],
        "Effective AF (v1.5)": ["af_effective", "af_source", "rarity_class_effective"],
        "Regulatory":          ["regulatory_element_type", "regulatory_element_id",
                                  "regulatory_source", "distance_to_regulatory"],
        "PRS / Ranking":       ["BETA", "IS_FLIP"],
    }

    for group, cols in GROUPS.items():
        print(f"── {group}")
        for col in cols:
            if col not in df.columns:
                print(f"   {col:<35}  [ABSENT]")
                continue
            s     = df[col]
            dtype = str(s.dtype)
            n_null = s.isna().sum()
            n_ok   = len(s) - n_null
            pct    = 100 * n_ok / len(s) if len(s) else 0
            # Sample non-null values
            sample = s.dropna().head(3).tolist()
            sample_str = str(sample)[:60]
            print(f"   {col:<35}  dtype={dtype:<10}  non-null={n_ok}/{len(s)} ({pct:.0f}%)  sample={sample_str}")
        print()

    # Summary flags
    print("── Summary flags")
    flags = {
        "CADD column present":     "cadd_phred" in df.columns,
        "CADD has values":         "cadd_phred" in df.columns and df["cadd_phred"].notna().any(),
        "REVEL has values":        "revel_score" in df.columns and df["revel_score"].notna().any(),
        "gnomAD AF available":     "af_gnomad" in df.columns and df["af_gnomad"].notna().any(),
        "MCPS AF available":       "af_mcps" in df.columns and df["af_mcps"].notna().any(),
        "Effective AF available":  "af_effective" in df.columns and df["af_effective"].notna().any(),
        "ClinVar data available":  "clinvar_clnsig" in df.columns and df["clinvar_clnsig"].notna().any(),
        "Regulatory BED annotated":"regulatory_element_type" in df.columns and df["regulatory_element_type"].notna().any(),
        "BETA present":            "BETA" in df.columns and df["BETA"].notna().any(),
    }
    for label, ok in flags.items():
        status = "OK " if ok else "---"
        print(f"   [{status}]  {label}")

    print()
    if "BETA" in df.columns:
        print(f"   Total Σ|BETA|: {df['BETA'].abs().sum():.6f}")
        print(f"   Max  |BETA|:   {df['BETA'].abs().max():.6f}")
        print(f"   Median |BETA|: {df['BETA'].abs().median():.6f}")

    print(f"\n{'='*60}\n")


if __name__ == "__main__":
    main()
