#!/usr/bin/env python3
"""
compute_prs_spark_gpu.py
========================
PRS por cromosoma usando PySpark para matching de variantes y CuPy para
scoring vectorizado en GPU.

CONTRATO DE SALIDA: idéntico a compute_prs.py (interfaz externa congelada).
  {OUT_BASE}/{pgs_id}/{pgs_id}_chr{chrom}_scores.tsv   — sample_id, PRS
  {OUT_BASE}/{pgs_id}/{pgs_id}_chr{chrom}_metadata.json
  OUT_BASE = /mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/scripts/resultados/

FUENTES (solo lectura, nunca se escribe en hot_nvme):
  Betamap : {SCORES_BASE}/{pgs_id}/{pgs_id}_hmPOS_GRCh38.betamap.tsv.gz
  Zarr    : {ZARR_BASE}/chr{chrom}.zarr  (zar_files/ por cromosoma)

Modelo aditivo:
    PRS = Σ_i  BETA_i × dosage_i
    dosage_i = call_DS (directo)       si IS_FLIP=0  (efecto=ALT)
    dosage_i = 2 − call_DS (invertido) si IS_FLIP=1  (efecto=REF)

Arquitectura en tres fases:
    Fase 1 – Spark  : broadcast join catálogo-zarr vs betamap;
                      resuelve IS_FLIP, strand-flip y mismatches de alelos.
    Fase 2 – CuPy   : scoring GPU por ventanas zarr; solo transfiere a VRAM
                      las filas matched (n_match << window_size), no la
                      ventana entera. Peak VRAM ≈ n_match × n_samples × 4 B.
    Fase 3 – Salida : TSV + JSON con el mismo esquema que compute_prs.py.

Integración con la infraestructura HPC:
    Spark:   /mnt/cephfs/biocontainers/spark/spark  (instalación real, sin pip)
    PySpark: importado vía PYTHONPATH=$SPARK_HOME/python (configurado en el sbatch)
    Java:    /usr/lib/jvm/java-17-openjdk-amd64 (expuesto al contenedor vía bind)
    Modo:    local[N] para nodo único GPU; spark://IP:7077 para cluster multi-nodo

Uso:
    python compute_prs_spark_gpu.py --pgs-id PGS000001 --chrom 1
    python compute_prs_spark_gpu.py --pgs-id PGS000001 --chrom X \\
        --window-size 20000 --missing-strategy mean \\
        --spark-master local[16]
"""

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import zarr

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stdout,
)
log = logging.getLogger(__name__)

# ── Rutas de ENTRADA (hot_nvme — solo lectura, nunca se escribe aquí) ────────
# Betamap pre-procesado del PGS Catalog
SCORES_BASE = Path("/mnt/cephfs/hot_nvme/pgscatalog/scores")
# Zarr por cromosoma (chr1.zarr … chrX.zarr) — validados en el pipeline previo
ZARR_BASE   = Path("/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files")

# ── Ruta de SALIDA (home del usuario — única zona de escritura) ───────────────
# TSV + JSON se guardan aquí; nunca en hot_nvme
OUT_BASE    = Path("/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/scripts/resultados")

# Complemento de bases para strand-flip
_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="PRS por cromosoma: Spark matching + CuPy GPU scoring."
    )
    p.add_argument("--pgs-id",  required=True, help="PGS Catalog ID, ej. PGS000001")
    p.add_argument("--chrom",   required=True, help="Cromosoma: 1-22 o X")
    p.add_argument(
        "--missing-strategy",
        choices=["mean", "zero", "skip"],
        default="mean",
        help=(
            "Estrategia para dosage faltante (NaN). "
            "mean=imputar con media, zero=contar como 0, skip=excluir variante. "
            "(default: mean)"
        ),
    )
    p.add_argument(
        "--window-size",
        type=int,
        default=10_000,
        help=(
            "Variantes zarr por ventana de lectura. "
            "Controla peak de RAM CPU: W × n_samples × dtype_bytes. "
            "La VRAM solo carga filas matched, no la ventana entera. "
            "(default: 10000)"
        ),
    )
    p.add_argument(
        "--spark-master",
        default="local[*]",
        help="Spark master URL. 'local[N]' usa N cores locales. (default: local[*])",
    )
    p.add_argument(
        "--spark-driver-memory",
        default="8g",
        help="Memoria del driver Spark (default: 8g)",
    )
    p.add_argument(
        "--spark-executor-memory",
        default="16g",
        help="Memoria de executors Spark (default: 16g)",
    )
    p.add_argument(
        "--gpu-batch-matches",
        type=int,
        default=4_000,
        help=(
            "Variantes matched que se acumulan en CPU antes de cada cp.dot. "
            "Controla el tamaño del bloque GPU: N × n_samples × 4 B de VRAM. "
            "RTX 3060 12 GB / 140k muestras: 4000→2.25 GB, 8000→4.5 GB. "
            "Si supera el límite seguro se reduce automáticamente. "
            "(default: 4000)"
        ),
    )
    return p.parse_args()


# =============================================================================
# Utilidades de alelos
# =============================================================================

def _decode_allele(val) -> str:
    if isinstance(val, (bytes, np.bytes_)):
        return val.decode("ascii", errors="replace").upper()
    return str(val).upper()


def _comp(allele: str) -> str:
    """Complemento de cadena de ADN (A↔T, C↔G). Maneja SNPs y multi-char."""
    return "".join(_COMPLEMENT.get(b, b) for b in allele.upper())


# =============================================================================
# FASE 1a — Cargar pesos del PRS (betamap)
# =============================================================================

def load_score_weights(pgs_id: str, chrom: str) -> pd.DataFrame:
    """
    Lee el betamap pre-procesado del PGS Catalog.
    Formato: CHROM, POS, IS_FLIP, BETA, EFFECT_ALLELE, OTHER_ALLELE, ID
    IS_FLIP=0 → efecto=ALT; IS_FLIP=1 → efecto=REF
    """
    path = SCORES_BASE / pgs_id / f"{pgs_id}_hmPOS_GRCh38.betamap.tsv.gz"
    if not path.exists():
        log.error("Archivo de pesos no encontrado: %s", path)
        sys.exit(1)
    log.info("Leyendo betamap desde %s", path)
    df = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        dtype={
            "CHROM":         str,
            "POS":           np.int64,
            "IS_FLIP":       np.int8,
            "BETA":          np.float64,
            "EFFECT_ALLELE": str,
            "OTHER_ALLELE":  str,
            "ID":            str,
        },
    )
    subset = df[df["CHROM"] == str(chrom)].copy().reset_index(drop=True)
    log.info("  Variantes en chr%s: %d", chrom, len(subset))
    return subset


# =============================================================================
# FASE 1b — Extraer catálogo de variantes del zarr (CPU, pandas)
# =============================================================================

def load_zarr_catalog(chrom: str) -> tuple[zarr.Group, pd.DataFrame]:
    """
    Abre el zarr del cromosoma y extrae el catálogo de variantes:
      zarr_idx (int64), POS (int64), REF (str), ALT (str)

    No carga los dosages: solo pos/alleles (metadata ligero).
    Retorna (store, catalog_df) para reutilizar el store en la Fase 2.
    """
    zarr_path = ZARR_BASE / f"chr{chrom}.zarr"
    if not zarr_path.exists():
        log.error("Zarr store no encontrado: %s", zarr_path)
        sys.exit(1)
    log.info("Abriendo zarr: %s", zarr_path)
    store = zarr.open(str(zarr_path), mode="r")

    pos     = store["variant_position"][:]
    alleles = store["variant_allele"][:]
    ref     = np.array([_decode_allele(v) for v in alleles[:, 0]])
    alt     = np.array([_decode_allele(v) for v in alleles[:, 1]])

    catalog = pd.DataFrame({
        "zarr_idx": np.arange(len(pos), dtype=np.int64),
        "POS":      pos.astype(np.int64),
        "REF":      ref,
        "ALT":      alt,
    })
    log.info("  Catálogo zarr chr%s: %d variantes", chrom, len(catalog))
    return store, catalog


# =============================================================================
# FASE 1c — Matching Spark: broadcast join con resolución de alelos
# =============================================================================

def spark_match(
    spark,
    df_weights_pd: pd.DataFrame,
    df_zarr_pd: pd.DataFrame,
) -> tuple[pd.DataFrame, list[dict]]:
    """
    Broadcast join Spark para alinear variantes PRS con el catálogo zarr.

    Lógica de matching (replica compute_prs.py / match_via_duckdb):
    ──────────────────────────────────────────────────────────────────
    IS_FLIP=0 (efecto=ALT): exp_ref=OA, exp_alt=EA
      Arm A_direct : zarr.REF==OA       AND zarr.ALT==EA       → flip_gpu=False
      Arm A_strand : zarr.REF==comp(OA) AND zarr.ALT==comp(EA) → flip_gpu=False

    IS_FLIP=1 (efecto=REF): exp_ref=EA, exp_alt=OA
      Arm B_direct : zarr.REF==EA       AND zarr.ALT==OA       → flip_gpu=True
      Arm B_strand : zarr.REF==comp(EA) AND zarr.ALT==comp(OA) → flip_gpu=True

    Deduplicación: para la misma variante PRS (ID), se prefiere:
      1. match directo (is_strand=0) sobre strand-flip (is_strand=1)
      2. zarr_idx menor si aún hay empate (primer candidato en el zarr)

    Retorna (matched_pd, excluded_list).
    matched_pd columnas: zarr_idx, ID, BETA, flip_gpu
    excluded_list: mismos campos que compute_prs.py excluded_allele_mismatch
    """
    from pyspark.sql import functions as F
    from pyspark.sql.functions import broadcast
    from pyspark.sql.window import Window

    # ── Añadir columnas derivadas al betamap ─────────────────────────────────
    wdf = df_weights_pd.copy()
    wdf["EA"]      = wdf["EFFECT_ALLELE"].str.upper()
    wdf["OA"]      = wdf["OTHER_ALLELE"].str.upper()
    wdf["EA_comp"] = wdf["EA"].map(_comp)
    wdf["OA_comp"] = wdf["OA"].map(_comp)

    # IS_FLIP asegurado como int (viene del betamap)
    wdf["IS_FLIP"]  = wdf["IS_FLIP"].astype(int)

    # Row ID sintético garantizado único: necesario cuando la columna ID tiene
    # nulos (ej. PGS000363 — toda la columna ID viene vacía del PGS Catalog).
    # Sin esto, Window.partitionBy("ID") agrupa TODOS los nulos en una sola
    # partición y row_number()==1 devuelve solo 1 variante de las 89k+.
    wdf["_betamap_rowid"] = np.arange(len(wdf), dtype=np.int64)

    # Spark DataFrames
    # df_zarr es el DataFrame grande (millones de variantes) → distribuido
    # df_w es el betamap pequeño (miles de variantes) → broadcast
    df_z = spark.createDataFrame(df_zarr_pd)
    df_w = broadcast(spark.createDataFrame(wdf))

    Z = df_z.alias("Z")

    join_pos = F.col("Z.POS") == F.col("W.POS")

    def arm(ref_col: str, alt_col: str, is_flip_val: int, is_strand_val: int, flip_gpu: bool):
        """Un brazo del join: filtra IS_FLIP y une por posición + alelos."""
        W = df_w.filter(F.col("IS_FLIP") == is_flip_val).alias("W")
        return (
            Z.join(
                W,
                join_pos
                & (F.col("Z.REF") == F.col(f"W.{ref_col}"))
                & (F.col("Z.ALT") == F.col(f"W.{alt_col}")),
                "inner",
            ).select(
                F.col("Z.zarr_idx").alias("zarr_idx"),
                F.col("W._betamap_rowid").alias("_betamap_rowid"),
                F.col("W.ID").alias("ID"),
                F.col("W.POS").alias("POS"),
                F.col("W.BETA").alias("BETA"),
                F.lit(int(flip_gpu)).cast("integer").alias("flip_gpu"),
                F.lit(is_strand_val).cast("integer").alias("is_strand"),
            )
        )

    # IS_FLIP=0: exp_ref=OA, exp_alt=EA
    a_direct = arm("OA",      "EA",      is_flip_val=0, is_strand_val=0, flip_gpu=False)
    a_strand = arm("OA_comp", "EA_comp", is_flip_val=0, is_strand_val=1, flip_gpu=False)

    # IS_FLIP=1: exp_ref=EA, exp_alt=OA
    b_direct = arm("EA",      "OA",      is_flip_val=1, is_strand_val=0, flip_gpu=True)
    b_strand = arm("EA_comp", "OA_comp", is_flip_val=1, is_strand_val=1, flip_gpu=True)

    all_matches = (
        a_direct
        .unionByName(a_strand)
        .unionByName(b_direct)
        .unionByName(b_strand)
    )

    # Deduplicar: por variante PRS (_betamap_rowid, siempre único),
    # preferir match directo (is_strand=0) sobre strand-flip (is_strand=1),
    # luego zarr_idx menor si aún hay empate.
    # NOTA: NO usar ID aquí — para muchos PGS (ej. PGS000363) la columna ID
    # viene vacía/null del Catalog, lo que colapsaría TODAS las variantes en
    # una sola partición y devolvería solo 1 ganadora de las 89k+.
    w_dedup = Window.partitionBy("_betamap_rowid").orderBy(
        F.col("is_strand").asc(),
        F.col("zarr_idx").asc(),
    )
    matched_spark = (
        all_matches
        .withColumn("_rn", F.row_number().over(w_dedup))
        .filter(F.col("_rn") == 1)
        .drop("_rn", "is_strand", "POS")
        # Conservar _betamap_rowid en el resultado para rastrear exactamente
        # qué filas del betamap emparejaron (necesario cuando ID es null).
    )

    matched_pd = matched_spark.toPandas()
    matched_pd["flip_gpu"] = matched_pd["flip_gpu"].astype(bool)
    log.info("  Spark: %d variantes emparejadas", len(matched_pd))

    # ── Diagnóstico de excluidos: desglose por razón de exclusión ───────────
    # Join solo por POS: identifica qué posiciones del betamap SÍ existen en
    # zarr pero fallan el filtro de alelos (vs. las que no están en zarr).
    pos_join_diag = (
        Z.join(
            df_w.select("_betamap_rowid", "POS").alias("Wdiag"),
            F.col("Z.POS") == F.col("Wdiag.POS"),
            "right",
        )
        .select(
            F.col("Wdiag._betamap_rowid").alias("_betamap_rowid"),
            (F.col("Z.zarr_idx").isNotNull()).alias("pos_in_zarr"),
        )
        .distinct()
        .toPandas()
    )
    rowids_pos_in_zarr: set = set(
        pos_join_diag.loc[pos_join_diag["pos_in_zarr"], "_betamap_rowid"].tolist()
    )

    matched_rowids: set = set(matched_pd["_betamap_rowid"].tolist())

    excluded: list[dict] = []
    n_excl_no_pos = 0
    n_excl_allele = 0
    for _, row in wdf.iterrows():
        rowid = int(row["_betamap_rowid"])
        if rowid in matched_rowids:
            continue
        pos = int(row["POS"])
        if rowid in rowids_pos_in_zarr:
            n_excl_allele += 1
            reason = "discordancia_de_alelos"
        else:
            n_excl_no_pos += 1
            reason = "posicion_no_encontrada_en_zarr"
        excluded.append({
            "ID":    row["ID"],
            "POS":   pos,
            "reason": reason,
        })

    log.info(
        "  Spark excluidos: %d total  |  %d posición-no-en-zarr  |  %d discordancia-alelos",
        len(excluded), n_excl_no_pos, n_excl_allele,
    )

    # Eliminar columna auxiliar del resultado final (no forma parte del contrato)
    matched_pd = matched_pd.drop(columns=["_betamap_rowid"])

    return matched_pd, excluded


# =============================================================================
# FASE 2 — Scoring GPU con CuPy por ventanas zarr
# =============================================================================

def _detect_dosage_key(store: zarr.Group) -> tuple[str, str]:
    """
    Detecta la clave del array de dosage y su layout.
    layout 'vs' → (n_variants, n_samples)
    layout 'sv' → (n_samples, n_variants)
    layout 'gt' → genotype (n_variants, n_samples, ploidy) → se convierte a dosage
    """
    for key in ("call_DS", "call_dosage", "calldata/DS", "dosage"):
        if key in store:
            shape = store[key].shape
            return key, "vs" if shape[0] >= shape[1] else "sv"
    for key in ("call_genotype", "calldata/GT"):
        if key in store:
            return key, "gt"
    raise KeyError(
        "Array de dosage/genotipo no encontrado. "
        "Candidatos: call_DS, call_dosage, calldata/DS, call_genotype, calldata/GT."
    )


def get_samples(store: zarr.Group) -> np.ndarray:
    """Extrae el array de IDs de muestra del zarr."""
    for key in ("samples", "sample_id", "sample/id", "calldata/samples"):
        if key in store:
            raw = store[key][:]
            return np.array([
                v.decode("ascii") if isinstance(v, bytes) else str(v)
                for v in raw
            ])
    raise KeyError("Array de muestras no encontrado en el zarr.")


def _read_window_cpu(
    store: zarr.Group,
    dosage_key: str,
    layout: str,
    w_start: int,
    w_end: int,
) -> np.ndarray:
    """
    Lee una ventana [w_start, w_end) del eje variantes.
    Retorna (w_size, n_samples) en float32.
    - layout 'vs': squeeze de la dim de ploidía como view (cero copias extra).
      Solo hace cast si el dtype no es ya float32.
    - layout 'gt': convierte genotype int8 → dosage float32 en CPU.
    """
    arr = store[dosage_key]

    if layout == "vs":
        raw = arr[w_start:w_end]
        if raw.ndim == 3:
            raw = raw[:, :, 0]   # squeeze ploidy trailing dim (view, sin copia)
        return raw if raw.dtype == np.float32 else raw.astype(np.float32)

    if layout == "sv":
        return arr[:, w_start:w_end].T.astype(np.float32)

    if layout == "gt":
        raw         = arr[w_start:w_end]              # (W, n_samples, ploidy) int8/int16
        missing_gt  = (raw < 0).any(axis=2)           # (W, n_samples) bool
        dosage      = raw.clip(min=0).sum(axis=2).astype(np.float32)
        dosage[missing_gt] = np.nan
        return dosage

    raise ValueError(f"Layout desconocido: {layout!r}")


def compute_prs_gpu(
    store: zarr.Group,
    dosage_key: str,
    layout: str,
    n_zarr_variants: int,
    zarr_indices: np.ndarray,       # int, ordenados ascendente
    betas: np.ndarray,              # float64
    flip_mask: np.ndarray,          # bool
    n_samples: int,
    window_size: int,
    missing_strategy: str,
    excluded_missing: list,
    gpu_batch_matches: int = 4_000, # variantes matched por bloque GPU
) -> np.ndarray:
    """
    Acumula PRS en GPU iterando sobre ventanas zarr con batching por bloque.

    Estrategia de throughput GPU:
    ──────────────────────────────
    - Loop externo → ventanas zarr [w_start, w_end)
    - Ventanas sin matches → skip (cero I/O zarr ni GPU)
    - Las filas matched se acumulan en buffers CPU (rows_buf, betas_buf)
      hasta alcanzar gpu_batch_matches filas o llegar a la última ventana.
    - Un solo cp.dot por flush: (gpu_batch_matches, n_samples) en lugar de
      ~N kernels de (n_match_per_window, n_samples), reduciendo lanzamientos
      de kernel y mejorando la ocupación de warps en la GPU.
    - Peak VRAM ≈ gpu_batch_matches × n_samples × 4 B (float32 G_gpu)
                + n_samples × 8 B (float64 prs_gpu, siempre residente)
    """
    import cupy as cp

    # ── Validación de VRAM en runtime ────────────────────────────────────────
    free_vram_bytes, total_vram_bytes = cp.cuda.runtime.memGetInfo()
    # Reservar: 512 MB para driver/pool de CuPy + prs_gpu (float64)
    _reserve   = 512 * 1024 * 1024 + n_samples * 8
    _safe_batch = max(1, int((free_vram_bytes - _reserve) / (n_samples * 4)))
    log.info(
        "  [VRAM] libre=%.2f GB  total=%.2f GB  |  "
        "límite seguro gpu_batch_matches=%d",
        free_vram_bytes / 1e9, total_vram_bytes / 1e9, _safe_batch,
    )
    if gpu_batch_matches > _safe_batch:
        log.warning(
            "  [VRAM] gpu_batch_matches=%d supera límite seguro=%d  "
            "(%.2f GB → %.2f GB por bloque) — reducido automáticamente.",
            gpu_batch_matches, _safe_batch,
            gpu_batch_matches * n_samples * 4 / 1e9,
            _safe_batch       * n_samples * 4 / 1e9,
        )
        gpu_batch_matches = _safe_batch

    prs_gpu = cp.zeros(n_samples, dtype=cp.float64)
    n_matched       = len(zarr_indices)
    n_windows_total = (n_zarr_variants + window_size - 1) // window_size
    n_windows_read  = 0
    n_windows_skip  = 0
    ptr_lo          = 0

    # Buffers CPU: acumulan filas matched de múltiples ventanas zarr
    rows_buf:  list[np.ndarray] = []   # (n_match_i, n_samples) float32 por ventana
    betas_buf: list[np.ndarray] = []   # (n_match_i,)            float32 por ventana
    n_buf         = 0
    n_gpu_flushes = 0

    def _flush_to_gpu() -> None:
        """Transfiere el buffer acumulado a GPU, ejecuta cp.dot y libera."""
        nonlocal n_buf, n_gpu_flushes, prs_gpu
        if n_buf == 0:
            return
        # Concatenar y garantizar shapes exactos antes de subir a GPU
        G_big = np.concatenate(rows_buf,  axis=0).reshape(n_buf, n_samples)
        W_big = np.concatenate(betas_buf, axis=0).reshape(n_buf)
        mem_gb = n_buf * n_samples * 4 / 1e9
        log.info(
            "  [GPU flush #%d] n_buf=%d  G_big.shape=%s  "
            "W_big.shape=%s  mem_bloque=%.3f GB",
            n_gpu_flushes + 1, n_buf, G_big.shape, W_big.shape, mem_gb,
        )
        G_gpu = cp.asarray(G_big, dtype=cp.float32);  del G_big
        W_gpu = cp.asarray(W_big, dtype=cp.float32);  del W_big
        prs_gpu += cp.dot(W_gpu, G_gpu)   # (n_samples,) acumulado en float64
        del G_gpu, W_gpu
        cp.get_default_memory_pool().free_all_blocks()
        rows_buf.clear()
        betas_buf.clear()
        n_buf         = 0
        n_gpu_flushes += 1

    for w_i in range(n_windows_total):
        w_start = w_i * window_size
        w_end   = min(w_start + window_size, n_zarr_variants)

        # Avanzar puntero al primer match no procesado en esta ventana
        while ptr_lo < n_matched and zarr_indices[ptr_lo] < w_start:
            ptr_lo += 1

        ptr_hi = ptr_lo
        while ptr_hi < n_matched and zarr_indices[ptr_hi] < w_end:
            ptr_hi += 1

        if ptr_lo == ptr_hi:
            n_windows_skip += 1
            continue   # ← cero I/O zarr ni GPU; flush se hace tras el loop

        # ── Leer ventana zarr → CPU RAM ──────────────────────────────────────
        window_data = _read_window_cpu(store, dosage_key, layout, w_start, w_end)
        n_windows_read += 1

        # ── Extraer solo filas matched (micro-slice en CPU) ──────────────────
        n_match       = ptr_hi - ptr_lo
        local_indices = zarr_indices[ptr_lo:ptr_hi] - w_start   # índices relativos
        chunk_rows    = window_data[local_indices]   # fancy-index produce array propio, .copy() redundante
        del window_data

        # ── Flip dosage (efecto=REF): 2 − DS ─────────────────────────────────
        flip_chunk = flip_mask[ptr_lo:ptr_hi]
        if flip_chunk.any():
            chunk_rows[flip_chunk] = 2.0 - chunk_rows[flip_chunk]

        # ── NaN / missing dosage (vectorizado en C, sin bucle Python) ──────────
        betas_chunk = betas[ptr_lo:ptr_hi].copy()
        nan_mask_2d = np.isnan(chunk_rows)

        if nan_mask_2d.any():
            if missing_strategy == "mean":
                # np.nanmean por fila en una sola pasada C; keepdims permite
                # broadcast directo contra la matriz (n_match, n_samples).
                row_means  = np.nanmean(chunk_rows, axis=1, keepdims=True)  # (n_match, 1)
                chunk_rows = np.where(nan_mask_2d, row_means, chunk_rows)   # (n_match, n_samples)

            elif missing_strategy == "zero":
                chunk_rows[nan_mask_2d] = 0.0

            elif missing_strategy == "skip":
                # 'skip' requiere registro individual por variante y modifica
                # betas_chunk → se mantiene el bucle solo para este caso.
                for m_local in range(n_match):
                    nm = nan_mask_2d[m_local]
                    if not nm.any():
                        continue
                    excluded_missing.append({
                        "zarr_index":        int(zarr_indices[ptr_lo + m_local]),
                        "n_missing_samples": int(nm.sum()),
                        "reason":            "skip_por_missing",
                    })
                    chunk_rows[m_local]  = 0.0
                    betas_chunk[m_local] = 0.0   # beta=0 → contribución nula

        # ── Normalizar shape y acumular en buffer CPU ─────────────────────────
        # reshape garantiza (n_match, n_samples) aunque chunk_rows llegue con
        # eje extra (ej. (n_match, n_samples, 1) en algunos layouts de zarr).
        rows_buf.append(chunk_rows.reshape(n_match, n_samples).astype(np.float32))
        betas_buf.append(betas_chunk.reshape(n_match).astype(np.float32))
        n_buf += n_match

        # ── Flush si el buffer alcanzó el umbral ─────────────────────────────
        if n_buf >= gpu_batch_matches:
            _flush_to_gpu()

        if n_windows_read % 50 == 0:
            pct = 100.0 * w_end / n_zarr_variants
            log.info(
                "  ventana %d/%d (%.1f%%) | leídas=%d  saltadas=%d  "
                "buf=%d  flushes=%d",
                w_i + 1, n_windows_total, pct,
                n_windows_read, n_windows_skip, n_buf, n_gpu_flushes,
            )

    # ── Flush del buffer residual (última(s) ventana(s) < gpu_batch_matches) ──
    _flush_to_gpu()

    log.info(
        "  Fin GPU: leídas=%d  saltadas=%d  flushes=%d  "
        "(%.1f%% de ventanas sin I/O)",
        n_windows_read, n_windows_skip, n_gpu_flushes,
        100.0 * n_windows_skip / max(n_windows_total, 1),
    )

    prs_cpu = cp.asnumpy(prs_gpu)
    del prs_gpu
    cp.get_default_memory_pool().free_all_blocks()
    return prs_cpu


# =============================================================================
# FASE 3 — Persistencia (contrato idéntico a compute_prs.py / save_results)
# =============================================================================

def save_results(
    pgs_id: str,
    chrom: str,
    samples: np.ndarray,
    prs: np.ndarray,
    metadata: dict,
) -> None:
    """
    Escribe los mismos dos archivos que compute_prs.py:
      {pgs_id}_chr{chrom}_scores.tsv   — sample_id, PRS (%.8f)
      {pgs_id}_chr{chrom}_metadata.json
    Destino: OUT_BASE/{pgs_id}/  (home del usuario, no hot_nvme)
    """
    out_dir = OUT_BASE / pgs_id
    out_dir.mkdir(parents=True, exist_ok=True)

    scores_path = out_dir / f"{pgs_id}_chr{chrom}_scores.tsv"
    pd.DataFrame({"sample_id": samples, "PRS": prs}).to_csv(
        scores_path, sep="\t", index=False, float_format="%.8f"
    )
    log.info("Scores guardados: %s", scores_path)

    meta_path = out_dir / f"{pgs_id}_chr{chrom}_metadata.json"
    with open(meta_path, "w") as fh:
        json.dump(metadata, fh, indent=2, default=str)
    log.info("Metadatos guardados: %s", meta_path)


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    args    = parse_args()
    pgs_id  = args.pgs_id
    chrom   = args.chrom
    t_start = datetime.now()

    # ── 1. Pesos del PRS (betamap) ────────────────────────────────────────────
    weights = load_score_weights(pgs_id, chrom)
    if weights.empty:
        log.warning("Sin variantes para chr%s en %s — saltando.", chrom, pgs_id)
        sys.exit(0)

    # ── 2. Catálogo zarr para este cromosoma ──────────────────────────────────
    log.info("[Fase 1] Extrayendo catálogo zarr para chr%s...", chrom)
    store, zarr_catalog = load_zarr_catalog(chrom)

    # ── 3. Sesión Spark ───────────────────────────────────────────────────────
    log.info("[Fase 1] Iniciando SparkSession (master=%s)...", args.spark_master)
    try:
        from pyspark.sql import SparkSession
    except ImportError:
        log.error(
            "PySpark no disponible. Instalar con: pip install pyspark\n"
            "O verificar que el entorno del contenedor incluya pyspark."
        )
        sys.exit(1)

    spark = (
        SparkSession.builder
        .appName(f"PRS_Spark_GPU_{pgs_id}_chr{chrom}")
        .master(args.spark_master)
        .config("spark.driver.memory",                    args.spark_driver_memory)
        .config("spark.executor.memory",                  args.spark_executor_memory)
        # Umbral de broadcast ampliado: el betamap es siempre pequeño (<< 512 MB)
        .config("spark.sql.autoBroadcastJoinThreshold",   str(512 * 1024 * 1024))
        # Shuffle partitions = CPUs Slurm disponibles
        .config("spark.sql.shuffle.partitions",           "16")
        # Puertos dinámicos: evita colisiones entre jobs paralelos en el mismo nodo
        .config("spark.ui.port",                          "0")
        .config("spark.driver.port",                      "0")
        .config("spark.port.maxRetries",                  "32")
        # Necesario en contenedores: el driver debe escuchar en loopback,
        # no intentar resolver el hostname del contenedor contra el DNS del cluster
        .config("spark.driver.host",                      "127.0.0.1")
        .config("spark.driver.bindAddress",               "127.0.0.1")
        .getOrCreate()
    )
    spark.sparkContext.setLogLevel("WARN")

    # ── 4. Matching Spark ─────────────────────────────────────────────────────
    log.info("[Fase 1] Ejecutando broadcast join (4 arms: IS_FLIP × strand)...")
    matched_pd, excluded_allele = spark_match(spark, weights, zarr_catalog)
    spark.stop()
    log.info("[Fase 1] Spark detenido.")

    if matched_pd.empty:
        log.error("Cero variantes emparejadas para chr%s en %s — abortando.", chrom, pgs_id)
        sys.exit(1)

    # Ordenar por zarr_idx: requerido por el loop de ventanas (acceso secuencial)
    matched_pd   = matched_pd.sort_values("zarr_idx").reset_index(drop=True)
    zarr_indices = matched_pd["zarr_idx"].values.astype(np.intp)
    betas        = matched_pd["BETA"].values.astype(np.float64)
    flip_mask    = matched_pd["flip_gpu"].values.astype(bool)

    log.info(
        "[Fase 1] Matching completado: %d / %d variantes emparejadas  |  %d excluidas",
        len(matched_pd), len(weights), len(excluded_allele),
    )

    # ── 5. Zarr: detección de layout y dimensiones ────────────────────────────
    samples          = get_samples(store)
    n_samples        = len(samples)
    dosage_key, layout = _detect_dosage_key(store)
    ax               = 0 if layout in ("vs", "gt") else 1
    n_zarr_variants  = store[dosage_key].shape[ax]
    dtype_bytes      = store[dosage_key].dtype.itemsize

    log.info(
        "Zarr: %d variantes × %d muestras  |  key='%s'  layout=%s",
        n_zarr_variants, n_samples, dosage_key, layout,
    )
    log.info(
        "Window W=%d  gpu_batch_matches=%d  |  "
        "Peak RAM CPU estimado: %.2f GB  |  "
        "Peak VRAM GPU estimado (bloque): %.2f GB",
        args.window_size,
        args.gpu_batch_matches,
        args.window_size * n_samples * dtype_bytes / 1e9,
        args.gpu_batch_matches * n_samples * 4 / 1e9,
    )

    # ── 6. Scoring GPU (CuPy) ─────────────────────────────────────────────────
    log.info("[Fase 2] Calculando PRS en GPU (CuPy)...")
    excluded_missing: list[dict] = []
    prs = compute_prs_gpu(
        store=store,
        dosage_key=dosage_key,
        layout=layout,
        n_zarr_variants=n_zarr_variants,
        zarr_indices=zarr_indices,
        betas=betas,
        flip_mask=flip_mask,
        n_samples=n_samples,
        window_size=args.window_size,
        missing_strategy=args.missing_strategy,
        excluded_missing=excluded_missing,
        gpu_batch_matches=args.gpu_batch_matches,
    )

    elapsed = (datetime.now() - t_start).total_seconds()

    # ── 7. Metadatos — MISMO ESQUEMA QUE compute_prs.py (contrato congelado) ──
    metadata = {
        # Identificadores
        "pgs_id":                       pgs_id,
        "chrom":                        chrom,
        # Timing
        "run_timestamp":                t_start.isoformat(),
        "elapsed_seconds":              round(elapsed, 2),
        # Parámetros de ejecución
        "missing_strategy":             args.missing_strategy,
        "window_size":                  args.window_size,
        "dosage_array_key":             dosage_key,
        "dosage_layout":                layout,
        "index_path":                   "spark_variant_catalog",  # antes: DuckDB path
        # Dimensiones
        "n_samples":                    int(n_samples),
        "n_zarr_variants":              int(n_zarr_variants),
        "n_variants_in_score_file":     int(len(weights)),
        "n_variants_matched":           int(len(matched_pd)),
        "n_variants_excluded_allele":   int(len(excluded_allele)),
        "n_variants_excluded_missing":  int(len(excluded_missing)),
        # Estadísticas del PRS (mismas que compute_prs.py)
        "prs_mean":                     float(np.mean(prs)),
        "prs_std":                      float(np.std(prs)),
        "prs_min":                      float(np.min(prs)),
        "prs_max":                      float(np.max(prs)),
        # Listas de excluidos (misma estructura que compute_prs.py)
        "excluded_allele_mismatch":     excluded_allele,
        "excluded_missing":             excluded_missing,
    }

    # ── 8. Guardar (Fase 3) ───────────────────────────────────────────────────
    log.info("[Fase 3] Guardando resultados...")
    save_results(pgs_id, chrom, samples, prs, metadata)
    log.info("chr%s completado en %.1f s ✓", chrom, elapsed)


if __name__ == "__main__":
    main()
