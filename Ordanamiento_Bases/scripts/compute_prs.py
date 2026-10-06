#!/usr/bin/env python3
"""
compute_prs.py — Cómputo de Polygenic Risk Score para un cromosoma.

Modelo aditivo:
    PRS = Σ_i  BETA_i × dosage_i

donde dosage_i es el número de copias del alelo de efecto (0, 1 ó 2):
  - IS_FLIP=0 → efecto=ALT  → dosage = dosage ALT (lectura directa del zarr)
  - IS_FLIP=1 → efecto=REF  → dosage = 2 − dosage_ALT

Estrategia de acceso al zarr:
    El loop externo itera en ventanas de W variantes zarr consecutivas.
    Solo se emite I/O para ventanas que contienen ≥1 variante del PRS.
    El matching posición-alelo se resuelve con una consulta al índice DuckDB
    (construido previamente por build_zarr_index.py); no se carga nunca el
    array completo pos/ref/alt del zarr en memoria.

Uso:
    python compute_prs.py --pgs-id PGS000001 --chrom 1
    python compute_prs.py --pgs-id PGS000001 --chrom X --window-size 20000
"""

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import zarr
import pdb

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stdout,
)
log = logging.getLogger(__name__)

SCORES_BASE = Path("/mnt/cephfs/hot_nvme/pgscatalog/scores")
ZARR_BASE = Path("/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zar_files")
INDEX_PATH = Path("/mnt/cephfs/hot_nvme/mcps/imputed-topmed/zarr_index.duckdb")


# =============================================================================
# Argumentos
# =============================================================================


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Calcula PRS por cromosoma usando modelo aditivo."
    )
    p.add_argument("--pgs-id", required=True, help="PGS Catalog ID, ej. PGS000001")
    p.add_argument("--chrom", required=True, help="Cromosoma: 1-22 o X")
    p.add_argument(
        "--missing-strategy",
        choices=["mean", "zero", "skip"],
        default="mean",
        help="Estrategia para dosage faltante (NaN). "
        "mean=imputar con media, zero=contar como 0, skip=excluir variante. "
        "(default: mean)",
    )
    p.add_argument(
        "--window-size",
        type=int,
        default=10_000,
        help="Variantes zarr por ventana de lectura. "
        "Controla el peak de RAM: peak ≈ W × n_muestras × bytes_dtype. "
        "Calcular con: W = RAM_disponible / (n_muestras × bytes_dtype). "
        "(default: 10000)",
    )
    p.add_argument(
        "--index-path",
        type=Path,
        default=INDEX_PATH,
        help="Ruta del índice DuckDB generado por build_zarr_index.py",
    )
    return p.parse_args()


# =============================================================================
# Pesos del PRS
# =============================================================================


def load_score_weights(pgs_id: str, chrom: str) -> pd.DataFrame:
    path = SCORES_BASE / pgs_id / f"{pgs_id}_hmPOS_GRCh38.betamap.tsv.gz"
    if not path.exists():
        log.error("Archivo de pesos no encontrado: %s", path)
        sys.exit(1)
    log.info("Leyendo pesos desde %s", path)
    df = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        dtype={"CHROM": str, "POS": np.int64, "IS_FLIP": np.int8, "BETA": np.float64},
    )
    subset = df[df["CHROM"] == str(chrom)].copy().reset_index(drop=True)
    log.info("  Variantes en chr%s: %d", chrom, len(subset))
    return subset


# =============================================================================
# Matching vía DuckDB
# =============================================================================


def _decode_allele(val) -> str:
    if isinstance(val, (bytes, np.bytes_)):
        return val.decode("ascii", errors="replace").upper()
    return str(val).upper()


def match_via_duckdb(
    weights: pd.DataFrame,
    chrom: str,
    index_path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[dict]]:
    """
    Resuelve el matching PRS-weights - zarr usando el índice DuckDB.

    1. Consulta DuckDB por las POS del PRS → devuelve solo las filas relevantes
       (n_PRS filas, no n_zarr filas): zarr_idx, ref, alt por posición.
    2. Verifica concordancia de alelos + strand-flip en Python sobre ese subset.
    3. Sitios con múltiples zarr_idx (duplicados en el zarr) se resuelven
       eligiendo el que tiene concordancia exacta de alelos.

    Retorna arrays (zarr_indices, betas, flip_mask, excluded) ordenados por
    zarr_idx ascendente, listos para el loop por ventanas zarr.
    """
    if not index_path.exists():
        log.error(
            "Índice DuckDB no encontrado: %s\n" "Ejecuta build_zarr_index.py primero.",
            index_path,
        )
        sys.exit(1)

    prs_positions = weights["POS"].tolist()

    con = duckdb.connect(str(index_path), read_only=True)
    rows = con.execute(
        """
        SELECT zarr_idx, pos, ref, alt
        FROM   zarr_variants
        WHERE  chrom = ?
          AND  pos   IN (SELECT UNNEST(?::BIGINT[]))
        ORDER  BY zarr_idx
    """,
        [chrom, prs_positions],
    ).fetchall()
    con.close()

    log.info("  Candidatos en DuckDB para chr%s: %d", chrom, len(rows))

    if not rows:
        log.warning("  Ninguna posición encontrada en el índice para chr%s", chrom)
        return np.array([], dtype=np.intp), np.array([]), np.array([], dtype=bool), []

    # Indexar candidatos zarr por posición para join rápido
    zarr_by_pos: dict[int, list[tuple]] = {}
    for zarr_idx, pos, ref, alt in rows:
        zarr_by_pos.setdefault(int(pos), []).append((int(zarr_idx), ref, alt))

    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}

    matched_zarr: list[int] = []
    matched_betas: list[float] = []
    matched_flip: list[bool] = []
    excluded: list[dict] = []

    for _, row in weights.iterrows():
        pos = int(row["POS"])
        ea = _decode_allele(row["EFFECT_ALLELE"])
        oa = _decode_allele(row["OTHER_ALLELE"])
        is_flip = int(row["IS_FLIP"])
        beta = float(row["BETA"])

        candidates = zarr_by_pos.get(pos)
        if candidates is None:
            excluded.append(
                {
                    "ID": row["ID"],
                    "POS": pos,
                    "reason": "posicion_no_encontrada_en_zarr",
                }
            )
            continue

        # IS_FLIP=0 → efecto=ALT → zarr_REF==OA, zarr_ALT==EA
        # IS_FLIP=1 → efecto=REF → zarr_REF==EA, zarr_ALT==OA
        exp_ref = oa if is_flip == 0 else ea
        exp_alt = ea if is_flip == 0 else oa

        best = None
        for z_idx, z_ref, z_alt in candidates:
            z_ref = _decode_allele(z_ref)
            z_alt = _decode_allele(z_alt)
            if z_ref == exp_ref and z_alt == exp_alt:
                best = z_idx
                break
            flip_ref = "".join(complement.get(b, b) for b in exp_ref)
            flip_alt = "".join(complement.get(b, b) for b in exp_alt)
            if z_ref == flip_ref and z_alt == flip_alt:
                log.debug("Strand flip en POS=%d — aceptado", pos)
                best = z_idx
                break

        if best is None:
            excluded.append(
                {
                    "ID": row["ID"],
                    "POS": pos,
                    "reason": "discordancia_de_alelos",
                    "zarr_candidates": [(z, r, a) for z, r, a in candidates],
                    "prs_EFFECT": ea,
                    "prs_OTHER": oa,
                    "IS_FLIP": is_flip,
                }
            )
            continue

        matched_zarr.append(best)
        matched_betas.append(beta)
        matched_flip.append(is_flip == 1)

    log.info(
        "  Emparejadas: %d / %d  |  Excluidas: %d",
        len(matched_zarr),
        len(weights),
        len(excluded),
    )

    order = np.argsort(matched_zarr)
    return (
        np.array(matched_zarr, dtype=np.intp)[order],
        np.array(matched_betas, dtype=np.float64)[order],
        np.array(matched_flip, dtype=bool)[order],
        excluded,
    )


# =============================================================================
# Acceso al zarr (solo dosage — no se cargan pos/ref/alt)
# =============================================================================


def open_zarr_store(chrom: str) -> zarr.Group:
    path = ZARR_BASE / f"chr{chrom}.zarr"
    if not path.exists():
        log.error("Zarr store no encontrado: %s", path)
        sys.exit(1)
    log.info("Abriendo zarr store: %s", path)
    return zarr.open(str(path), mode="r")


def get_samples(store: zarr.Group) -> np.ndarray:
    candidates = ["samples", "sample_id", "sample/id", "calldata/samples"]
    for key in candidates:
        if key in store:
            return np.array(store[key])
    raise KeyError(f"Array de muestras no encontrado. Candidatos: {candidates}")


def _detect_dosage_key(store: zarr.Group) -> tuple[str, str]:
    """Devuelve (key, layout) donde layout ∈ {'vs', 'sv', 'gt'}."""
    for key in ("call_dosage", "calldata/DS", "dosage"):
        if key in store:
            shape = store[key].shape
            return key, "vs" if shape[0] >= shape[1] else "sv"
    for key in ("call_genotype", "calldata/GT"):
        if key in store:
            return key, "gt"
    raise KeyError(
        "Array de dosage/genotipo no encontrado. "
        "Candidatos: call_dosage, calldata/DS, call_genotype, calldata/GT."
    )


def read_zarr_window(
    store: zarr.Group,
    dosage_key: str,
    layout: str,
    w_start: int,
    w_end: int,
) -> np.ndarray:
    """
    Lee arr[w_start:w_end] del zarr en su dtype nativo.
    Retorna (w_end−w_start, n_samples).

    Slice contiguo = máxima eficiencia zarr: descomprime exactamente los
    chunks internos que cubren el rango, sin tocar ningún otro.
    """
    arr = store[dosage_key]

    if layout == "vs":
        return arr[w_start:w_end]

    if layout == "sv":
        return arr[:, w_start:w_end].T

    if layout == "gt":
        raw = arr[w_start:w_end]  # (W, n_samples, ploidy) dtype nativo int8/int16
        # raw.astype(float64) ANTES de sumar crearía (W, n_samples, 2)×8B → OOM.
        # Orden correcto: marcar missing → sumar en int → convertir a float64 al final.
        missing_mask = (raw < 0).any(axis=2)  # (W, n_samples) bool — 1 B/elemento
        allele_sum = raw.clip(min=0).sum(axis=2)  # (W, n_samples) int  — 2-4 B/elemento
        dosage = allele_sum.astype(np.float16)  # (W, n_samples) float64 — solo aquí
        dosage[missing_mask] = np.nan
        return dosage

    raise ValueError(f"Layout desconocido: {layout!r}")


# =============================================================================
# Cómputo del PRS por ventanas zarr
# =============================================================================


def compute_prs_windowed(
    store: zarr.Group,
    dosage_key: str,
    layout: str,
    n_zarr_variants: int,
    zarr_indices: np.ndarray,  # ordenados ascendente
    betas: np.ndarray,
    flip_mask: np.ndarray,
    n_samples: int,
    window_size: int,
    missing_strategy: str,
    excluded_from_missing: list,
) -> np.ndarray:
    """
    PRS aditivo acumulado iterando sobre ventanas fijas del zarr.

    Loop externo  → ventanas [w, w+W) del eje variantes zarr.
    Loop interno  → solo las variantes del PRS que caen en esa ventana.

    Peak de RAM ≈ W × n_samples × bytes_dtype_fuente
    Ventanas sin matches del PRS se saltan sin ningún I/O al zarr.
    """
    # pdb.set_trace()
    prs = np.zeros(n_samples, dtype=np.float64)
    n_matched = len(zarr_indices)
    n_windows_total = (n_zarr_variants + window_size - 1) // window_size
    n_windows_read = 0
    n_windows_skip = 0
    ptr_lo = 0  # cursor al primer match no procesado

    for w_i in range(n_windows_total):
        w_start = w_i * window_size
        w_end = min(w_start + window_size, n_zarr_variants)

        # Avanzar cursor: descartar matches de ventanas anteriores ya procesadas
        while ptr_lo < n_matched and zarr_indices[ptr_lo] < w_start:
            ptr_lo += 1

        # Encontrar el límite superior de matches en esta ventana
        ptr_hi = ptr_lo
        while ptr_hi < n_matched and zarr_indices[ptr_hi] < w_end:
            ptr_hi += 1

        if ptr_lo == ptr_hi:
            n_windows_skip += 1
            continue  # ← cero I/O: ningún match en esta ventana

        # ── Lectura zarr: un slice contiguo [w_start, w_end) ──────────────────
        window_data = read_zarr_window(store, dosage_key, layout, w_start, w_end)
        n_windows_read += 1
        # pdb.set_trace()

        # ── Acumular contribuciones de los matches en esta ventana ─────────────
        for m in range(ptr_lo, ptr_hi):
            local_idx = int(zarr_indices[m]) - w_start
            d = window_data[local_idx].astype(np.float16)

            if flip_mask[m]:
                d = 2.0 - d

            nan_mask = np.isnan(d)
            n_missing = int(nan_mask.sum())

            if n_missing > 0:
                if missing_strategy == "mean":
                    d = d.copy()
                    d[nan_mask] = float(np.nanmean(d))
                elif missing_strategy == "zero":
                    d = d.copy()
                    d[nan_mask] = 0.0
                elif missing_strategy == "skip":
                    excluded_from_missing.append(
                        {
                            "zarr_index": int(zarr_indices[m]),
                            "n_missing_samples": n_missing,
                            "reason": "skip_por_missing",
                        }
                    )
                    continue

            prs += betas[m] * d

        # window_data sale de scope → GC libera la RAM de la ventana

        if n_windows_read == 0 or n_windows_read % 50 == 0:
            pct = 100.0 * w_end / n_zarr_variants
            log.info(
                "  ventana %d/%d  (%.1f%% zarr)  leídas=%d  saltadas=%d",
                w_i + 1,
                n_windows_total,
                pct,
                n_windows_read,
                n_windows_skip,
            )

    log.info(
        "  Fin: leídas=%d  saltadas=%d  (%.1f%% skip)",
        n_windows_read,
        n_windows_skip,
        100.0 * n_windows_skip / max(n_windows_total, 1),
    )
    return prs


# =============================================================================
# Guardar resultados
# =============================================================================


def save_results(
    pgs_id: str,
    chrom: str,
    samples: np.ndarray,
    prs: np.ndarray,
    metadata: dict,
) -> None:
    out_dir = SCORES_BASE / pgs_id
    out_dir.mkdir(parents=True, exist_ok=True)

    scores_path = out_dir / f"{pgs_id}_chr{chrom}_scores.tsv"
    pd.DataFrame({"sample_id": samples, "PRS": prs}).to_csv(
        scores_path, sep="\t", index=False, float_format="%.8f"
    )
    log.info("Scores guardados en: %s", scores_path)

    meta_path = out_dir / f"{pgs_id}_chr{chrom}_metadata.json"
    with open(meta_path, "w") as fh:
        json.dump(metadata, fh, indent=2, default=str)
    log.info("Metadatos guardados en: %s", meta_path)


# =============================================================================
# Main
# =============================================================================


def main() -> None:
    args = parse_args()
    pgs_id = args.pgs_id
    chrom = args.chrom
    t_start = datetime.now()

    # ── 1. Pesos del PRS ──────────────────────────────────────────────────────
    weights = load_score_weights(pgs_id, chrom)
    if weights.empty:
        log.warning("Sin variantes para chr%s en %s — saltando.", chrom, pgs_id)
        sys.exit(0)

    # ── 2. Matching vía DuckDB ────────────────────────────────────────────────
    # No se cargan los arrays pos/ref/alt del zarr: la consulta retorna solo
    # las filas del PRS (O(n_PRS) << O(n_zarr)).
    log.info("Consultando índice DuckDB: %s", args.index_path)
    z_idx_sorted, betas_sorted, flip_sorted, excluded = match_via_duckdb(
        weights, chrom, args.index_path
    )

    if len(z_idx_sorted) == 0:
        log.error("Cero variantes emparejadas — abortando.")
        sys.exit(1)

    # ── 3. Zarr: apertura y detección (solo dosage) ───────────────────────────
    store = open_zarr_store(chrom)
    samples = get_samples(store)

    dosage_key, layout = _detect_dosage_key(store)
    ax = 0 if layout in ("vs", "gt") else 1
    n_zarr_variants = store[dosage_key].shape[ax]
    n_samples = len(samples)

    log.info(
        "Zarr: %d variantes × %d muestras  |  key='%s'  layout=%s",
        n_zarr_variants,
        n_samples,
        dosage_key,
        layout,
    )
    dtype_bytes = store[dosage_key].dtype.itemsize
    log.info(
        "Window size W=%d  |  Peak RAM estimado: %.2f GB",
        args.window_size,
        args.window_size * n_samples * dtype_bytes / 1e9,
    )

    # ── 4. PRS por ventanas zarr ──────────────────────────────────────────────
    log.info(
        "Calculando PRS: %d matches, W=%d, missing=%s...",
        len(z_idx_sorted),
        args.window_size,
        args.missing_strategy,
    )
    excluded_missing: list[dict] = []
    prs = compute_prs_windowed(
        store=store,
        dosage_key=dosage_key,
        layout=layout,
        n_zarr_variants=n_zarr_variants,
        zarr_indices=z_idx_sorted,
        betas=betas_sorted,
        flip_mask=flip_sorted,
        n_samples=n_samples,
        window_size=args.window_size,
        missing_strategy=args.missing_strategy,
        excluded_from_missing=excluded_missing,
    )

    elapsed = (datetime.now() - t_start).total_seconds()

    # ── 5. Metadatos ──────────────────────────────────────────────────────────
    metadata = {
        "pgs_id": pgs_id,
        "chrom": chrom,
        "run_timestamp": t_start.isoformat(),
        "elapsed_seconds": round(elapsed, 2),
        "missing_strategy": args.missing_strategy,
        "window_size": args.window_size,
        "dosage_array_key": dosage_key,
        "dosage_layout": layout,
        "index_path": str(args.index_path),
        "n_samples": int(n_samples),
        "n_zarr_variants": int(n_zarr_variants),
        "n_variants_in_score_file": int(len(weights)),
        "n_variants_matched": int(len(z_idx_sorted)),
        "n_variants_excluded_allele": int(len(excluded)),
        "n_variants_excluded_missing": int(len(excluded_missing)),
        "prs_mean": float(np.mean(prs)),
        "prs_std": float(np.std(prs)),
        "prs_min": float(np.min(prs)),
        "prs_max": float(np.max(prs)),
        "excluded_allele_mismatch": excluded,
        "excluded_missing": excluded_missing,
    }

    # ── 6. Guardar ────────────────────────────────────────────────────────────
    save_results(pgs_id, chrom, samples, prs, metadata)
    log.info("chr%s completado en %.1f s ✓", chrom, elapsed)


if __name__ == "__main__":
    main()
