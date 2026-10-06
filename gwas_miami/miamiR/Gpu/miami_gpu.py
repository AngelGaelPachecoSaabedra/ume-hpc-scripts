#!/usr/bin/env python3
"""
Miami Plot Generator - GPU Accelerated Version
===============================================
Usa CuPy para cómputo en GPU + Datashader para rasterización eficiente
+ matplotlib solo para decoración (ejes, títulos, anotaciones)

Requisitos (disponibles en BioContainer jupyter-biotools-1.2.sif):
  - Entorno jax-gpu: cupy, numpy, scipy
  - Entorno jupyter: datashader, matplotlib, pandas, seaborn, adjustText
  
Ejecución:
  biocupy miami_gpu.py
  
  O con Slurm:
  sbatch miami_gpu_slurm.sh
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mimage
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from pathlib import Path
import gzip
import time
from typing import Tuple, Optional, Dict, List
import warnings
from collections import defaultdict
import sys
import os

warnings.filterwarnings('ignore')

# ============================================================
# Intentar importar CuPy (GPU). Si no hay GPU, fallback a NumPy
# ============================================================
try:
    import cupy as cp
    GPU_AVAILABLE = True
    print(f"✅ CuPy {cp.__version__} detectado — GPU: {cp.cuda.runtime.getDevice()}")
    mempool = cp.get_default_memory_pool()
except ImportError:
    cp = np  # fallback transparente
    GPU_AVAILABLE = False
    print("⚠️  CuPy no disponible, usando NumPy (CPU fallback)")

# ============================================================
# Intentar importar Datashader
# ============================================================
try:
    import datashader as ds
    import datashader.transfer_functions as tf
    from datashader.mpl_ext import dsshow
    DATASHADER_AVAILABLE = True
    print(f"✅ Datashader {ds.__version__} disponible")
except ImportError:
    DATASHADER_AVAILABLE = False
    print("⚠️  Datashader no disponible, usando matplotlib scatter (más lento)")

try:
    from adjustText import adjust_text
    ADJUSTTEXT_AVAILABLE = True
except ImportError:
    ADJUSTTEXT_AVAILABLE = False
    print("⚠️  adjustText no disponible, anotaciones sin ajuste automático")


# ============================================================
# Configuración de estilo
# ============================================================
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")


class GTFParserGPU:
    """
    Parser GTF optimizado con búsqueda vectorizada en GPU.
    Carga genes una vez, luego anota miles de SNPs en segundos
    usando operaciones matriciales en GPU (CuPy).
    """

    def __init__(self, gtf_path: str):
        self.gtf_path = gtf_path
        self.genes_by_chr = defaultdict(list)  # {chr: [(start, end, gene_name)]}
        # Arrays GPU precalculados por cromosoma
        self.gpu_starts = {}   # {chr: cupy.array}
        self.gpu_ends = {}     # {chr: cupy.array}
        self.gene_names = {}   # {chr: list[str]}
        self.loaded = False

    def load_genes(self) -> bool:
        """Carga genes desde archivo GTF y prepara arrays GPU"""
        if not Path(self.gtf_path).exists():
            print(f"❌ Archivo GTF no encontrado: {self.gtf_path}")
            return False

        print(f"📋 Cargando genes desde {self.gtf_path}...")
        start_time = time.time()

        valid_chroms = set([str(i) for i in range(1, 23)] + ['X', 'Y', 'MT'])
        gene_count = 0

        try:
            with gzip.open(self.gtf_path, 'rt') as f:
                for line_num, line in enumerate(f):
                    if line_num % 200000 == 0 and line_num > 0:
                        print(f"   {line_num:,} líneas procesadas, {gene_count:,} genes")

                    if line.startswith('#'):
                        continue

                    fields = line.strip().split('\t')
                    if len(fields) < 9 or fields[2] != 'gene':
                        continue

                    chrom = fields[0].replace('chr', '')
                    if chrom not in valid_chroms:
                        continue

                    start = int(fields[3])
                    end = int(fields[4])
                    gene_name = self._extract_attr(fields[8], 'gene_name')

                    if gene_name:
                        self.genes_by_chr[chrom].append((start, end, gene_name))
                        gene_count += 1

            # Ordenar y crear arrays GPU por cromosoma
            self._build_gpu_arrays()

            elapsed = time.time() - start_time
            print(f"✅ GTF cargado en {elapsed:.2f}s: {gene_count:,} genes, "
                  f"{len(self.genes_by_chr)} cromosomas")
            if GPU_AVAILABLE:
                print(f"   GPU arrays creados para búsqueda vectorizada")
            self.loaded = True
            return True

        except Exception as e:
            print(f"❌ Error cargando GTF: {e}")
            import traceback
            traceback.print_exc()
            return False

    def _extract_attr(self, attributes: str, key: str) -> Optional[str]:
        """Extrae un atributo del campo de atributos GTF"""
        for attr in attributes.split(';'):
            attr = attr.strip()
            if attr.startswith(key):
                parts = attr.split('"')
                if len(parts) >= 2:
                    return parts[1]
        return None

    def _build_gpu_arrays(self):
        """Construye arrays en GPU para búsqueda vectorizada"""
        for chrom in self.genes_by_chr:
            genes = sorted(self.genes_by_chr[chrom], key=lambda x: x[0])
            self.genes_by_chr[chrom] = genes

            starts = np.array([g[0] for g in genes], dtype=np.int64)
            ends = np.array([g[1] for g in genes], dtype=np.int64)
            names = [g[2] for g in genes]

            if GPU_AVAILABLE:
                self.gpu_starts[chrom] = cp.asarray(starts)
                self.gpu_ends[chrom] = cp.asarray(ends)
            else:
                self.gpu_starts[chrom] = starts
                self.gpu_ends[chrom] = ends

            self.gene_names[chrom] = names

    def annotate_snps_vectorized(self, snps_df: pd.DataFrame) -> pd.DataFrame:
        """
        Anota TODOS los SNPs de un cromosoma de forma vectorizada en GPU.
        En vez de iterar SNP por SNP, calcula distancias con operaciones matriciales.
        
        Para cromosomas con muchos genes (>50k), usa batching para no exceder VRAM.
        """
        if not self.loaded:
            print("⚠️  GTF no cargado")
            return self._add_empty_annotations(snps_df)

        print(f"🧬 Anotando {len(snps_df):,} SNPs con {'GPU' if GPU_AVAILABLE else 'CPU'} vectorizado...")
        start_time = time.time()

        all_gene_names = []
        all_distances = []

        chromosomes = snps_df['CHR'].unique()

        for chrom in chromosomes:
            chrom_str = str(chrom)
            mask = snps_df['CHR'] == chrom
            snp_positions = snps_df.loc[mask, 'POS'].values

            if chrom_str not in self.gpu_starts or len(snp_positions) == 0:
                all_gene_names.extend([f"UNKNOWN_chr{chrom}"] * len(snp_positions))
                all_distances.extend([-1] * len(snp_positions))
                continue

            gene_starts = self.gpu_starts[chrom_str]
            gene_ends = self.gpu_ends[chrom_str]
            names = self.gene_names[chrom_str]

            n_snps = len(snp_positions)
            n_genes = len(names)

            print(f"   Chr{chrom_str}: {n_snps:,} SNPs × {n_genes:,} genes", end="")

            # Decidir si necesitamos batching
            # Cada par SNP-gen usa ~8 bytes (int64), la matriz es n_snps × n_genes
            matrix_size_bytes = n_snps * n_genes * 8
            max_gpu_matrix = 2 * 1024**3  # 2 GB límite por batch

            if matrix_size_bytes > max_gpu_matrix:
                # Procesar en batches
                batch_size = max(1, int(max_gpu_matrix / (n_genes * 8)))
                print(f" [batched: {batch_size} SNPs/batch]")
                g_names, g_dists = self._annotate_batched(
                    snp_positions, gene_starts, gene_ends, names, batch_size
                )
            else:
                print(f" [vectorizado completo]")
                g_names, g_dists = self._annotate_full_vectorized(
                    snp_positions, gene_starts, gene_ends, names
                )

            all_gene_names.extend(g_names)
            all_distances.extend(g_dists)

        # Reconstruir DataFrame respetando el orden original
        result_df = snps_df.copy()
        result_df['nearest_gene'] = all_gene_names
        result_df['distance_to_gene'] = all_distances

        elapsed = time.time() - start_time
        annotated = sum(1 for g in all_gene_names if not g.startswith('UNKNOWN'))
        within = sum(1 for d in all_distances if d == 0)

        print(f"✅ Anotación completada en {elapsed:.2f}s")
        print(f"   Anotados: {annotated}/{len(snps_df)} "
              f"({annotated/max(len(snps_df),1)*100:.1f}%)")
        print(f"   Dentro de genes: {within} ({within/max(len(snps_df),1)*100:.1f}%)")

        if GPU_AVAILABLE:
            mempool.free_all_blocks()

        return result_df

    def _annotate_full_vectorized(self, snp_positions, gene_starts, gene_ends, names):
        """Anotación vectorizada completa — toda la matriz en GPU de una vez"""
        xp = cp if GPU_AVAILABLE else np
        snp_pos = xp.asarray(snp_positions, dtype=xp.int64)

        # Matrices de distancia: (n_snps, n_genes)
        # snp_pos[:, None] -> columna, gene_starts[None, :] -> fila
        pos_col = snp_pos[:, None]  # (N, 1)
        starts_row = gene_starts[None, :]  # (1, M)
        ends_row = gene_ends[None, :]  # (1, M)

        # Distancia: 0 si dentro del gen, positiva si fuera
        dist_upstream = starts_row - pos_col  # positivo si SNP está antes del gen
        dist_downstream = pos_col - ends_row  # positivo si SNP está después del gen
        inside = (pos_col >= starts_row) & (pos_col <= ends_row)

        # Distancia final: max(dist_upstream, dist_downstream, 0) pero 0 si inside
        distances = xp.where(
            inside, 0,
            xp.where(
                dist_upstream > 0, dist_upstream,
                xp.where(dist_downstream > 0, dist_downstream, 0)
            )
        )

        # Gen más cercano por SNP
        nearest_idx = xp.argmin(distances, axis=1)
        min_distances = distances[xp.arange(len(snp_pos)), nearest_idx]

        # Traer a CPU
        if GPU_AVAILABLE:
            nearest_idx_cpu = cp.asnumpy(nearest_idx)
            min_distances_cpu = cp.asnumpy(min_distances)
        else:
            nearest_idx_cpu = nearest_idx
            min_distances_cpu = min_distances

        result_names = [names[i] for i in nearest_idx_cpu]
        result_dists = min_distances_cpu.tolist()

        return result_names, result_dists

    def _annotate_batched(self, snp_positions, gene_starts, gene_ends, names, batch_size):
        """Anotación por batches cuando la matriz es demasiado grande para GPU"""
        all_names = []
        all_dists = []

        n_snps = len(snp_positions)
        for i in range(0, n_snps, batch_size):
            batch_pos = snp_positions[i:i+batch_size]
            b_names, b_dists = self._annotate_full_vectorized(
                batch_pos, gene_starts, gene_ends, names
            )
            all_names.extend(b_names)
            all_dists.extend(b_dists)

            if GPU_AVAILABLE:
                mempool.free_all_blocks()

        return all_names, all_dists

    def _add_empty_annotations(self, df):
        df = df.copy()
        df['nearest_gene'] = 'UNKNOWN'
        df['distance_to_gene'] = -1
        return df


class MiamiPlotGPU:
    """
    Generador de Miami Plot acelerado por GPU.
    
    Pipeline:
    1. Carga datos GWAS (pandas, con chunks para archivos grandes)
    2. Transformaciones numéricas en GPU (CuPy)
    3. Anotación GTF vectorizada en GPU
    4. Rasterización con Datashader (millones de puntos → imagen)
    5. Decoración con matplotlib (ejes, títulos, líneas, anotaciones)
    """

    def __init__(self, gtf_path: str = "Homo_sapiens.GRCh38.105.gtf.gz"):
        self.gtf_path = gtf_path
        self.chr_labels = [str(i) for i in range(1, 23)] + ['X']
        self.colors_male = ['#4A90E2', '#7B68EE']    # Azules alternos
        self.colors_female = ['#E24A6A', '#EE687B']   # Rojos alternos
        self.gtf_parser = GTFParserGPU(gtf_path)

    # ================================================================
    # CARGA DE DATOS
    # ================================================================
    def load_gwas_data(self, filepath: str) -> pd.DataFrame:
        """Carga datos GWAS completos con chunks optimizados"""
        print(f"📂 Cargando datos de {filepath}...")
        start_time = time.time()

        chunk_size = 200_000
        chunks = []
        total_rows = 0
        required_cols = ['CHROM', 'GENPOS', 'P']

        try:
            # Verificar columnas
            if filepath.endswith('.gz'):
                with gzip.open(filepath, 'rt') as f:
                    header = f.readline().strip().split('\t')
            else:
                with open(filepath, 'rt') as f:
                    header = f.readline().strip().split('\t')

            missing = [c for c in required_cols if c not in header]
            if missing:
                print(f"❌ Columnas faltantes: {missing}")
                print(f"   Columnas disponibles: {header}")
                return pd.DataFrame()

            # Leer por chunks
            for i, chunk in enumerate(pd.read_csv(filepath, sep='\t',
                                                   chunksize=chunk_size,
                                                   usecols=required_cols)):
                processed = self._process_chunk(chunk)
                if len(processed) > 0:
                    chunks.append(processed)
                    total_rows += len(processed)

                if (i + 1) % 5 == 0:
                    print(f"   Chunk {i+1}: {total_rows:,} SNPs válidos")

        except Exception as e:
            print(f"❌ Error cargando {filepath}: {e}")
            import traceback
            traceback.print_exc()
            return pd.DataFrame()

        if not chunks:
            print(f"⚠️  No se encontraron datos válidos")
            return pd.DataFrame()

        df = pd.concat(chunks, ignore_index=True)
        df = self._optimize_dtypes(df)

        elapsed = time.time() - start_time
        mem_mb = df.memory_usage(deep=True).sum() / 1024**2
        print(f"✅ Cargado en {elapsed:.2f}s: {len(df):,} SNPs ({mem_mb:.1f} MB)")
        return df

    def _process_chunk(self, chunk: pd.DataFrame) -> pd.DataFrame:
        """Procesa un chunk de datos GWAS"""
        chunk = chunk.dropna(subset=['CHROM', 'GENPOS', 'P'])
        chunk['CHROM'] = chunk['CHROM'].astype(str).replace('23', 'X')

        valid_chrs = [str(i) for i in range(1, 23)] + ['23', 'X']
        chunk = chunk[chunk['CHROM'].isin(valid_chrs)]
        if len(chunk) == 0:
            return pd.DataFrame()

        chunk['CHR'] = pd.Categorical(chunk['CHROM'], categories=self.chr_labels)
        chunk['POS'] = pd.to_numeric(chunk['GENPOS'], errors='coerce')
        chunk['pvalue'] = pd.to_numeric(chunk['P'], errors='coerce')
        chunk = chunk.dropna(subset=['POS', 'pvalue'])
        chunk = chunk[chunk['pvalue'] > 0]

        if len(chunk) == 0:
            return pd.DataFrame()

        chunk['SNP'] = 'rs' + chunk['POS'].astype(int).astype(str)
        return chunk[['SNP', 'CHR', 'POS', 'pvalue']].copy()

    def _optimize_dtypes(self, df: pd.DataFrame) -> pd.DataFrame:
        """Optimiza tipos de datos para memoria"""
        if df['POS'].max() < 2**31:
            df['POS'] = df['POS'].astype('int32')
        if df['pvalue'].min() > 1e-38:
            df['pvalue'] = df['pvalue'].astype('float32')
        return df

    # ================================================================
    # TRANSFORMACIONES GPU
    # ================================================================
    def compute_plot_coordinates_gpu(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Calcula coordenadas del plot en GPU:
        - Posiciones acumulativas por cromosoma
        - -log10(p-value)
        """
        print(f"🎮 Calculando coordenadas con {'GPU' if GPU_AVAILABLE else 'CPU'}...")
        start_time = time.time()

        df = df.copy()
        xp = cp if GPU_AVAILABLE else np

        # 1. Calcular offsets acumulativos por cromosoma
        chr_max_pos = {}
        for chr_name in self.chr_labels:
            chr_data = df[df['CHR'] == chr_name]
            chr_max_pos[chr_name] = chr_data['POS'].max() if len(chr_data) > 0 else 0

        cumulative_offset = 0
        chr_offsets = {}
        chr_centers = {}
        for chr_name in self.chr_labels:
            chr_offsets[chr_name] = cumulative_offset
            chr_data = df[df['CHR'] == chr_name]
            if len(chr_data) > 0:
                chr_centers[chr_name] = cumulative_offset + chr_data['POS'].mean()
            else:
                chr_centers[chr_name] = cumulative_offset
            cumulative_offset += chr_max_pos[chr_name] + 10_000_000  # gap entre cromosomas

        # 2. Aplicar transformaciones en GPU/CPU
        positions = xp.asarray(df['POS'].values, dtype=xp.float64)
        pvalues = xp.asarray(df['pvalue'].values, dtype=xp.float64)

        # -log10(p)
        log10p = -xp.log10(pvalues)

        # Posiciones acumulativas (mapear offsets por cromosoma)
        offsets_array = xp.zeros(len(df), dtype=xp.float64)
        for chr_name in self.chr_labels:
            mask = (df['CHR'] == chr_name).values
            if GPU_AVAILABLE:
                mask_gpu = cp.asarray(mask)
                offsets_array[mask_gpu] = chr_offsets[chr_name]
            else:
                offsets_array[mask] = chr_offsets[chr_name]

        pos_cumulative = positions + offsets_array

        # Traer a CPU para pandas
        if GPU_AVAILABLE:
            df['pos_cumulative'] = cp.asnumpy(pos_cumulative)
            df['-log10p'] = cp.asnumpy(log10p)
            mempool.free_all_blocks()
        else:
            df['pos_cumulative'] = pos_cumulative
            df['-log10p'] = log10p

        # Color index por cromosoma (para colores alternos)
        df['chr_idx'] = df['CHR'].cat.codes

        elapsed = time.time() - start_time
        print(f"   Coordenadas calculadas en {elapsed:.2f}s")

        return df, chr_centers

    # ================================================================
    # BÚSQUEDA DE SNPs SIGNIFICATIVOS
    # ================================================================
    def find_significant_snps(self, df_males: pd.DataFrame, df_females: pd.DataFrame,
                              threshold: float = 5e-8) -> pd.DataFrame:
        """Encuentra todos los SNPs significativos"""
        print(f"🔍 Buscando SNPs significativos (p < {threshold})...")
        males_sig = df_males[df_males['pvalue'] < threshold].copy()
        males_sig['dataset'] = 'Males'
        females_sig = df_females[df_females['pvalue'] < threshold].copy()
        females_sig['dataset'] = 'Females'

        all_sig = pd.concat([males_sig, females_sig], ignore_index=True)
        print(f"   Males: {len(males_sig):,}, Females: {len(females_sig):,}, "
              f"Total: {len(all_sig):,}")
        return all_sig

    def select_top_snp_per_gene(self, annotated_snps: pd.DataFrame,
                                max_distance: int = 500_000) -> pd.DataFrame:
        """
        Selecciona el SNP más significativo por gen, separado por sexo.
        Filtra SNPs que estén a más de max_distance del gen.
        """
        if annotated_snps.empty:
            return pd.DataFrame()

        # Filtrar por distancia
        in_range = annotated_snps[
            (annotated_snps['distance_to_gene'] >= 0) &
            (annotated_snps['distance_to_gene'] <= max_distance)
        ]

        if in_range.empty:
            return pd.DataFrame()

        # Top SNP por gen y sexo
        idx = in_range.groupby(['nearest_gene', 'dataset'])['pvalue'].idxmin()
        top = in_range.loc[idx].reset_index(drop=True)

        print(f"   SNPs únicos por gen: {len(top)} "
              f"(Males: {(top['dataset']=='Males').sum()}, "
              f"Females: {(top['dataset']=='Females').sum()})")
        return top

    def select_best_for_annotation(self, top_snps: pd.DataFrame,
                                    max_males: int = 50,
                                    max_females: int = 50) -> pd.DataFrame:
        """Selecciona los N mejores SNPs para anotar en el plot"""
        if top_snps.empty:
            return pd.DataFrame()

        males = top_snps[top_snps['dataset'] == 'Males'].nsmallest(max_males, 'pvalue')
        females = top_snps[top_snps['dataset'] == 'Females'].nsmallest(max_females, 'pvalue')

        selected = pd.concat([males, females], ignore_index=True)
        return selected.sort_values('pvalue').reset_index(drop=True)

    # ================================================================
    # MIAMI PLOT — DATASHADER + MATPLOTLIB
    # ================================================================
    def create_miami_plot(self, df_males: pd.DataFrame, df_females: pd.DataFrame,
                          chr_centers_males: dict, chr_centers_females: dict,
                          top_snps: pd.DataFrame = None,
                          output_file: str = "miami_plot_gpu.png",
                          figsize: Tuple[int, int] = (18, 12),
                          max_annotations_males: int = 50,
                          max_annotations_females: int = 50) -> None:
        """
        Genera Miami plot:
        - Datashader para rasterizar millones de puntos (rápido)
        - matplotlib para ejes, títulos, líneas de significancia, anotaciones
        """
        print("📊 Generando Miami plot...")
        start_time = time.time()

        # Offset vertical para separar del eje central
        y_offset = 5

        fig, ax = plt.subplots(1, 1, figsize=figsize, facecolor='white')

        # ============================================================
        # Preparar datos para Miami (males arriba, females abajo)
        # ============================================================
        males_plot = df_males.copy()
        males_plot['y'] = males_plot['-log10p'] + y_offset

        females_plot = df_females.copy()
        females_plot['y'] = -(females_plot['-log10p'] + y_offset)

        # Color por cromosoma (0 o 1 para colores alternos)
        males_plot['color_cat'] = males_plot['chr_idx'] % 2
        females_plot['color_cat'] = females_plot['chr_idx'] % 2

        # ============================================================
        # DATASHADER o MATPLOTLIB SCATTER
        # ============================================================
        if DATASHADER_AVAILABLE:
            self._render_datashader(ax, males_plot, females_plot, figsize)
        else:
            self._render_matplotlib_scatter(ax, males_plot, females_plot, y_offset)

        # ============================================================
        # Líneas de significancia
        # ============================================================
        sig_line = -np.log10(5e-8) + y_offset
        sug_line = -np.log10(1e-5) + y_offset

        # Males (arriba)
        ax.axhline(y=sig_line, color='red', linestyle='--', alpha=0.8, linewidth=1.5,
                    label='Genome-wide (5e-8)')
        ax.axhline(y=sug_line, color='blue', linestyle=':', alpha=0.5, linewidth=1)
        # Females (abajo)
        ax.axhline(y=-sig_line, color='red', linestyle='--', alpha=0.8, linewidth=1.5)
        ax.axhline(y=-sug_line, color='blue', linestyle=':', alpha=0.5, linewidth=1)

        # ============================================================
        # Línea central
        # ============================================================
        ax.axhline(y=0, color='black', linewidth=2, alpha=0.9)

        # ============================================================
        # Anotaciones de genes
        # ============================================================
        if top_snps is not None and len(top_snps) > 0:
            self._add_gene_annotations(ax, top_snps, males_plot, females_plot, y_offset)

        # ============================================================
        # Configurar ejes
        # ============================================================
        self._configure_axes(ax, males_plot, females_plot, chr_centers_males, y_offset)

        # ============================================================
        # Etiquetas y título
        # ============================================================
        ax.text(0.02, 0.88, 'MALES', transform=ax.transAxes,
                fontsize=16, fontweight='bold', color='darkblue',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue', alpha=0.7))

        ax.text(0.02, 0.05, 'FEMALES', transform=ax.transAxes,
                fontsize=16, fontweight='bold', color='darkred',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='lightpink', alpha=0.7))

        ax.set_title('Miami Plot — GWAS Results\nMales vs Females (GPU Accelerated)',
                      fontsize=18, fontweight='bold', pad=20)

        # Guardar
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"💾 Guardando {output_file}...")

        # También guardar versión SVG para mayor calidad
        svg_file = output_file.replace('.png', '.svg')
        try:
            plt.savefig(svg_file, format='svg', bbox_inches='tight', facecolor='white')
            print(f"💾 SVG guardado: {svg_file}")
        except Exception:
            pass

        plt.close()

        elapsed = time.time() - start_time
        print(f"✅ Miami plot generado en {elapsed:.2f}s → {output_file}")

    def _render_datashader(self, ax, males_plot, females_plot, figsize):
        """Rasteriza puntos con Datashader — extremadamente rápido para millones de puntos"""
        print("   🎯 Rasterizando con Datashader...")

        x_range = (
            min(males_plot['pos_cumulative'].min(), females_plot['pos_cumulative'].min()),
            max(males_plot['pos_cumulative'].max(), females_plot['pos_cumulative'].max())
        )
        y_range = (
            females_plot['y'].min() * 1.1,
            males_plot['y'].max() * 1.1
        )

        plot_width = int(figsize[0] * 100)   # píxeles
        plot_height = int(figsize[1] * 100)

        canvas = ds.Canvas(
            plot_width=plot_width,
            plot_height=plot_height,
            x_range=x_range,
            y_range=y_range
        )

        # Males: azules
        males_df = males_plot[['pos_cumulative', 'y', 'color_cat']].copy()
        males_df.columns = ['x', 'y', 'cat']
        males_df['cat'] = males_df['cat'].astype('category')

        agg_m = canvas.points(males_df, 'x', 'y', ds.count_cat('cat'))
        img_m = tf.shade(agg_m, color_key={0: '#4A90E2', 1: '#7B68EE'}, how='log')

        # Females: rojos
        females_df = females_plot[['pos_cumulative', 'y', 'color_cat']].copy()
        females_df.columns = ['x', 'y', 'cat']
        females_df['cat'] = females_df['cat'].astype('category')

        agg_f = canvas.points(females_df, 'x', 'y', ds.count_cat('cat'))
        img_f = tf.shade(agg_f, color_key={0: '#E24A6A', 1: '#EE687B'}, how='log')

        # Combinar imágenes
        combined = tf.stack(img_m, img_f)

        # Convertir a numpy array para matplotlib
        img_array = np.array(combined.to_pil())

        ax.imshow(img_array, extent=[x_range[0], x_range[1], y_range[0], y_range[1]],
                  aspect='auto', origin='lower', interpolation='nearest')

    def _render_matplotlib_scatter(self, ax, males_plot, females_plot, y_offset):
        """Fallback: matplotlib scatter con rasterizado (más lento)"""
        print("   📊 Renderizando con matplotlib scatter (fallback)...")

        for i, chr_name in enumerate(self.chr_labels):
            # Males
            chr_m = males_plot[males_plot['CHR'] == chr_name]
            if len(chr_m) > 0:
                color = self.colors_male[i % 2]
                alpha = min(0.8, max(0.1, 10000 / max(len(chr_m), 1)))
                ax.scatter(chr_m['pos_cumulative'], chr_m['y'],
                          c=color, s=1.5, alpha=alpha, rasterized=True)

            # Females
            chr_f = females_plot[females_plot['CHR'] == chr_name]
            if len(chr_f) > 0:
                color = self.colors_female[i % 2]
                alpha = min(0.8, max(0.1, 10000 / max(len(chr_f), 1)))
                ax.scatter(chr_f['pos_cumulative'], chr_f['y'],
                          c=color, s=1.5, alpha=alpha, rasterized=True)

    def _add_gene_annotations(self, ax, top_snps, males_data, females_data, y_offset):
        """Añade anotaciones de genes usando adjustText para evitar solapamientos"""
        print(f"   📝 Añadiendo {len(top_snps)} anotaciones de genes...")

        texts = []
        males_snps = top_snps[top_snps['dataset'] == 'Males']
        females_snps = top_snps[top_snps['dataset'] == 'Females']

        for dataset_snps, data, is_male in [
            (males_snps, males_data, True),
            (females_snps, females_data, False)
        ]:
            if dataset_snps.empty:
                continue

            color = 'darkblue' if is_male else 'darkred'

            for _, snp in dataset_snps.iterrows():
                match = data[(data['CHR'] == snp['CHR']) & (data['POS'] == snp['POS'])]
                if len(match) == 0:
                    continue

                x = match.iloc[0]['pos_cumulative']
                y_base = match.iloc[0]['-log10p']

                if is_male:
                    y = y_base + y_offset
                else:
                    y = -(y_base + y_offset)

                gene_label = snp['nearest_gene']

                text = ax.text(
                    x, y, gene_label,
                    fontsize=5.5,
                    fontweight='bold',
                    color=color,
                    ha='center',
                    va='bottom' if is_male else 'top'
                )
                texts.append(text)

        if texts and ADJUSTTEXT_AVAILABLE:
            adjust_text(
                texts, ax=ax,
                arrowprops=dict(
                    arrowstyle='-', color='gray', alpha=0.6,
                    shrinkA=3, shrinkB=3, linewidth=0.8
                ),
                expand_text=(1.5, 1.5),
                expand_points=(2, 2),
                force_text=0.7,
                force_points=0.3,
                only_move={'text': 'xy'},
                autoalign=True,
                precision=0.01
            )

    def _configure_axes(self, ax, males_data, females_data, chr_centers, y_offset):
        """Configura ejes estilo Miami"""
        # Eje X: cromosomas
        chr_positions = [chr_centers.get(c, 0) for c in self.chr_labels]
        ax.set_xticks(chr_positions)
        ax.set_xticklabels(self.chr_labels, fontsize=11, fontweight='bold')

        # Eje Y: simétrico
        max_y_m = males_data['y'].max() if len(males_data) > 0 else 15
        max_y_f = abs(females_data['y'].min()) if len(females_data) > 0 else 15
        y_lim = max(max_y_m, max_y_f) * 1.15
        ax.set_ylim(-y_lim, y_lim)

        # Etiquetas Y absolutas
        y_ticks = ax.get_yticks()
        ax.set_yticklabels([f'{abs(t):.0f}' for t in y_ticks])
        ax.set_ylabel('-log₁₀(p-value)', fontsize=14, fontweight='bold')

        # Mover eje X al centro
        ax.spines['bottom'].set_position('zero')
        ax.xaxis.set_ticks_position('bottom')

        # Estilo
        ax.grid(True, alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)


# ================================================================
# MAIN
# ================================================================
def main():
    """Pipeline principal"""

    # ── Archivos de entrada ──
    male_file = "combined-gwas-results-all-ctrl-males.txt.gz"
    female_file = "combined-gwas-results-all-ctrl-females.txt.gz"
    gtf_file = "Homo_sapiens.GRCh38.105.gtf.gz"

    # ── Parámetros ──
    SIGNIFICANCE_THRESHOLD = 5e-8
    MAX_GENE_DISTANCE = 500_000       # bp
    MAX_ANNOTATIONS_MALES = 50
    MAX_ANNOTATIONS_FEMALES = 50
    OUTPUT_PLOT = "miami_plot_gpu.png"
    OUTPUT_CSV = "top_snps_gpu_annotated.csv"

    # ── Verificar archivos ──
    for f in [male_file, female_file, gtf_file]:
        if not Path(f).exists():
            print(f"❌ Archivo no encontrado: {f}")
            print(f"   Directorio actual: {os.getcwd()}")
            print(f"   Archivos disponibles: {list(Path('.').glob('*'))}")
            sys.exit(1)

    total_start = time.time()
    print("=" * 70)
    print("🚀 MIAMI PLOT GENERATOR — GPU ACCELERATED")
    print(f"   GPU: {'✅ CuPy' if GPU_AVAILABLE else '❌ CPU fallback'}")
    print(f"   Datashader: {'✅' if DATASHADER_AVAILABLE else '❌ matplotlib fallback'}")
    print("=" * 70)

    # ── Inicializar ──
    generator = MiamiPlotGPU(gtf_file)

    # ── 1. Cargar datos GWAS ──
    print("\n" + "─" * 50)
    print("PASO 1: Carga de datos GWAS")
    print("─" * 50)
    df_males = generator.load_gwas_data(male_file)
    df_females = generator.load_gwas_data(female_file)

    if df_males.empty or df_females.empty:
        print("❌ Error: No se pudieron cargar los datos")
        sys.exit(1)

    # ── 2. Calcular coordenadas en GPU ──
    print("\n" + "─" * 50)
    print("PASO 2: Cálculo de coordenadas (GPU)")
    print("─" * 50)
    df_males, chr_centers_m = generator.compute_plot_coordinates_gpu(df_males)
    df_females, chr_centers_f = generator.compute_plot_coordinates_gpu(df_females)

    # ── 3. Encontrar SNPs significativos ──
    print("\n" + "─" * 50)
    print("PASO 3: Identificación de SNPs significativos")
    print("─" * 50)
    all_sig = generator.find_significant_snps(df_males, df_females, SIGNIFICANCE_THRESHOLD)

    # ── 4. Anotar genes con GTF (GPU vectorizado) ──
    top_snps_for_plot = pd.DataFrame()

    if not all_sig.empty:
        print("\n" + "─" * 50)
        print("PASO 4: Anotación GTF (GPU vectorizado)")
        print("─" * 50)

        # Cargar GTF
        generator.gtf_parser.load_genes()

        # Anotar SNPs significativos
        annotated = generator.gtf_parser.annotate_snps_vectorized(all_sig)

        # Top SNP por gen
        top_snps = generator.select_top_snp_per_gene(annotated, MAX_GENE_DISTANCE)

        if not top_snps.empty:
            # Seleccionar mejores para anotación
            top_snps_for_plot = generator.select_best_for_annotation(
                top_snps,
                max_males=MAX_ANNOTATIONS_MALES,
                max_females=MAX_ANNOTATIONS_FEMALES
            )

            # Guardar CSV
            top_snps.to_csv(OUTPUT_CSV, index=False)
            print(f"💾 CSV guardado: {OUTPUT_CSV}")

            # Mostrar tabla
            print(f"\n🧬 Top SNPs para anotación ({len(top_snps_for_plot)}):")
            display_cols = ['CHR', 'POS', 'pvalue', 'nearest_gene', 'distance_to_gene', 'dataset']
            available_cols = [c for c in display_cols if c in top_snps_for_plot.columns]
            print(top_snps_for_plot[available_cols].to_string(index=False))
    else:
        print("\n⚠️  No se encontraron SNPs significativos")

    # ── 5. Generar Miami plot ──
    print("\n" + "─" * 50)
    print("PASO 5: Generación del Miami Plot")
    print("─" * 50)
    generator.create_miami_plot(
        df_males, df_females,
        chr_centers_m, chr_centers_f,
        top_snps=top_snps_for_plot if not top_snps_for_plot.empty else None,
        output_file=OUTPUT_PLOT,
        figsize=(18, 12),
        max_annotations_males=MAX_ANNOTATIONS_MALES,
        max_annotations_females=MAX_ANNOTATIONS_FEMALES
    )

    # ── Resumen final ──
    total_elapsed = time.time() - total_start
    print("\n" + "=" * 70)
    print("📋 RESUMEN FINAL")
    print("=" * 70)
    print(f"  Aceleración:     {'GPU (CuPy)' if GPU_AVAILABLE else 'CPU'}")
    print(f"  Rasterización:   {'Datashader' if DATASHADER_AVAILABLE else 'matplotlib'}")
    print(f"  Males:           {len(df_males):,} SNPs")
    print(f"  Females:         {len(df_females):,} SNPs")
    print(f"  Total SNPs:      {len(df_males) + len(df_females):,}")
    print(f"  Significativos:  {len(all_sig):,}")
    print(f"  Genes anotados:  {len(top_snps_for_plot) if not top_snps_for_plot.empty else 0}")
    print(f"  Plot:            {OUTPUT_PLOT}")
    print(f"  CSV:             {OUTPUT_CSV}")
    print(f"  Tiempo total:    {total_elapsed:.2f}s ({total_elapsed/60:.1f} min)")
    print("=" * 70)

    # Estadísticas de GPU si disponible
    if GPU_AVAILABLE:
        print(f"\n🎮 GPU Memory Peak: {mempool.total_bytes() / 1024**2:.1f} MB")
        mempool.free_all_blocks()


if __name__ == "__main__":
    main()
