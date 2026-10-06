#!/usr/bin/env python3
"""
Miami Plot — FASE 1: Cómputo GPU
================================================
Fixes:
  - Overflow en cumulative_offset (float64).
  - Ordenamiento estricto por CHR/POS.
"""

import numpy as np
import gzip
import time
import sys
import os
from collections import defaultdict

# CuPy GPU
try:
    import cupy as cp
    GPU = True
    print(f"CuPy {cp.__version__} — GPU device {cp.cuda.runtime.getDevice()}")
except Exception as e:
    cp = np
    GPU = False
    print(f"uPy no disponible ({e}), usando NumPy CPU")


def get_chr_order(chr_array):
    chr_map = {str(i): i for i in range(1, 23)}
    chr_map.update({'X': 23, 'Y': 24, 'MT': 25, 'M': 25, '23': 23})
    numeric_chrs = np.array([chr_map.get(str(c).replace('chr', ''), 99) for c in chr_array])
    return numeric_chrs

def load_gwas_numpy(filepath: str) -> dict:
    print(f"Cargando {filepath}...")
    t0 = time.time()
    opener = gzip.open if filepath.endswith('.gz') else open

    chroms = []
    positions = []
    pvalues = []
    
    try:
        with opener(filepath, 'rt') as f:
            header = f.readline().strip().split('\t')
            idx_chr = header.index('CHROM')
            idx_pos = header.index('GENPOS')
            idx_p = header.index('P')
            valid_chrs = set([str(i) for i in range(1, 23)] + ['23', 'X', 'Y', 'MT'])

            for line in f:
                fields = line.strip().split('\t')
                if len(fields) <= max(idx_chr, idx_pos, idx_p): continue
                c = fields[idx_chr].replace('chr', '')
                if c not in valid_chrs: continue
                try:
                    p = float(fields[idx_p])
                    pos = int(fields[idx_pos])
                    if p <= 0: continue
                except: continue
                if c == '23': c = 'X'
                chroms.append(c)
                positions.append(pos)
                pvalues.append(p)
    except ValueError as e:
        print(f"Error columnas: {e}")
        return {}

    data = {
        'CHR': np.array(chroms, dtype='U2'),
        'POS': np.array(positions, dtype=np.int32),
        'pvalue': np.array(pvalues, dtype=np.float64),
    }
    
    print("Ordenando datos por CHR y POS...")
    chr_nums = get_chr_order(data['CHR'])
    sort_idx = np.lexsort((data['POS'], chr_nums))
    
    data['CHR'] = data['CHR'][sort_idx]
    data['POS'] = data['POS'][sort_idx]
    data['pvalue'] = data['pvalue'][sort_idx]
    
    print(f"{len(data['POS']):,} SNPs cargados y ordenados en {time.time()-t0:.1f}s")
    return data


def compute_coordinates(data: dict, chr_labels: list) -> dict:
    print(f"🎮 Calculando coordenadas...")
    xp = cp if GPU else np

    pvals = xp.asarray(data['pvalue'])
    log10p = -xp.log10(pvals)

    # Convertir a float64 para evitar overflow ---
    pos = data['POS'].astype(np.float64) 
    chrs = data['CHR']

    cumulative_offset = 0.0 # Float explícito
    chr_centers = {}
    gap = 10_000_000.0

    pos_cumulative = np.zeros(len(pos), dtype=np.float64)

    for chr_name in chr_labels:
        mask = (chrs == chr_name)
        if np.any(mask):
            current_pos = pos[mask]
            min_p, max_p = np.min(current_pos), np.max(current_pos)
            chr_centers[chr_name] = cumulative_offset + (min_p + max_p) / 2
            pos_cumulative[mask] = current_pos + cumulative_offset
            cumulative_offset += max_p + gap
        else:
            chr_centers[chr_name] = cumulative_offset

    if GPU: log10p = cp.asnumpy(log10p)

    data['log10p'] = log10p
    data['pos_cumulative'] = pos_cumulative
    data['chr_centers'] = chr_centers
    return data

class GTFAnnotatorGPU:
    def __init__(self, gtf_path: str):
        self.gtf_path = gtf_path
        self.genes_by_chr = defaultdict(list)
        self.loaded = False

    def load(self) -> bool:
        if not os.path.exists(self.gtf_path): return False
        print(f"Cargando GTF: {self.gtf_path}")
        t0 = time.time()
        with gzip.open(self.gtf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'): continue
                fields = line.strip().split('\t')
                if len(fields) < 9 or fields[2] != 'gene': continue
                c = fields[0].replace('chr', '')
                try:
                    s, e = int(fields[3]), int(fields[4])
                    attr = fields[8]
                    name_idx = attr.find('gene_name "')
                    if name_idx != -1:
                        name = attr[name_idx+11:].split('"')[0]
                        self.genes_by_chr[c].append((s, e, name))
                except: continue
        for c in self.genes_by_chr: self.genes_by_chr[c].sort(key=lambda x: x[0])
        self.loaded = True
        print(f"Genes cargados en {time.time()-t0:.1f}s")
        return True

    def annotate(self, chrs, positions):
        if not self.loaded: return None, None
        print(f"Anotando {len(chrs):,} SNPs...")
        xp = cp if GPU else np
        res_names = np.full(len(chrs), 'UNKNOWN', dtype='U20')
        res_dists = np.full(len(chrs), -1, dtype=np.int64)
        unique_chrs = np.unique(chrs)
        for chrom in unique_chrs:
            if chrom not in self.genes_by_chr: continue
            mask = (chrs == chrom)
            snp_pos = xp.asarray(positions[mask])
            genes = self.genes_by_chr[chrom]
            g_starts = xp.asarray([g[0] for g in genes])
            g_ends = xp.asarray([g[1] for g in genes])
            g_names = np.array([g[2] for g in genes])
            
            p = snp_pos[:, None]
            s = g_starts[None, :]
            e = g_ends[None, :]
            d_up = s - p
            d_down = p - e
            dist_mat = xp.where((p >= s) & (p <= e), 0,
                        xp.where(d_up > 0, d_up,
                        xp.where(d_down > 0, d_down, 2e9)))
            
            min_idx = xp.argmin(dist_mat, axis=1)
            min_vals = xp.min(dist_mat, axis=1)
            idx_cpu = cp.asnumpy(min_idx) if GPU else min_idx
            vals_cpu = cp.asnumpy(min_vals) if GPU else min_vals
            mask_idx = np.where(mask)[0]
            res_names[mask_idx] = g_names[idx_cpu]
            res_dists[mask_idx] = vals_cpu.astype(np.int64)
            if GPU: cp.get_default_memory_pool().free_all_blocks()
        return res_names, res_dists

def main():
    male_file = "combined-gwas-results-all-ctrl-males.txt.gz"
    female_file = "combined-gwas-results-all-ctrl-females.txt.gz"
    gtf_file = "Homo_sapiens.GRCh38.105.gtf.gz"
    chr_labels = [str(i) for i in range(1, 23)] + ['X']

    males = load_gwas_numpy(male_file)
    females = load_gwas_numpy(female_file)
    if not males or not females: return

    males = compute_coordinates(males, chr_labels)
    females = compute_coordinates(females, chr_labels)

    THRESHOLD = 5e-8
    gtf = GTFAnnotatorGPU(gtf_file)
    gtf.load()
    
    m_mask = males['pvalue'] < THRESHOLD
    f_mask = females['pvalue'] < THRESHOLD
    
    s_chr = np.concatenate([males['CHR'][m_mask], females['CHR'][f_mask]])
    s_pos = np.concatenate([males['POS'][m_mask], females['POS'][f_mask]])
    s_log = np.concatenate([males['log10p'][m_mask], females['log10p'][f_mask]])
    s_cum = np.concatenate([males['pos_cumulative'][m_mask], females['pos_cumulative'][f_mask]])
    s_set = np.concatenate([np.full(m_mask.sum(), 'Males'), np.full(f_mask.sum(), 'Females')])
    
    names, dists = gtf.annotate(s_chr, s_pos)
    
    np.savez_compressed("miami_intermediate_males.npz", **males)
    np.savez_compressed("miami_intermediate_females.npz", **females)
    np.savez_compressed("miami_intermediate_annotations.npz", 
                        dataset=s_set, pos_cumulative=s_cum, log10p=s_log, 
                        gene_name=names, gene_dist=dists, pvalue=10**(-s_log))
    with open("miami_intermediate_chr_centers.csv", 'w') as f:
        for c, val in males['chr_centers'].items():
            f.write(f"{c},{val}\n")
    print("Fase 1 Completada.")

if __name__ == "__main__":
    main()
