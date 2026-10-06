#!/usr/bin/env python3
"""
Miami Plot — FASE 2: Plotting
=========================================================
Características V9:
  - Estilo limpio y preciso.
  - Conexión natural de etiquetas (sin huecos).
  - Sin bordes dorados distractores.
  - Línea de referencia (threshold) configurable.
  - Dibuja la línea punteada roja en -log10(p) para ambos lados.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import sys

# --- CONFIGURACIÓN ---
GAP_SIZE = 2.0
SIGNIFICANCE_THRESHOLD = 5e-8  # <--- Modifica esto si cambia tu corte (ej. 1e-6)

try:
    from adjustText import adjust_text
    HAS_ADJUST_TEXT = True
except ImportError:
    HAS_ADJUST_TEXT = False
    print("djustText no encontrado. Se usará posicionamiento estático.")

def load_data(prefix="miami_intermediate"):
    print("Cargando datos intermedios...")
    males = dict(np.load(f"{prefix}_males.npz", allow_pickle=True))
    females = dict(np.load(f"{prefix}_females.npz", allow_pickle=True))
    annot = dict(np.load(f"{prefix}_annotations.npz", allow_pickle=True))
    
    centers = {}
    with open(f"{prefix}_chr_centers.csv") as f:
        for l in f:
            k,v = l.strip().split(',')
            centers[k] = float(v)
    return males, females, centers, annot

def add_smart_annotations(ax, annot, dataset_name, color, limit_y, barrier_line):
    mask = annot['dataset'] == dataset_name
    if not np.any(mask): return

    valid = (mask) & (annot['gene_name'] != 'UNKNOWN') & (annot['gene_dist'] < 500000)
    indices = np.where(valid)[0]
    if len(indices) == 0: return

    TOP_N = 18
    
    sig_order = np.argsort(annot['log10p'][indices])[::-1]
    sorted_indices = indices[sig_order]
    
    final_indices = []
    seen_genes = set()
    
    for idx in sorted_indices:
        gene = str(annot['gene_name'][idx])
        if gene not in seen_genes:
            seen_genes.add(gene)
            final_indices.append(idx)
        if len(final_indices) >= TOP_N:
            break
            
    indices = np.array(final_indices)
    
    texts = []
    
    if dataset_name == 'Males':
        direction = 1
        va_setting = 'bottom'
    else:
        direction = -1
        va_setting = 'top'
    
    print(f"Etiquetando {len(indices)} genes TOP para {dataset_name}...")

    # Coordenadas Y reales
    y_points = (annot['log10p'][indices] + GAP_SIZE) * direction

    # 1. RESALTADO SUTIL
    ax.scatter(annot['pos_cumulative'][indices], 
               y_points,
               color=color, 
               s=25,        
               zorder=10,
               edgecolors='none')

    for i, idx in enumerate(indices):
        x_real = float(annot['pos_cumulative'][idx])
        y_real = y_points[i]
        name = str(annot['gene_name'][idx])
        
        # 2. INICIALIZAR TEXTO
        y_bias = y_real + (limit_y * 0.01 * direction)
        
        t = ax.text(x_real, y_bias, name,
                    fontsize=9, fontweight='bold', color=color,
                    ha='center', va=va_setting,
                    zorder=11,
                    bbox=dict(boxstyle='round,pad=0.1', fc='white', alpha=0.6, ec='none'))
        texts.append(t)

    # 3. CONEXIÓN NATURAL
    if HAS_ADJUST_TEXT:
        arrow_props = dict(
            arrowstyle='-',      
            color='#444444',     
            alpha=0.7,
            lw=0.8,              
            shrinkA=0, shrinkB=0 
        )
        
        print(f"      ...optimizando posiciones...")
        
        adjust_text(texts, ax=ax, 
                    only_move={'text': 'y', 'points': 'y'}, 
                    arrowprops=arrow_props,
                    add_objects=[barrier_line],
                    expand_points=(1.4, 1.6), 
                    avoid_self=True,
                    lim=1000)

def main():
    try:
        males, females, centers, annot = load_data()
    except Exception as e:
        print(f"Error cargando datos: {e}")
        return

    chr_labels = [str(i) for i in range(1, 23)] + ['X']
    
    colors_m = ['#1f77b4', '#6baed6']
    colors_f = ['#d62728', '#fb6a4a']

    fig, ax = plt.subplots(figsize=(22, 14))
    
    print("Generando scatter plot...")
    
    for i, chr_name in enumerate(chr_labels):
        c_idx = i % 2
        mask_m = (males['CHR'] == chr_name)
        mask_f = (females['CHR'] == chr_name)
        
        if np.any(mask_m):
            ax.scatter(males['pos_cumulative'][mask_m], 
                       males['log10p'][mask_m] + GAP_SIZE, 
                       s=3, c=colors_m[c_idx], alpha=0.8, rasterized=True, edgecolors='none')
            
        if np.any(mask_f):
            ax.scatter(females['pos_cumulative'][mask_f], 
                       -females['log10p'][mask_f] - GAP_SIZE, 
                       s=3, c=colors_f[c_idx], alpha=0.8, rasterized=True, edgecolors='none')

    
    # --- LÍNEAS DE REFERENCIA (GAP Y SIGNIFICANCIA) ---
    barrier_m = ax.axhline(GAP_SIZE, color='gray', lw=0.8, ls='-', alpha=0.3, zorder=1)
    barrier_f = ax.axhline(-GAP_SIZE, color='gray', lw=0.8, ls='-', alpha=0.3, zorder=1)

    # Calcular altura visual del threshold
    thresh_log = -np.log10(SIGNIFICANCE_THRESHOLD)
    thresh_y_m = thresh_log + GAP_SIZE
    thresh_y_f = -thresh_log - GAP_SIZE

    print(f"Dibujando línea de corte p={SIGNIFICANCE_THRESHOLD} en Y=±{thresh_log:.2f}")

    # Línea Hombres
    ax.axhline(thresh_y_m, color='#CC0000', linestyle='--', linewidth=1.0, alpha=0.6, zorder=2)
    # Línea Mujeres
    ax.axhline(thresh_y_f, color='#CC0000', linestyle='--', linewidth=1.0, alpha=0.6, zorder=2)

    # Etiquetas centrales de cromosomas
    ax.set_xticks([])
    for chr_name in chr_labels:
        x_pos = centers.get(chr_name, 0)
        ax.text(x_pos, 0, chr_name, 
                ha='center', va='center', 
                fontsize=12, fontweight='bold', color='#333333')

    # Configuración Eje Y
    max_val = max(males['log10p'].max(), females['log10p'].max())
    limit_y = max_val + GAP_SIZE + (max_val * 0.2)
    ax.set_ylim(-limit_y, limit_y)
    
    tick_step = 10
    if max_val > 100: tick_step = 20
    
    ticks_raw = np.arange(tick_step, max_val + tick_step, tick_step)
    pos_ticks_loc = ticks_raw + GAP_SIZE
    neg_ticks_loc = -ticks_raw - GAP_SIZE
    
    ax.set_yticks(np.concatenate([neg_ticks_loc, pos_ticks_loc]))
    ax.set_yticklabels(np.concatenate([ticks_raw.astype(str), ticks_raw.astype(str)]), fontsize=11)
    ax.set_ylabel(r"$-\log_{10}(P)$", fontsize=16, fontweight='bold', labelpad=15)
    
    # Anotaciones
    add_smart_annotations(ax, annot, 'Males', '#0c2c52', limit_y, barrier_m)
    add_smart_annotations(ax, annot, 'Females', '#5e0d0d', limit_y, barrier_f)
    
    # Estilo final
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_position(('outward', 10))
    
    ax.text(0.01, 0.97, 'MALES', transform=ax.transAxes, fontweight='bold', color=colors_m[0], fontsize=20)
    ax.text(0.01, 0.03, 'FEMALES', transform=ax.transAxes, fontweight='bold', color=colors_f[0], fontsize=20)
    
    plt.title("Miami Plot - GWAS Results", fontsize=22, pad=30, fontweight='bold')
    plt.tight_layout()
    output_file = "miami_plot_final_v9_threshold.png"
    plt.savefig(output_file, dpi=300)
    print(f"Plot generado: {output_file}")

if __name__ == "__main__":
    main()
