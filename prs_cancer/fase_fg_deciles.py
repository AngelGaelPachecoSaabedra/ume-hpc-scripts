"""
FASE F+G — 6 scores: integración PRS + caso-control, deciles, OR, gráficas
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import statsmodels.formula.api as smf
import warnings
warnings.filterwarnings('ignore')

PRS_PATH = "/mnt/cephfs/orgs/home/angel.pacheco/prs_cancer/prs_cancer_mcps_all_samples.tsv"
CC_PATH  = "/mnt/cephfs/orgs/home/angel.pacheco/Ordanamiento_Bases/MCPS_Cohorte_Casos_Controles_Oncologicos_F145K.csv"
OUT_DIR  = "/mnt/cephfs/orgs/home/angel.pacheco/prs_cancer"

# ================================================================
# FASE F: JOIN
# ================================================================
df_prs = pd.read_csv(PRS_PATH, sep='\t', na_values='NA')
df_cc  = pd.read_csv(CC_PATH, low_memory=False)

Z_COLS = [c for c in df_prs.columns if c.startswith('PRS_Z_')]

df = df_cc.merge(df_prs[['IID'] + Z_COLS],
                 left_on='IID_EXOME', right_on='IID', how='inner')
df.drop(columns='IID', inplace=True)

print(f"Join: {len(df):,} participantes | Casos: {(df['PHENOTYPE_CANCER_STATUS']=='CASO').sum():,} | "
      f"Controles: {(df['PHENOTYPE_CANCER_STATUS']=='CONTROL').sum():,}")

df.to_csv(f"{OUT_DIR}/prs_cancer_integrated.tsv", sep='\t', index=False,
          float_format='%.6f', na_rep='NA')

# ================================================================
# CONFIGURACIÓN DE LOS 6 SCORES
# ================================================================
SCORE_CONFIGS = [
    {'pgs_id':'PGS003766','z_col':'PRS_Z_PGS003766_PROSTATE', 'site':'PROSTATE',
     'sex_val':1,    'sex_label':'Hombres',      'color':'#1565C0',
     'label':'Próstata (PGS003766)\nConti 2021 — 386 SNPs',    'n_snps':386},
    {'pgs_id':'PGS000028','z_col':'PRS_Z_PGS000028_BREAST',   'site':'BREAST',
     'sex_val':0,    'sex_label':'Mujeres',       'color':'#C62828',
     'label':'Mama (PGS000028)\nMavaddat 2019 — 79 SNPs',      'n_snps':79},
    {'pgs_id':'PGS000078','z_col':'PRS_Z_PGS000078_LUNG',     'site':'LUNG',
     'sex_val':None, 'sex_label':'Ambos sexos',  'color':'#2E7D32',
     'label':'Pulmón (PGS000078)\nMcKay 2017 — 83 SNPs',      'n_snps':83},
    {'pgs_id':'PGS000073','z_col':'PRS_Z_PGS000073_CERVICAL', 'site':'CERVICAL',
     'sex_val':0,    'sex_label':'Mujeres',       'color':'#6A1B9A',
     'label':'Cervical (PGS000073)\nChen 2024 — 10 SNPs',     'n_snps':10},
    {'pgs_id':'PGS002265','z_col':'PRS_Z_PGS002265_GI',       'site':'GI',
     'sex_val':None, 'sex_label':'Ambos sexos',  'color':'#E65100',
     'label':'GI/Colorrectal (PGS002265)\nHuyghe 2020 — 125 SNPs', 'n_snps':125},
    {'pgs_id':'PGS000874','z_col':'PRS_Z_PGS000874_HEMATOLOGIC','site':'HEMATOLOGIC',
     'sex_val':None, 'sex_label':'Ambos sexos',  'color':'#00695C',
     'label':'Hematológico/CLL (PGS000874)\nBerndt 2018 — 36 SNPs', 'n_snps':36},
]

# ================================================================
# FUNCIONES DE ANÁLISIS
# ================================================================
def compute_decile_table(df_sub, z_col, site):
    df_w = df_sub[df_sub[z_col].notna()].copy()
    df_w['DECIL'] = pd.qcut(df_w[z_col], q=10, labels=False) + 1

    rows = []
    for d in range(1, 11):
        sub = df_w[df_w['DECIL'] == d]
        n_total      = len(sub)
        n_caso_site  = (sub['CANCER_SITE'] == site).sum()
        pct_site     = n_caso_site / n_total * 100 if n_total > 0 else np.nan
        rows.append({'DECIL': d, 'N_TOTAL': n_total,
                     'N_CASO_SITE': n_caso_site, 'PCT_CASO_SITE': pct_site})

    df_d = pd.DataFrame(rows)

    # OR vs decil 5
    ref = df_d[df_d['DECIL'] == 5].iloc[0]
    ref_c = ref['N_CASO_SITE']
    ref_n = ref['N_TOTAL'] - ref_c
    ors, lo_, hi_ = [], [], []
    for _, row in df_d.iterrows():
        c = row['N_CASO_SITE']; n = row['N_TOTAL'] - c
        if row['DECIL'] == 5:
            ors.append(1.0); lo_.append(np.nan); hi_.append(np.nan); continue
        if ref_c == 0 or ref_n == 0 or n == 0 or c == 0:
            ors.append(np.nan); lo_.append(np.nan); hi_.append(np.nan); continue
        or_v = (c * ref_n) / (n * ref_c)
        se   = np.sqrt(1/c + 1/n + 1/ref_c + 1/ref_n)
        ors.append(round(or_v, 3))
        lo_.append(round(np.exp(np.log(or_v) - 1.96*se), 3))
        hi_.append(round(np.exp(np.log(or_v) + 1.96*se), 3))

    df_d['OR_vs_D5'] = ors; df_d['OR_IC95_LOW'] = lo_; df_d['OR_IC95_HI'] = hi_
    return df_d, df_w

def trend_test(df_w, z_col, site):
    df_t = df_w[[z_col, 'CANCER_SITE']].dropna().copy()
    df_t['Y'] = (df_t['CANCER_SITE'] == site).astype(int)
    df_t = df_t.rename(columns={z_col: 'PRS_Z'})
    if df_t['Y'].sum() < 5: return None
    try:
        mod   = smf.logit('Y ~ PRS_Z', data=df_t).fit(disp=False)
        beta  = mod.params['PRS_Z']
        pval  = mod.pvalues['PRS_Z']
        ci    = mod.conf_int().loc['PRS_Z']
        return {'OR_per_1SD': round(np.exp(beta), 3),
                'CI_lo':      round(np.exp(ci[0]), 3),
                'CI_hi':      round(np.exp(ci[1]), 3),
                'p_value':    pval,
                'n_cases':    int(df_t['Y'].sum()),
                'n_total':    len(df_t)}
    except: return None

# ================================================================
# FASE G: LOOP POR SCORE
# ================================================================
print("\n" + "="*65)
print("FASE G — DECILES, OR Y GRÁFICAS (6 SCORES)")
print("="*65)

all_results  = {}
all_trends   = {}
summary_rows = []

for cfg in SCORE_CONFIGS:
    pgs_id = cfg['pgs_id']; z_col = cfg['z_col']
    site   = cfg['site'];   color = cfg['color']

    if cfg['sex_val'] is not None:
        df_sub = df[df['MALE'] == cfg['sex_val']].copy()
    else:
        df_sub = df.copy()

    n_valid  = df_sub[z_col].notna().sum()
    n_casos  = (df_sub['CANCER_SITE'] == site).sum()

    print(f"\n>> {pgs_id} — {site} ({cfg['sex_label']})")
    print(f"   Subcohorte: {n_valid:,} | Casos {site}: {n_casos:,}")

    df_d, df_w = compute_decile_table(df_sub, z_col, site)
    trend      = trend_test(df_w, z_col, site)

    all_results[pgs_id] = df_d
    all_trends[pgs_id]  = trend

    # Tabla en consola
    sig = ''
    if trend:
        sig = '***' if trend['p_value']<0.001 else ('**' if trend['p_value']<0.01
              else ('*' if trend['p_value']<0.05 else 'ns'))
        print(f"   OR/1SD={trend['OR_per_1SD']} [{trend['CI_lo']}–{trend['CI_hi']}] "
              f"p={trend['p_value']:.2e} {sig}")

    print(f"   {'D':>3} {'N':>8} {'Casos':>6} {'%Casos':>8} {'OR_D5':>7} {'IC95':>18}")
    print(f"   {'-'*55}")
    for _, row in df_d.iterrows():
        ic = (f"[{row['OR_IC95_LOW']:.2f}–{row['OR_IC95_HI']:.2f}]"
              if pd.notna(row['OR_IC95_LOW']) else "    [ref]")
        or_s = f"{row['OR_vs_D5']:.2f}" if pd.notna(row['OR_vs_D5']) else "ref"
        print(f"   {int(row['DECIL']):>3} {int(row['N_TOTAL']):>8} "
              f"{int(row['N_CASO_SITE']):>6} {row['PCT_CASO_SITE']:>7.3f}% "
              f"{or_s:>7} {ic:>18}")

    df_d.to_csv(f"{OUT_DIR}/deciles_{pgs_id}_{site}.tsv",
                sep='\t', index=False, float_format='%.4f', na_rep='NA')

    # Summary row
    d1  = df_d[df_d['DECIL']==1 ]['PCT_CASO_SITE'].values[0]
    d5  = df_d[df_d['DECIL']==5 ]['PCT_CASO_SITE'].values[0]
    d10 = df_d[df_d['DECIL']==10]['PCT_CASO_SITE'].values[0]
    or10 = df_d[df_d['DECIL']==10]['OR_vs_D5'].values[0]
    or10_lo = df_d[df_d['DECIL']==10]['OR_IC95_LOW'].values[0]
    or10_hi = df_d[df_d['DECIL']==10]['OR_IC95_HI'].values[0]
    summary_rows.append({
        'CANCER_SITE': site, 'PGS_ID': pgs_id,
        'N_CASOS': n_casos, 'N_COHORTE': n_valid,
        'PCT_D1': round(d1, 3), 'PCT_D5': round(d5, 3), 'PCT_D10': round(d10, 3),
        'OR_D10_vs_D5': round(or10, 3) if pd.notna(or10) else np.nan,
        'OR_D10_IC_LO': round(or10_lo, 3) if pd.notna(or10_lo) else np.nan,
        'OR_D10_IC_HI': round(or10_hi, 3) if pd.notna(or10_hi) else np.nan,
        'OR_per_1SD': trend['OR_per_1SD'] if trend else np.nan,
        'OR_1SD_CI_lo': trend['CI_lo'] if trend else np.nan,
        'OR_1SD_CI_hi': trend['CI_hi'] if trend else np.nan,
        'p_value': trend['p_value'] if trend else np.nan,
        'significance': sig if trend else 'NA',
    })

    # ---- Barplot %casos por decil ----
    fig, ax = plt.subplots(figsize=(11, 6))
    x = df_d['DECIL'].values; y = df_d['PCT_CASO_SITE'].values
    avg = np.nanmean(y)
    bars = ax.bar(x, y, color=color, alpha=0.82, width=0.7, edgecolor='white', linewidth=0.5)
    bars[9].set_edgecolor('black'); bars[9].set_linewidth(1.5)
    ax.axhline(avg, linestyle='--', color='gray', linewidth=0.9, alpha=0.8)
    ax.text(10.55, avg, f'Media\n{avg:.2f}%', va='center', fontsize=8, color='gray')
    for bar, yv in zip(bars, y):
        if not np.isnan(yv):
            ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+avg*0.03,
                    f'{yv:.2f}%', ha='center', va='bottom', fontsize=7.5, color='#333')
    ax.set_xlabel('Decil de PRS (1 = menor riesgo · 10 = mayor riesgo)', fontsize=11)
    ax.set_ylabel(f'% Casos ({site})', fontsize=11)
    ax.set_title(cfg['label'].replace('\n',' — '), fontsize=12, fontweight='bold', pad=10)
    ax.set_xticks(x); ax.set_xlim(0.3, 11.4)
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y', linestyle=':', linewidth=0.5, alpha=0.6)
    if trend:
        ax.annotate(
            f"OR/1SD = {trend['OR_per_1SD']} [{trend['CI_lo']}–{trend['CI_hi']}]\np = {trend['p_value']:.2e}  {sig}",
            xy=(0.02,0.97), xycoords='axes fraction', fontsize=8.5, va='top',
            bbox=dict(boxstyle='round,pad=0.4', fc='white', ec='lightgray', alpha=0.9))
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/barplot_deciles_{pgs_id}_{site}.png", dpi=200, bbox_inches='tight')
    plt.savefig(f"{OUT_DIR}/barplot_deciles_{pgs_id}_{site}.pdf", bbox_inches='tight')
    plt.close()

    # ---- OR por decil ----
    fig2, ax2 = plt.subplots(figsize=(11, 5))
    df_p = df_d.copy()
    x2 = df_p['DECIL'].values; ors = df_p['OR_vs_D5'].values
    lo = df_p['OR_IC95_LOW'].values; hi = df_p['OR_IC95_HI'].values
    ax2.axhline(1.0, linestyle='--', color='gray', linewidth=1.0, alpha=0.8)
    valid = ~np.isnan(ors) & ~np.isnan(lo) & ~np.isnan(hi)
    ax2.errorbar(x2[valid], ors[valid],
                 yerr=[ors[valid]-lo[valid], hi[valid]-ors[valid]],
                 fmt='o', color=color, markersize=7, capsize=4, linewidth=1.5)
    ax2.plot(5, 1.0, 'D', color='gray', markersize=8, label='Referencia D5')
    ax2.set_xlabel('Decil de PRS', fontsize=11); ax2.set_ylabel('OR vs Decil 5', fontsize=11)
    ax2.set_title(f"OR por Decil — {cfg['label'].split(chr(10))[0]}", fontsize=12, fontweight='bold')
    ax2.set_xticks(range(1,11)); ax2.set_xlim(0.3,10.7)
    ax2.spines[['top','right']].set_visible(False)
    ax2.grid(axis='y', linestyle=':', linewidth=0.5, alpha=0.6); ax2.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/OR_deciles_{pgs_id}_{site}.png", dpi=200, bbox_inches='tight')
    plt.savefig(f"{OUT_DIR}/OR_deciles_{pgs_id}_{site}.pdf", bbox_inches='tight')
    plt.close()

    print(f"   Gráficas guardadas: barplot + OR deciles")

# ================================================================
# FIGURA COMBINADA 2×3
# ================================================================
print("\n>> Generando figura combinada 2×3...")
fig, axes = plt.subplots(2, 3, figsize=(20, 11))
axes = axes.flatten()

for i, cfg in enumerate(SCORE_CONFIGS):
    ax = axes[i]; pgs_id = cfg['pgs_id']
    df_d = all_results[pgs_id]
    x = df_d['DECIL'].values; ors = df_d['OR_vs_D5'].values
    lo = df_d['OR_IC95_LOW'].values; hi = df_d['OR_IC95_HI'].values
    trend = all_trends[pgs_id]

    ax.axhline(1.0, linestyle='--', color='gray', linewidth=0.9, alpha=0.7)
    valid = ~np.isnan(ors) & ~np.isnan(lo) & ~np.isnan(hi)
    ax.errorbar(x[valid], ors[valid],
                yerr=[ors[valid]-lo[valid], hi[valid]-ors[valid]],
                fmt='o-', color=cfg['color'], markersize=5.5, capsize=3, linewidth=1.5)
    ax.plot(5, 1.0, 'D', color='gray', markersize=7, zorder=5)

    # Anotar D1 y D10
    for d_idx in [0, 9]:
        if pd.notna(ors[d_idx]):
            ax.annotate(f'{ors[d_idx]:.2f}',
                        xy=(x[d_idx], ors[d_idx]),
                        xytext=(x[d_idx]+(-0.5 if d_idx==0 else 0.3), ors[d_idx]),
                        fontsize=7.5, color=cfg['color'], fontweight='bold',
                        ha='right' if d_idx==0 else 'left')

    if trend:
        sig = '***' if trend['p_value']<0.001 else ('**' if trend['p_value']<0.01
              else ('*' if trend['p_value']<0.05 else 'ns'))
        ax.set_title(f"{cfg['label'].replace(chr(10),' | ')}\nOR/SD={trend['OR_per_1SD']} p={trend['p_value']:.1e} {sig}",
                     fontsize=8.5, fontweight='bold', pad=6, color=cfg['color'])
    else:
        ax.set_title(cfg['label'].replace(chr(10),' | '), fontsize=8.5, fontweight='bold', pad=6)

    ax.set_xlabel('Decil PRS', fontsize=9); ax.set_ylabel('OR vs D5', fontsize=9)
    ax.set_xticks(range(1,11)); ax.tick_params(labelsize=8)
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y', linestyle=':', linewidth=0.5, alpha=0.5)

fig.suptitle("OR por Decil de PRS — 6 Tipos de Cáncer (Cohorte MCPS ~138K)\nReferencia: Decil 5 | Solo betas PGS Catalog",
             fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/OR_deciles_combined_6cancer.png", dpi=200, bbox_inches='tight')
plt.savefig(f"{OUT_DIR}/OR_deciles_combined_6cancer.pdf", bbox_inches='tight')
plt.close()
print("Figura 2×3 guardada.")

# ================================================================
# TABLA RESUMEN FINAL
# ================================================================
df_sum = pd.DataFrame(summary_rows)
df_sum.to_csv(f"{OUT_DIR}/resumen_PRS_6cancer.tsv", sep='\t', index=False, float_format='%.4f', na_rep='NA')

print("\n" + "="*72)
print("TABLA RESUMEN EJECUTIVO — 6 SCORES")
print("="*72)
print(f"{'Site':<12} {'PGS':<12} {'N_casos':>8} {'%D1':>7} {'%D5':>7} {'%D10':>7} "
      f"{'OR_D10':>7} {'OR/SD':>7} {'p':>10} {'Sig':>4}")
print("-"*72)
for r in summary_rows:
    p_str = f"{r['p_value']:.2e}" if pd.notna(r['p_value']) else 'NA'
    or10  = f"{r['OR_D10_vs_D5']:.2f}" if pd.notna(r['OR_D10_vs_D5']) else 'NA'
    orsd  = f"{r['OR_per_1SD']:.3f}"   if pd.notna(r['OR_per_1SD'])   else 'NA'
    print(f"{r['CANCER_SITE']:<12} {r['PGS_ID']:<12} {r['N_CASOS']:>8} "
          f"{r['PCT_D1']:>7.3f} {r['PCT_D5']:>7.3f} {r['PCT_D10']:>7.3f} "
          f"{or10:>7} {orsd:>7} {p_str:>10} {r['significance']:>4}")

print("\n==== FASE F+G COMPLETADA ====")
