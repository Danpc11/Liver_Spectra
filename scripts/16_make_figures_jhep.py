"""16 — Figures 1–6 for the Journal of Hepatology version (clinical framing). Requires scripts 01–15 outputs.
Also writes Table 1 (top lineage-intrinsic targets with drug landscape) as CSV for the Word draft.
"""
import os, numpy as np, pandas as pd, statsmodels.api as sm, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, scienceplots, seaborn as sns
from sklearn.mixture import GaussianMixture
from adjustText import adjust_text
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from common import *

plt.style.use(['science', 'nature', 'no-latex'])
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial', 'Helvetica'], 'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6,
                     'xtick.top': False, 'ytick.right': False, 'xtick.minor.visible': False, 'ytick.minor.visible': False, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False, 'pdf.fonttype': 42})
MM = 1 / 25.4; W2 = 183 * MM; FO = os.path.join(FIG, 'jhep'); os.makedirs(FO, exist_ok=True)
ST = ['Normal', 'F0', 'F1', 'F2', 'F3', 'F4']; pal = dict(zip(ST, sns.color_palette('viridis', 6)))
OI = {'blue': '#0072B2', 'orange': '#E69F00', 'green': '#009E73', 'red': '#D55E00', 'purple': '#CC79A7', 'sky': '#56B4E9', 'grey': '#7F7F7F', 'yellow': '#F0E442'}
def save(fig, name): fig.savefig(f'{FO}/{name}.pdf', bbox_inches='tight'); fig.savefig(f'{FO}/{name}.png', dpi=400, bbox_inches='tight'); fig.savefig(f'{FO}/{name}.tiff', dpi=300, bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'}); plt.close(fig)
def lab(ax, s, dx=-0.2, dy=1.14): ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')

C = pd.read_csv(tab('Table_S6b_composition_scores_per_sample.csv'), index_col=0); C = C[C.estadio != 'Control']
B = pd.read_csv(tab('Table_S7a_threshold_vs_linear.csv')); S = pd.read_csv(tab('Table_S15b_mechanics_scores_per_sample.csv'), index_col=0)
EN = {'Hepatocito': 'Hepatocyte', 'HSC': 'Stellate / myofibroblast', 'Colangiocito': 'Cholangiocyte', 'Macrofago': 'Macrophage', 'Linfocito': 'Lymphocyte', 'Endotelio': 'Endothelium'}
CT_COL = {'Hepatocito': OI['green'], 'HSC': OI['red'], 'Colangiocito': OI['purple'], 'Macrofago': OI['orange'], 'Linfocito': OI['blue'], 'Endotelio': OI['grey']}

# ================= Fig 1
fig, axs = plt.subplots(2, 2, figsize=(W2 * 0.85, W2 * 0.7)); axs = axs.ravel(); plt.subplots_adjust(hspace=0.55, wspace=0.45)
ax = axs[0]
for ct in CT_COL: g = C.groupby('estadio')[ct].agg(['mean', 'sem']).reindex(ST); ax.errorbar(range(6), g['mean'], yerr=g['sem'], fmt='o-', color=CT_COL[ct], ms=2.5, lw=0.9, capsize=1.5, elinewidth=0.6, label=EN[ct])
ax.set_xticks(range(6)); ax.set_xticklabels(ST); ax.axhline(0, color='0.7', lw=0.5); ax.set_ylabel('Composition score (marker z)'); ax.set_xlabel('Fibrosis stage'); ax.set_ylim(-0.9, 1.6); ax.legend(ncol=2, handlelength=1.2, fontsize=5.5, loc='upper left'); lab(ax, 'A')
ax = axs[1]
for st in ST: g = C[C.estadio == st]; ax.scatter(g.Hepatocito, g.offset, s=6, color=pal[st], alpha=.85, lw=0, label=st)
ax.set_xlabel('Hepatocyte score'); ax.set_ylabel('Global transcriptome variance\n(aperiodic offset)'); ax.set_xlim(-2.3, 2.1); ax.set_ylim(3.09, 3.41); ax.text(0.03, 0.98, 'Stage effect: t = −8.0 → +1.4\nafter adjusting for composition', transform=ax.transAxes, fontsize=6, va='top'); ax.legend(loc='lower right', title='Stage', title_fontsize=6, markerscale=1.5); lab(ax, 'B')
ax = axs[2]; H = B.set_index('variable')[[c for c in B.columns if c.startswith('dAIC_step_at_')]]; H.columns = ['≥F1', '≥F2', '≥F3', '≥F4']; H = H.loc[['Hepatocito', 'HSC', 'Colangiocito', 'Macrofago', 'Linfocito', 'Endotelio']]; H.index = [EN[i] for i in H.index]
im = ax.imshow(H.values, cmap='RdBu_r', vmin=-15, vmax=15, aspect='auto'); ax.set_xticks(range(4)); ax.set_xticklabels(H.columns); ax.set_yticks(range(len(H))); ax.set_yticklabels(H.index); ax.set_xlabel('Step from stage')
for i in range(H.shape[0]):
    for j in range(H.shape[1]): ax.text(j, i, f'{H.values[i, j]:.0f}', ha='center', va='center', fontsize=6, color='white' if abs(H.values[i, j]) > 9 else 'k')
cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03); cb.set_label('ΔAIC (step − linear)', fontsize=6); cb.ax.tick_params(labelsize=5.5); ax.spines[['left', 'bottom']].set_visible(False); ax.tick_params(length=0); lab(ax, 'C', -0.5)
ax = axs[3]; y3 = C[C.estadio == 'F3'].Hepatocito.values; y4 = C[C.estadio == 'F4'].Hepatocito.values; yN = C[C.estadio == 'Normal'].Hepatocito.values
bins = np.linspace(-2.2, 1.2, 28); ax.hist(y3, bins=bins, color=pal['F3'], alpha=0.8, lw=0, label=f'F3 (n={len(y3)})', density=True); ax.hist(y4, bins=bins, color=pal['F4'], alpha=0.45, lw=0, label=f'F4 (n={len(y4)})', density=True)
g2 = GaussianMixture(2, random_state=0, n_init=5).fit(y3.reshape(-1, 1)); xs = np.linspace(-2.2, 1.2, 200)
for k in range(2): ax.plot(xs, g2.weights_[k] * (1 / np.sqrt(2 * np.pi * g2.covariances_[k, 0, 0])) * np.exp(-(xs - g2.means_[k, 0]) ** 2 / (2 * g2.covariances_[k, 0, 0])), color='k', lw=0.8, ls='--' if g2.weights_[k] < 0.5 else '-')
ax.set_xlabel('Hepatocyte score'); ax.set_ylabel('Density'); ax.set_ylim(0, 2.1); ax.legend(loc='upper left'); ax.text(0.98, 0.95, f'F3 minor component: {100*g2.weights_.min():.0f}%\nΔBIC(1−2) = 12.8', transform=ax.transAxes, ha='right', va='top', fontsize=6); lab(ax, 'D')
save(fig, 'Fig1_linear_replacement_thresholds')

# ================= Fig 2 (mechanics)
PROG = {'matrix_crosslinking': ('Matrix crosslinking', OI['grey']), 'YAP_TAZ_targets': ('YAP/TAZ targets', OI['red']), 'mechanosensing': ('Mechanosensing machinery', OI['orange']), 'TGFb_SMAD': ('TGF-β / SMAD', OI['blue']), 'senescence_SASP': ('Senescence / SASP', OI['purple'])}
R = pd.read_csv(tab('Table_S15_mechanics_scores_thresholds.csv')).set_index('score')
fig, axs = plt.subplots(1, 3, figsize=(W2, W2 * 0.34)); plt.subplots_adjust(wspace=0.5)
ax = axs[0]
for k, (nm, col) in PROG.items(): g = S.groupby('estadio')[k].agg(['mean', 'sem']).reindex(ST); ax.errorbar(range(6), g['mean'], yerr=g['sem'], fmt='o-', color=col, ms=2.5, lw=0.9, capsize=1.5, elinewidth=0.6, label=nm)
ax.set_xticks(range(6)); ax.set_xticklabels(ST); ax.axhline(0, color='0.7', lw=0.5); ax.set_ylabel('Programme score (marker z)'); ax.set_xlabel('Fibrosis stage'); ax.legend(fontsize=5.5, loc='upper left'); lab(ax, 'A')
ax = axs[1]; x = np.arange(len(PROG)); ax.bar(x - 0.2, [R.loc[k, 't_step_F4'] for k in PROG], 0.4, color='0.6', lw=0, label='unadjusted'); ax.bar(x + 0.2, [R.loc[k, 't_step_F4_comp_adj'] for k in PROG], 0.4, color=OI['orange'], lw=0, label='adjusted for cell composition')
ax.axhline(2, color='0.4', ls='--', lw=0.6); ax.axhline(0, color='k', lw=0.5); ax.set_xticks(x); ax.set_xticklabels(['Matrix\ncrosslinking', 'YAP/TAZ\ntargets', 'Mechano-\nsensing', 'TGF-β/\nSMAD', 'Senescence/\nSASP'], fontsize=6); ax.set_ylabel('t of step at F4'); ax.set_ylim(0, 7.2); ax.legend(fontsize=6, loc='upper right'); lab(ax, 'B')
ax = axs[2]
for st in ST: g = S[S.estadio == st]; ax.scatter(g.HSC, g.YAP_TAZ_targets, s=7, color=pal[st], alpha=.85, lw=0, label=st)
sub = S[S.estadio != 'F4']; b = np.polyfit(sub.HSC, sub.YAP_TAZ_targets, 1); xs = np.linspace(S.HSC.min(), S.HSC.max(), 50); ax.plot(xs, np.polyval(b, xs), color='k', lw=0.8, ls='--', label='fit, Normal–F3')
ax.set_xlabel('Stellate / myofibroblast score'); ax.set_ylabel('YAP/TAZ target score'); ax.legend(fontsize=5.5, loc='upper left', ncol=2, markerscale=1.4); ax.text(0.98, 0.04, 'F4 lies above the\nNormal–F3 relation', transform=ax.transAxes, ha='right', fontsize=6); lab(ax, 'C')
save(fig, 'Fig2_mechanotransduction_switch')

# ================= Fig 3 (neighbourhood coupling)
A = pd.read_csv(tab('Table_S5d_acf_with_without_composition.csv')); D = pd.read_csv(tab('Table_S5g_concordance_by_distance.csv'), index_col=0); SC = pd.read_csv(tab('Table_S8c_within_type_spatial_autocorrelation.csv'))
en = {'Fagocito_mononuclear': 'Mono. phagocyte', 'Endotelio': 'Endothelium', 'Mesenquima_HSC': 'Mesenchyme/HSC', 'Colangiocito': 'Cholangiocyte', 'Plasma': 'Plasma', 'T_NK': 'T/NK', 'B': 'B', 'pDC': 'pDC'}
fig, axs = plt.subplots(1, 3, figsize=(W2, W2 * 0.34)); plt.subplots_adjust(wspace=0.55)
ax = axs[0]; ax.axhspan(-A['null_p97.5'].max(), A['null_p97.5'].max(), color='0.9', lw=0, label='95% permutation null'); ax.plot(A.lag, A.acf_unadjusted, 'o-', color='k', ms=3, lw=0.9, label='Stage effect per gene'); ax.plot(A.lag, A.acf_composition_adjusted, 's--', color=OI['orange'], ms=3, lw=0.9, label='Adjusted for cell composition'); ax.set_xscale('log'); ax.set_xlabel('Distance (genes)'); ax.set_ylabel('Spatial autocorrelation'); ax.legend(); lab(ax, 'A')
ax = axs[1]; xx = range(len(D)); ax.plot(xx, D.r_t, 'o-', color='k', ms=3, lw=0.9, label='Concordance of stage effect'); ax.plot(xx, D.r_t_comp_adj, 'o--', color=OI['orange'], ms=3, lw=0.9, label='…adjusted for composition'); ax.plot(xx, D.coexpr, 's-', color=OI['blue'], ms=3, lw=0.9, label='Co-expression across biopsies')
ax.set_xticks(list(xx)); ax.set_xticklabels(['overlap', '<1', '1–10', '10–50', '50–100', '100–500', '>500'], rotation=45, ha='right'); ax.set_xlabel('Intergenic distance (kb)'); ax.set_ylabel('Correlation'); ax.axhline(0, color='0.6', lw=0.5); ax.legend(loc='upper right'); ax.set_ylim(-0.08, 0.5); lab(ax, 'B')
ax = axs[2]; SC = SC.sort_values('acf_lag1'); ax.barh(range(len(SC)), SC.acf_lag1, color=OI['red'], lw=0, height=0.65); ax.errorbar([0] * len(SC), range(len(SC)), xerr=SC['null_p97.5'], fmt='none', color='k', capsize=2, lw=0.6, label='95% null'); ax.set_yticks(range(len(SC))); ax.set_yticklabels([en[t] for t in SC.type]); ax.set_xlabel('Within-cell-type autocorrelation\n(cirrhosis vs healthy, scRNA-seq)'); ax.legend(loc='lower right'); lab(ax, 'C', -0.6)
save(fig, 'Fig3_neighbourhood_coupling')

# ================= Fig 4 (two layers: GSEA + TF)
G1 = pd.read_csv(tab('Table_S10a_GSEA_hallmark_ordinal.csv')).set_index('Term'); G2 = pd.read_csv(tab('Table_S10b_GSEA_hallmark_composition_adjusted.csv')).set_index('Term'); TF = pd.read_csv(tab('Table_S10d_DoRothEA_TF_activity_vs_stage.csv')).set_index('tf')
paths = ['TNFA_SIGNALING_VIA_NFKB', 'EPITHELIAL_MESENCHYMAL_TRANSITION', 'ALLOGRAFT_REJECTION', 'INFLAMMATORY_RESPONSE', 'ANGIOGENESIS', 'P53_PATHWAY', 'IL6_JAK_STAT3_SIGNALING', 'TGF_BETA_SIGNALING', 'ADIPOGENESIS', 'FATTY_ACID_METABOLISM', 'OXIDATIVE_PHOSPHORYLATION', 'PEROXISOME', 'XENOBIOTIC_METABOLISM', 'BILE_ACID_METABOLISM']
paths = [p for p in paths if p in G1.index and p in G2.index]; pretty = lambda p: p.replace('_', ' ').title().replace('Tnfa Signaling Via Nfkb', 'TNFα / NF-κB').replace('Il6 Jak Stat3 Signaling', 'IL6–JAK–STAT3').replace('Tgf Beta', 'TGF-β').replace('P53', 'p53')
fig, axs = plt.subplots(1, 2, figsize=(W2, W2 * 0.4), gridspec_kw={'width_ratios': [1.15, 1]}); plt.subplots_adjust(wspace=0.65)
ax = axs[0]; y = np.arange(len(paths))
for i, p in enumerate(paths): ax.plot([G1.loc[p, 'NES'], G2.loc[p, 'NES']], [i, i], color='0.7', lw=1)
ax.scatter([G1.loc[p, 'NES'] for p in paths], y, s=18, color='0.35', zorder=3, label='Raw stage statistic'); ax.scatter([G2.loc[p, 'NES'] for p in paths], y, s=18, color=OI['orange'], zorder=3, label='Composition-adjusted')
ax.axvline(0, color='k', lw=0.5); ax.set_yticks(y); ax.set_yticklabels([pretty(p) for p in paths]); ax.set_xlabel('GSEA normalised enrichment score'); ax.invert_yaxis(); ax.legend(loc='lower right'); lab(ax, 'A', -0.75)
ax = axs[1]; tfs = ['HNF4A', 'HNF1A', 'FOXA1', 'FOXA2', 'RXRA', 'PPARA', 'NFKB1', 'RELA', 'TP53', 'SMAD3', 'SMAD4', 'KLF4', 'FOXM1', 'JUN']; tfs = [t for t in tfs if t in TF.index]; y = np.arange(len(tfs))
for i, t in enumerate(tfs): ax.plot([TF.loc[t, 't_stage'], TF.loc[t, 't_stage_comp_adjusted']], [i, i], color='0.7', lw=1)
ax.scatter([TF.loc[t, 't_stage'] for t in tfs], y, s=18, color='0.35', zorder=3); ax.scatter([TF.loc[t, 't_stage_comp_adjusted'] for t in tfs], y, s=18, color=OI['orange'], zorder=3)
ax.axvline(0, color='k', lw=0.5); ax.set_yticks(y); ax.set_yticklabels(tfs, style='italic'); ax.set_xlabel('t of stage on inferred TF activity'); ax.invert_yaxis(); ax.axhline(4.5, color='0.8', lw=0.5, ls=':'); ax.text(0.98, 0.9, 'Hepatocyte identity', transform=ax.transAxes, ha='right', fontsize=6, color='0.4'); ax.text(0.02, 0.3, 'Intrinsic\nresponse', transform=ax.transAxes, fontsize=6, color='0.4'); lab(ax, 'B', -0.35)
save(fig, 'Fig4_two_layers_pathways_TFs')

# ================= Fig 5 (targets) + Table 1
Cn = pd.read_csv(tab('Table_S9a_candidate_targets_lineage_intrinsic.csv')); Cn['development_stage'] = Cn.development_stage.fillna('not curated')
colmap = {'clinical, liver': OI['red'], 'clinical, inflammasome': OI['orange'], 'clinical, other fibrosis': OI['orange'], 'clinical, other': OI['orange'], 'approved, other': OI['orange'], 'preclinical': OI['sky'], 'tool compounds': OI['sky'], 'biomarker': OI['green'], 'none': '0.7', 'not curated': '0.85'}
fig, ax = plt.subplots(figsize=(W2 * 0.62, W2 * 0.5)); top = Cn.head(40)
ax.scatter(top.t_bulk_comp_adjusted, top.t_sc_max, s=24, c=[colmap.get(s, '0.8') for s in top.development_stage], lw=0.3, edgecolor='k', zorder=3)
texts = [ax.text(r.t_bulk_comp_adjusted, r.t_sc_max, r.gene_name, fontsize=5.5, style='italic') for _, r in top.iterrows()]; ax.set_xlim(0.8, 8.5); ax.set_ylim(1.7, 5.6)
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='0.5', lw=0.3), expand=(1.4, 1.6), force_text=(0.5, 0.8))
hd = [Line2D([], [], marker='o', ls='', color=c, markeredgecolor='k', markeredgewidth=0.3, ms=4, label=l) for l, c in [('Clinical, liver', OI['red']), ('Clinical, other indication', OI['orange']), ('Preclinical / tool compound', OI['sky']), ('Biomarker only', OI['green']), ('No modulator', '0.7'), ('Not curated', '0.85')]]
ax.legend(handles=hd, loc='upper right', ncol=2); ax.set_xlabel('Bulk stage effect, composition-adjusted (t)'); ax.set_ylabel('Maximal within-cell-type effect (t, scRNA-seq)'); save(fig, 'Fig5_lineage_intrinsic_targets')
T1 = Cn.head(25)[['gene_name', 'cell_type_up_sc', 't_bulk_stage', 't_bulk_comp_adjusted', 't_sc_max', 'coupled_neighbour', 'target_class', 'drug_landscape', 'development_stage']].copy()
T1.columns = ['Gene', 'Cell type (scRNA-seq)', 'Bulk t (stage)', 'Bulk t (composition-adjusted)', 'Max within-cell-type t', 'Cis-coupled neighbour', 'Target class', 'Agents / status', 'Development stage']
T1['Cell type (scRNA-seq)'] = T1['Cell type (scRNA-seq)'].str.replace('Mesenquima_HSC', 'stellate').str.replace('Fagocito_mononuclear', 'macrophage').str.replace('Endotelio', 'endothelium').str.replace('Colangiocito', 'cholangiocyte')
T1.round(1).to_csv(os.path.join(FO, 'Table1_targets.csv'), index=False)

# ================= Fig 6 (GWAS)
Gw = pd.read_csv(tab('Table_S16c_GWAS_catalog_enrichment_by_trait.csv')); order = ['liver fibrosis/cirrhosis', 'MASLD/NAFLD', 'ALT (liver enzyme)', 'height (control)', 'educational attainment (control)']; Gw = Gw.set_index('trait_set').loc[order]
lbl = ['Fibrosis/\ncirrhosis', 'MASLD', 'ALT', 'Height', 'Education']; cols = [OI['red'], OI['orange'], OI['sky'], '0.6', '0.6']
fig, axs = plt.subplots(1, 2, figsize=(W2 * 0.85, W2 * 0.36)); plt.subplots_adjust(wspace=0.45)
ax = axs[0]; ax.bar(range(5), 100 * Gw.frac_in_coupled_pair, color=cols, lw=0, width=0.65); ax.axhline(100 * Gw.bg_coupled.iloc[0], color='k', ls='--', lw=0.7, label='All genes')
for i, (f, p) in enumerate(zip(Gw.frac_in_coupled_pair, Gw.P_coupled)): ax.text(i, 100 * f + 0.8, f'P={p:.3f}' if p < 0.1 else 'n.s.', ha='center', fontsize=5.5)
ax.set_xticks(range(5)); ax.set_xticklabels(lbl, fontsize=6.5); ax.set_xlabel('GWAS trait (grey = non-liver controls)'); ax.set_ylabel('% of locus genes in a\ncis-coupled neighbour pair'); ax.set_ylim(0, 40); ax.legend(loc='upper right'); lab(ax, 'A')
ax = axs[1]; ax.bar(range(5), Gw.acf_windows, color=cols, lw=0, width=0.65); ax.errorbar(range(5), Gw.acf_random, yerr=2 * Gw.acf_random_sd, fmt='o', color='k', ms=3, capsize=3, lw=0.6, label='Random windows ± 2 s.d.')
for i, p in enumerate(Gw.P_acf): ax.text(i, Gw.acf_windows.iloc[i] + 0.008, f'P={p:.3f}' if p < 0.1 else 'n.s.', ha='center', fontsize=5.5)
ax.set_xticks(range(5)); ax.set_xticklabels(lbl, fontsize=6.5); ax.set_xlabel('GWAS trait'); ax.set_ylabel('Spatial autocorrelation of stage\neffect, ±5 genes around loci'); ax.legend(loc='upper right'); lab(ax, 'B')
save(fig, 'Fig6_GWAS_loci_neighbourhoods'); print('figures ok')
