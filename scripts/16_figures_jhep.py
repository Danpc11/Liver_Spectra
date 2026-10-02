"""16 — Journal of Hepatology figure set (final numbering). Run after scripts 01–15.
Fig 1  Study design and the positional framework (methodology, invariance, GTEx tissue specificity)
Fig 2  Linear cell replacement with thresholds at F3–F4
Fig 3  Mechanotransduction switch at F4
Fig 4  Cis-coupled gene neighbourhoods within each cell type
Fig 5  MASLD and cirrhosis GWAS loci in coupled neighbourhoods
Fig 6  Integrated molecular map (circos + GSEA + TF activity + lineage-intrinsic targets)
Table 1 Top lineage-intrinsic targets with drug landscape.
Outputs: results/figures/jhep/{pdf,png,tiff} and Table1_targets.csv
"""
import os, numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, scienceplots, seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib import image as mpimg
from adjustText import adjust_text
from pycirclize import Circos
from sklearn.mixture import GaussianMixture
from common import *

plt.style.use(['science', 'nature', 'no-latex'])
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial', 'Helvetica'], 'font.size': 7, 'axes.labelsize': 7,
                     'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6, 'xtick.top': False, 'ytick.right': False,
                     'xtick.minor.visible': False, 'ytick.minor.visible': False, 'axes.spines.top': False, 'axes.spines.right': False,
                     'legend.frameon': False, 'pdf.fonttype': 42})
MM = 1 / 25.4; W2 = 183 * MM; FO = os.path.join(FIG, 'jhep'); os.makedirs(FO, exist_ok=True)
ST = ['Normal', 'F0', 'F1', 'F2', 'F3', 'F4']; pal = dict(zip(ST, sns.color_palette('viridis', 6)))
OI = {'blue': '#0072B2', 'orange': '#E69F00', 'green': '#009E73', 'red': '#D55E00', 'purple': '#CC79A7', 'sky': '#56B4E9', 'grey': '#7F7F7F', 'yellow': '#F0E442'}
def save(fig, name):
    fig.savefig(f'{FO}/{name}.pdf', bbox_inches='tight'); fig.savefig(f'{FO}/{name}.png', dpi=400, bbox_inches='tight')
    fig.savefig(f'{FO}/{name}.tiff', dpi=300, bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'}); plt.close(fig)
def lab(ax, s, dx=-0.2, dy=1.15): ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')
bandlab = lambda cs: [str(c).replace('(', '').replace(']', '').replace(', ', '–').replace('.0', '') for c in cs]

grid = load_grid(); N = chrom_lengths(grid); M = pd.read_pickle(inter('meta.pkl')); Wz, Wl = pd.read_pickle(inter('W_adj.pkl'))

C = pd.read_csv(tab('Table_S6b_composition_scores_per_sample.csv'), index_col=0); C = C[C.estadio != 'Control']
B = pd.read_csv(tab('Table_S7a_threshold_vs_linear.csv')); S = pd.read_csv(tab('Table_S15b_mechanics_scores_per_sample.csv'), index_col=0)
EN = {'Hepatocito': 'Hepatocyte', 'HSC': 'Stellate / myofibroblast', 'Colangiocito': 'Cholangiocyte', 'Macrofago': 'Macrophage', 'Linfocito': 'Lymphocyte', 'Endotelio': 'Endothelium'}
CT_COL = {'Hepatocito': OI['green'], 'HSC': OI['red'], 'Colangiocito': OI['purple'], 'Macrofago': OI['orange'], 'Linfocito': OI['blue'], 'Endotelio': OI['grey']}

# ================= Fig 1
# ======================================================= FIGURE 1
U = pd.read_csv(tab('Table_S2a_consensus_spectrum_universal_peaks.csv')); U['chr'] = U.chr.astype(str)
S14 = pd.read_csv(tab('Table_S14c_gtex_universal_peaks_by_tissue.csv')); S14 = S14[S14.gene_set == 'own']
S14e = pd.read_csv(tab('Table_S14e_biopsy_peak_enrichment_by_tissue.csv'))
S11 = pd.read_csv(tab('Table_S11_spectral_fingerprint_real_vs_shuffled.csv'))
pretty = {t: t.replace('_', ' ').replace('brain cerebellar hemisphere', 'brain (cerebellum)').replace('adipose visceral omentum', 'adipose (visceral)')
          .replace('adipose subcutaneous', 'adipose (subcut.)').replace('colon transverse', 'colon').replace('kidney cortex', 'kidney')
          .replace('muscle skeletal', 'muscle').replace('artery aorta', 'aorta').replace('whole blood', 'blood') for t in S14.tissue}

fig = plt.figure(figsize=(W2, W2 * 0.86)); gs = fig.add_gridspec(3, 6, height_ratios=[1.15, 1, 1], hspace=0.85, wspace=1.5)
# a — genome-wide consensus spectrum
ax = fig.add_subplot(gs[0, :])
m = U.set_index(['chr', 'k']).mean_log_w; frac = U.set_index(['chr', 'k']).frac_samples_w_gt3
idx = [(c, k) for c, k in zip(U.chr, U.k)]; x = np.arange(len(U)); chr_of = U.chr.values; u = U.universal.values
for i, c in enumerate(CHR):
    sel = chr_of == c; ax.plot(x[sel], np.exp(U.mean_log_w.values[sel] + EULER), color=OI['grey'] if i % 2 else '#4d4d4d', lw=0.35)
ax.axhline(1, color='k', lw=0.5, ls=':'); ax.scatter(x[u], np.exp(U.mean_log_w.values[u] + EULER), s=3.5, color=OI['red'], lw=0, zorder=5,
                                                     label=f'Universal peaks: >3× background in ≥90% of 437 biopsies (n={u.sum()})')
for c, k, nm, dy in [('7', 131, 'chr7 k=131\nperiod ≈ 7 genes', 12), ('1', 117, 'chr1 k=117\n≈18 genes', 12), ('19', 565, 'chr19 k=565\n≈2.6', 12), ('12', 201, 'chr12 k=201\n≈5', 12), ('2', 168, 'chr2 k=168\n≈7.5', 12)]:
    i = idx.index((c, k)); ax.annotate(nm, (x[i], np.exp(U.mean_log_w.values[i] + EULER)), fontsize=5.2, xytext=(0, dy), textcoords='offset points',
                                       ha='center', arrowprops=dict(arrowstyle='-', lw=0.4, color='0.4'))
xt = [np.where(chr_of == c)[0].mean() for c in CHR]
ax.set_xticks(xt); ax.set_xticklabels([c if (c == 'X' or int(c) % 2 == 1) else '' for c in CHR]); ax.set_yscale('log'); ax.set_ylim(0.3, 110)
ax.set_ylabel('Power / 1/f background'); ax.set_xlabel('Chromosome (spatial frequency increases within each; odd chromosomes labelled)', labelpad=2)
ax.legend(loc='lower left', bbox_to_anchor=(0, 1.04)); lab(ax, 'A', -0.05, 1.38)
# b — period density
ax = fig.add_subplot(gs[1, 0:2]); per = U.period.values
ax.hexbin(np.log10(per), np.log10(np.exp(U.mean_log_w.values + EULER)), gridsize=38, cmap='Greys', mincnt=1, linewidths=0)
ax.scatter(np.log10(per[u]), np.log10(np.exp(U.mean_log_w.values[u] + EULER)), s=3.5, color=OI['red'], lw=0, zorder=5)
ax.set_xticks([0.3, 1, 2, 3]); ax.set_xticklabels(['2', '10', '100', '1000']); ax.set_yticks([0, 0.5, 1]); ax.set_yticklabels(['1', '3', '10'])
ax.set_xlabel('Period (genes)'); ax.set_ylabel('Mean power / background'); lab(ax, 'B', -0.3)
# c — reproducibility
ax = fig.add_subplot(gs[1, 2:4]); ax.hist(U.frac_samples_w_gt3, bins=50, color=OI['grey'], lw=0); ax.set_yscale('log')
ax.axvline(0.9, color=OI['red'], ls='--', lw=0.8); ax.text(0.88, ax.get_ylim()[1] * 0.5, f'n = {u.sum()}', ha='right', fontsize=6, color=OI['red'])
ax.set_xlabel('Fraction of biopsies with power >3×'); ax.set_ylabel('Number of frequencies'); lab(ax, 'C', -0.3)
# d — chr7 peak by stage
ax = fig.add_subplot(gs[1, 4:6]); d = Wl.loc[('7', 131)]
bp = ax.boxplot([d[M.index[M.estadio == s]].values for s in ST], tick_labels=ST, patch_artist=True, showfliers=False, widths=0.62,
                medianprops=dict(color='white', lw=1.2), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6), boxprops=dict(lw=0.5))
for p, s in zip(bp['boxes'], ST): p.set_facecolor(pal[s]); p.set_edgecolor('none')
ax.axhline(4.6, color='0.4', ls=':', lw=0.6); ax.set_ylim(3.3, 10.6); ax.text(0.98, 0.02, 'single-sample P < 0.01', transform=ax.transAxes, ha='right', va='bottom', fontsize=5.3, color='0.45')
ax.set_ylabel('Power / background\n(chr7 k=131)'); ax.set_xlabel('Fibrosis stage'); ax.tick_params(axis='x', labelsize=5.8); lab(ax, 'D', -0.32)
# e — GTEx tissue specificity
ax = fig.add_subplot(gs[2, 0:3]); e = S14e.sort_values('odds_vs_background')
ax.barh(range(len(e)), np.log2(e.odds_vs_background), color=[OI['red'] if t == 'liver' else '0.6' for t in e.tissue], lw=0, height=0.7)
ax.axvline(0, color='k', lw=0.5); ax.set_yticks(range(len(e))); ax.set_yticklabels([pretty[t] for t in e.tissue], fontsize=6)
ax.set_xlabel('Enrichment of the 186 liver-biopsy peaks among\neach GTEx tissue’s own universal peaks (log2 odds)'); lab(ax, 'E', -0.42)
# f — shuffle control
ax = fig.add_subplot(gs[2, 3:6]); tasks = [('early_vs_advanced', 'F0–F1 vs\nF3–F4 (AUC)'), ('normal_vs_advanced', 'Normal vs\nF3–F4 (AUC)'), ('ordinal', 'Stage\n(Spearman ρ)')]
feats = [('spectrum', 'Positional spectrum', OI['blue']), ('spectrum_shuffled_order', 'Spectrum, gene order shuffled', OI['sky']), ('expression_2000', 'Expression (2,000 genes)', '0.6')]
w = 0.26; xx = np.arange(len(tasks))
for i, (f, nm, col) in enumerate(feats):
    vals = [S11[(S11.features == f) & (S11.task == t)].value.iloc[0] for t, _ in tasks]
    ax.bar(xx + (i - 1) * w, vals, w, color=col, lw=0, label=nm)
ax.set_xticks(xx); ax.set_xticklabels([n for _, n in tasks], fontsize=6); ax.set_ylim(0, 1.45); ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax.axhline(0.5, color='0.5', ls=':', lw=0.6); ax.set_ylabel('Leave-one-cohort-out performance'); ax.legend(loc='upper center', ncol=1, bbox_to_anchor=(0.5, 1.03), handlelength=1.1)
lab(ax, 'F', -0.22)
save(fig, 'Fig1_design_and_positional_framework')
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
save(fig, 'Fig2_linear_replacement_thresholds')
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
save(fig, 'Fig3_mechanical_switch')
# ======================================================= FIGURE 4
A = pd.read_csv(tab('Table_S5d_acf_with_without_composition.csv')); D = pd.read_csv(tab('Table_S5g_concordance_by_distance.csv'), index_col=0)
SC = pd.read_csv(tab('Table_S8c_within_type_spatial_autocorrelation.csv')); CO = pd.read_csv(tab('Table_S8b_within_type_neighbour_coexpression.csv'))
EQ = pd.read_csv(tab('Table_S13e_shared_GTEx_eQTL_coupled_vs_uncoupled.csv'))
en = {'Fagocito_mononuclear': 'Mono. phagocyte', 'Endotelio': 'Endothelium', 'Mesenquima_HSC': 'Mesenchyme/HSC', 'Colangiocito': 'Cholangiocyte',
      'Plasma': 'Plasma', 'T_NK': 'T/NK', 'B': 'B', 'pDC': 'pDC', 'Hepatocito': 'Hepatocyte'}
fig, axs = plt.subplots(2, 2, figsize=(W2 * 0.82, W2 * 0.68)); axs = axs.ravel(); plt.subplots_adjust(hspace=0.62, wspace=0.5)
ax = axs[0]; ax.axhspan(-A['null_p97.5'].max(), A['null_p97.5'].max(), color='0.9', lw=0, label='95% permutation null')
ax.plot(A.lag, A.acf_unadjusted, 'o-', color='k', ms=3, lw=0.9, label='Stage effect per gene')
ax.plot(A.lag, A.acf_composition_adjusted, 's--', color=OI['orange'], ms=3, lw=0.9, label='Adjusted for cell composition')
ax.set_xscale('log'); ax.set_xlabel('Distance (genes)'); ax.set_ylabel('Spatial autocorrelation'); ax.legend(); lab(ax, 'A')
ax = axs[1]; xx = range(len(D))
ax.plot(xx, D.r_t, 'o-', color='k', ms=3, lw=0.9, label='Concordance of stage effect')
ax.plot(xx, D.r_t_comp_adj, 'o--', color=OI['orange'], ms=3, lw=0.9, label='…adjusted for composition')
ax.plot(xx, D.coexpr, 's-', color=OI['blue'], ms=3, lw=0.9, label='Co-expression across biopsies')
ax.set_xticks(list(xx)); ax.set_xticklabels(['overlap', '<1', '1–10', '10–50', '50–100', '100–500', '>500'], rotation=45, ha='right')
ax.set_xlabel('Intergenic distance (kb)'); ax.set_ylabel('Correlation'); ax.axhline(0, color='0.6', lw=0.5); ax.set_ylim(-0.08, 0.52); ax.legend(loc='upper right'); lab(ax, 'B')
ax = axs[2]; SCs = SC.sort_values('acf_lag1'); co = CO.groupby('type')[['r_neighbours', 'r_random']].mean().reindex(SCs.type)
y = np.arange(len(SCs)); ax.barh(y, SCs.acf_lag1, color=OI['red'], lw=0, height=0.68)
ax.errorbar([0] * len(SCs), y, xerr=SCs['null_p97.5'], fmt='none', color='k', capsize=2, lw=0.6, label='95% null')
ax.set_yticks(y); ax.set_yticklabels([en[t] for t in SCs.type]); ax.set_xlabel('Autocorrelation of the cirrhosis-vs-healthy\neffect within each cell type'); ax.legend(loc='lower right'); lab(ax, 'C', -0.55)
ax = axs[3]; y = np.arange(len(co))
ax.barh(y + 0.19, co.r_neighbours, 0.38, color=OI['blue'], lw=0, label='Neighbouring genes')
ax.barh(y - 0.19, co.r_random, 0.38, color='0.65', lw=0, label='Random gene pairs')
ax.set_yticks(y); ax.set_yticklabels([en[t] for t in co.index]); ax.set_xlabel('Co-expression across single cells'); ax.legend(loc='lower right'); lab(ax, 'D', -0.55)
r = EQ.iloc[1]
ax.text(1.0, -0.42, f'Coupled neighbours share a same-direction liver eQTL more often than uncoupled pairs\nat matched distance: {100*r.frac_coupled:.0f}% vs {100*r.frac_uncoupled:.0f}%, odds ratio {r.MH_OR_distance_stratified:.2f} (95% CI {r.CI_low:.2f}–{r.CI_high:.2f}), P = {r.MH_p:.3f}',
        transform=ax.transAxes, ha='right', va='top', fontsize=6)
save(fig, 'Fig4_neighbourhood_coupling')
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
save(fig, 'Fig5_GWAS_loci_neighbourhoods')
# ======================================================= FIGURE 5 (composite with circos)
X, keep = pd.read_pickle(inter('expr.pkl')); de, de2 = pd.read_pickle(inter('de.pkl'))
T = pd.read_csv(tab('Table_S6c_stage_effect_per_gene_with_without_composition.csv')).set_index('gene_id')
V = pd.read_csv(tab('Table_S5c_contiguous_DE_neighbourhoods.csv')); V['chr'] = V.chr.astype(str)
Dor = pd.read_csv(tab('Table_S10c_dorothea_ABC_regulons_used.csv')); G1 = pd.read_csv(tab('Table_S10a_GSEA_hallmark_ordinal.csv')).set_index('Term')
G2 = pd.read_csv(tab('Table_S10b_GSEA_hallmark_composition_adjusted.csv')).set_index('Term'); TF = pd.read_csv(tab('Table_S10d_DoRothEA_TF_activity_vs_stage.csv')).set_index('tf')
Cn = pd.read_csv(tab('Table_S9a_candidate_targets_lineage_intrinsic.csv')); Cn['development_stage'] = Cn.development_stage.fillna('not curated')
K = keep.set_index('gene_id').join(de2[['stat']]).join(T[['t_stage_comp_adjusted']]); K['stat'] = K.stat.fillna(0); K['gene_name'] = K.gene_name.astype(str)

# -- circos rendered to a raster, then embedded
hubs = ['TNFA_SIGNALING_VIA_NFKB', 'EPITHELIAL_MESENCHYMAL_TRANSITION', 'INFLAMMATORY_RESPONSE', 'ANGIOGENESIS', 'P53_PATHWAY', 'OXIDATIVE_PHOSPHORYLATION', 'BILE_ACID_METABOLISM']
short = {'TNFA_SIGNALING_VIA_NFKB': 'TNFα/NF-κB', 'EPITHELIAL_MESENCHYMAL_TRANSITION': 'EMT', 'INFLAMMATORY_RESPONSE': 'Inflammation',
         'ANGIOGENESIS': 'Angiogenesis', 'P53_PATHWAY': 'p53', 'OXIDATIVE_PHOSPHORYLATION': 'OXPHOS', 'BILE_ACID_METABOLISM': 'Bile acid'}
lead = {r.Term: set(str(r.Lead_genes).split(';')) for r in G1.reset_index().itertuples()}; nes = G1.NES.to_dict()
hcol = dict(zip(hubs, [OI['orange'], OI['blue'], OI['green'], OI['red'], OI['purple'], OI['yellow'], OI['sky']]))
links = []
for _, v in V.iterrows():
    genes = set(str(v.genes).split(',')); mid = (v.start_grid + v.end_grid) / 2
    for h in hubs:
        if genes & lead.get(h, set()): links.append({'chr': v.chr, 'pos': mid, 'hub': h})
L = pd.DataFrame(links)
tot = sum(N.values()); hubsize = int(0.075 * tot)
nfkb = set(Dor[Dor.tf.isin(['NFKB1', 'RELA'])].target); cand = set(Cn.gene_name)
circos = Circos({**{c: N[c] for c in CHR}, **{short[h]: hubsize for h in hubs}}, space=1.4, start=0, end=352)
for sec in circos.sectors:
    c = sec.name
    if c in CHR:
        sec.text(c, r=110, size=6.5)
        k = K[K.chr == c].sort_values('grid_index'); xx = k.grid_index.values.astype(float); yy = k.stat.values
        tr = sec.add_track((84, 100), r_pad_ratio=0.05); tr.axis(fc='#f4f4f4', ec='none')
        tr.bar(xx[yy > 0], np.clip(yy[yy > 0], 0, 12), width=1, color=OI['red'], lw=0, vmin=0, vmax=12, bottom=0)
        tr.bar(xx[yy < 0], np.clip(-yy[yy < 0], 0, 12), width=1, color=OI['blue'], lw=0, vmin=0, vmax=12, bottom=0)
        ya = k.t_stage_comp_adjusted.fillna(0).values
        tr2 = sec.add_track((72, 82), r_pad_ratio=0.05); tr2.axis(fc='#fbfbfb', ec='none')
        tr2.bar(xx[ya > 0], np.clip(ya[ya > 0], 0, 6), width=1, color='#b05a2c', lw=0, vmin=0, vmax=6, bottom=0)
        tr2.bar(xx[ya < 0], np.clip(-ya[ya < 0], 0, 6), width=1, color='#2c5f8c', lw=0, vmin=0, vmax=6, bottom=0)
        tr3 = sec.add_track((64, 70)); tr3.axis(fc='white', ec='0.85', lw=0.3)
        for _, v in V[V.chr == c].iterrows(): tr3.rect(v.start_grid, v.end_grid + 1, fc=OI['red'] if v.direction == 'up' else OI['blue'], ec='none')
        gi = k[k.gene_name.isin(nfkb)].grid_index.values; tr3.scatter(gi, [0.5] * len(gi), s=0.7, color=OI['orange'], vmin=0, vmax=1, lw=0)
        for g_ in k[k.gene_name.isin(cand)].grid_index.values: tr3.line([g_, g_], [0, 1], color=OI['green'], lw=0.55, vmin=0, vmax=1)
    else:
        h = [hh for hh in hubs if short[hh] == c][0]
        sec.text(f'{c}\nNES {nes.get(h, 0):+.1f}', r=113, size=6, color=hcol[h], fontweight='bold')
        tr = sec.add_track((64, 100)); tr.axis(fc=hcol[h], ec='none', alpha=0.32)
for _, r in L.iterrows():
    hs = short[r.hub]; tsec = circos.get_sector(hs); e = tsec.size * ((CHR.index(r.chr) + 0.5) / len(CHR))
    circos.link((r.chr, r.pos - 1.5, r.pos + 1.5), (hs, e - hubsize * 0.014, e + hubsize * 0.014), color=hcol[r.hub], alpha=0.45, lw=0.25, r1=63, r2=63)


fig = plt.figure(figsize=(W2, W2 * 1.08)); gs = fig.add_gridspec(2, 6, height_ratios=[1.45, 1.0], hspace=0.30, wspace=2.1)
axc = fig.add_subplot(gs[0, :], polar=True); circos.plotfig(ax=axc)
hdc = [Patch(fc=OI['red'], label='Stage effect per gene, up'), Patch(fc=OI['blue'], label='…down (outer: DESeq2;\ninner: composition-adjusted)'),
       Line2D([], [], marker='o', ls='', color=OI['orange'], ms=3.5, label='NF-κB target'), Line2D([], [], color=OI['green'], lw=1.2, label='Lineage-intrinsic target'),
       Patch(fc='0.5', label='DE neighbourhood (≥3 genes)'), Line2D([], [], color='0.5', lw=0.9, label='Link: neighbourhood →\nHallmark pathway (NES)')]
axc.legend(handles=hdc, loc='center left', bbox_to_anchor=(0.80, 0.90), ncol=1, fontsize=5.2, columnspacing=1.0, handlelength=1.3, labelspacing=0.32)
axc.text(-0.02, 1.02, 'A', transform=axc.transAxes, fontsize=9, fontweight='bold', va='top')
# B — GSEA two layers
ax = fig.add_subplot(gs[1, 0:2])
paths = ['TNFA_SIGNALING_VIA_NFKB', 'EPITHELIAL_MESENCHYMAL_TRANSITION', 'INFLAMMATORY_RESPONSE', 'ANGIOGENESIS', 'P53_PATHWAY', 'ADIPOGENESIS', 'FATTY_ACID_METABOLISM', 'OXIDATIVE_PHOSPHORYLATION', 'BILE_ACID_METABOLISM']
paths = [p for p in paths if p in G1.index and p in G2.index]
nm = lambda p: p.replace('_', ' ').title().replace('Tnfa Signaling Via Nfkb', 'TNFα / NF-κB').replace('Epithelial Mesenchymal Transition', 'EMT').replace('P53', 'p53').replace('Oxidative Phosphorylation', 'OXPHOS')
y = np.arange(len(paths))
for i, p in enumerate(paths): ax.plot([G1.loc[p, 'NES'], G2.loc[p, 'NES']], [i, i], color='0.75', lw=1)
ax.scatter([G1.loc[p, 'NES'] for p in paths], y, s=16, color='0.35', zorder=3, label='Raw stage statistic')
ax.scatter([G2.loc[p, 'NES'] for p in paths], y, s=16, color=OI['orange'], zorder=3, label='Composition-adjusted')
ax.axvline(0, color='k', lw=0.5); ax.set_yticks(y); ax.set_yticklabels([nm(p) for p in paths], fontsize=6); ax.invert_yaxis()
ax.set_xlabel('GSEA normalised enrichment score'); ax.set_xlim(-3.6, 4.6); ax.legend(loc='lower left', fontsize=5.3, handletextpad=0.4); lab(ax, 'B', -0.72, 1.1)
# C — TF activity
ax = fig.add_subplot(gs[1, 2:4]); tfs = [t for t in ['HNF4A', 'HNF1A', 'FOXA1', 'PPARA', 'NFKB1', 'RELA', 'TP53', 'SMAD3', 'SMAD4', 'KLF4', 'FOXM1'] if t in TF.index]
y = np.arange(len(tfs))
for i, t in enumerate(tfs): ax.plot([TF.loc[t, 't_stage'], TF.loc[t, 't_stage_comp_adjusted']], [i, i], color='0.75', lw=1)
ax.scatter([TF.loc[t, 't_stage'] for t in tfs], y, s=16, color='0.35', zorder=3)
ax.scatter([TF.loc[t, 't_stage_comp_adjusted'] for t in tfs], y, s=16, color=OI['orange'], zorder=3)
ax.axvline(0, color='k', lw=0.5); ax.set_yticks(y); ax.set_yticklabels(tfs, style='italic', fontsize=6); ax.invert_yaxis()
ax.axhline(3.5, color='0.8', lw=0.5, ls=':'); ax.set_xlabel('t of stage on inferred TF activity')
ax.set_xlim(-9.5, 20.5)
ax.text(0.99, 0.97, 'Hepatocyte identity', transform=ax.transAxes, ha='right', va='top', fontsize=5.8, color='0.45')
ax.text(0.99, 0.05, 'Intrinsic stress response', transform=ax.transAxes, ha='right', va='bottom', fontsize=5.8, color='0.45'); lab(ax, 'C', -0.35, 1.12)
# D — targets
ax = fig.add_subplot(gs[1, 4:6]); colmap = {'clinical, liver': OI['red'], 'clinical, inflammasome': OI['orange'], 'clinical, other fibrosis': OI['orange'],
                                          'clinical, other': OI['orange'], 'approved, other': OI['orange'], 'preclinical': OI['sky'], 'tool compounds': OI['sky'],
                                          'biomarker': OI['green'], 'none': '0.7', 'not curated': '0.85'}
top = Cn.head(22)
ax.scatter(top.t_bulk_comp_adjusted, top.t_sc_max, s=22, c=[colmap.get(s, '0.8') for s in top.development_stage], lw=0.3, edgecolor='k', zorder=3)
texts = [ax.text(r.t_bulk_comp_adjusted, r.t_sc_max, r.gene_name, fontsize=5.2, style='italic') for _, r in top.iterrows()]
ax.set_xlim(0.4, 10.0); ax.set_ylim(1.4, 6.5); adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='0.5', lw=0.3), expand=(1.6, 1.9), force_text=(0.6, 0.9))
hd2 = [Line2D([], [], marker='o', ls='', color=c, markeredgecolor='k', markeredgewidth=0.3, ms=4, label=l) for l, c in
       [('Clinical, liver', OI['red']), ('Clinical, other indication', OI['orange']), ('Preclinical / tool compound', OI['sky']), ('Biomarker only', OI['green']), ('No modulator', '0.7')]]
ax.legend(handles=hd2, loc='upper left', ncol=1, fontsize=5.0, handletextpad=0.4, labelspacing=0.28)
ax.set_xlabel('Bulk stage effect, composition-adjusted (t)'); ax.set_ylabel('Maximal within-cell-type effect (t)'); lab(ax, 'D', -0.3, 1.1)
save(fig, 'Fig6_integrated_molecular_map')
T1 = Cn.head(25)[['gene_name', 'cell_type_up_sc', 't_bulk_stage', 't_bulk_comp_adjusted', 't_sc_max', 'coupled_neighbour', 'target_class', 'drug_landscape', 'development_stage']].copy()
T1.columns = ['Gene', 'Cell type (scRNA-seq)', 'Bulk t (stage)', 'Bulk t (composition-adjusted)', 'Max within-cell-type t', 'Cis-coupled neighbour', 'Target class', 'Agents / status', 'Development stage']
T1['Cell type (scRNA-seq)'] = T1['Cell type (scRNA-seq)'].str.replace('Mesenquima_HSC', 'stellate').str.replace('Fagocito_mononuclear', 'macrophage').str.replace('Endotelio', 'endothelium').str.replace('Colangiocito', 'cholangiocyte')
T1.round(1).to_csv(os.path.join(FO, 'Table1_targets.csv'), index=False)
print('JHEP figure set written to', FO)
