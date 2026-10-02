"""16 — Journal of Hepatology main figures: Fig 1 (composition, thresholds, F3 bimodality),
Fig 2 (mechanotransduction switch at F4), Fig 6 (GWAS loci in coupled neighbourhoods).
Figures 3-5 are adapted from scripts 09/10. Run after 01-08, 15.
"""
import os, numpy as np, pandas as pd, statsmodels.api as sm, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, scienceplots, seaborn as sns
from sklearn.mixture import GaussianMixture
from common import *

plt.style.use(['science', 'nature', 'no-latex'])
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial'], 'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
                     'xtick.top': False, 'ytick.right': False, 'xtick.minor.visible': False, 'ytick.minor.visible': False, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False, 'pdf.fonttype': 42})
MM = 1 / 25.4; W = 180 * MM
OI = {'blue': '#0072B2', 'orange': '#E69F00', 'green': '#009E73', 'red': '#D55E00', 'purple': '#CC79A7', 'sky': '#56B4E9', 'grey': '#7F7F7F'}
ST = STAGES; pal = dict(zip(ST, sns.color_palette('viridis', 6)))
def lab(ax, t, dx=-0.17, dy=1.12): ax.text(dx, dy, t, transform=ax.transAxes, fontsize=10, fontweight='bold', va='top')
def save(fig, name): fig.savefig(os.path.join(FIG, name + '.pdf'), bbox_inches='tight'); fig.savefig(os.path.join(FIG, name + '.png'), dpi=400, bbox_inches='tight'); plt.close(fig)

C = pd.read_csv(tab('Table_S6b_composition_scores_per_sample.csv'), index_col=0)
S = pd.read_csv(tab('Table_S15b_mechanics_scores_per_sample.csv'), index_col=0)
B = pd.read_csv(tab('Table_S7a_threshold_vs_linear.csv')); Mx = pd.read_csv(tab('Table_S15_mechanics_scores_thresholds.csv'))
EN = ['Hepatocyte', 'Stellate/myofibroblast', 'Cholangiocyte', 'Macrophage', 'Lymphocyte', 'Endothelium']
cmap_ct = dict(zip(MARKERS, [OI['green'], OI['red'], OI['purple'], OI['orange'], OI['blue'], OI['grey']]))

# ---------------- Figure 1
fig = plt.figure(figsize=(W, W * 0.60)); gs = fig.add_gridspec(2, 2, hspace=0.6, wspace=0.42)
ax = fig.add_subplot(gs[0, 0])
for ct, en in zip(MARKERS, EN):
    g = C[C.estadio != 'Control'].groupby('estadio')[ct].agg(['mean', 'sem']).reindex(ST)
    ax.errorbar(range(6), g['mean'], yerr=g['sem'], fmt='o-', color=cmap_ct[ct], ms=3, lw=1.1, capsize=1.5, elinewidth=.7, label=en)
ax.axhline(0, color='0.75', lw=.6); ax.set_xticks(range(6)); ax.set_xticklabels(ST); ax.set_xlabel('Fibrosis stage'); ax.set_ylabel('Composition score (marker z)')
ax.set_ylim(-1.0, 1.6); ax.legend(ncol=2, handlelength=1.1, columnspacing=.8, fontsize=6, loc='upper left'); lab(ax, 'A')
ax = fig.add_subplot(gs[0, 1])
for s in ST:
    g = C[C.estadio == s]; ax.scatter(g.Hepatocito, g.offset, s=7, color=pal[s], alpha=.85, lw=0, label=s)
ax.set_xlabel('Hepatocyte score'); ax.set_ylabel('Global transcriptome variance'); ax.set_xlim(-2.3, 2.0)
ax.text(.03, .97, 'Stage effect: t = \u22128.0 \u2192 +1.4\nafter composition adjustment', transform=ax.transAxes, va='top', fontsize=6.5)
ax.legend(loc='lower right', ncol=1, markerscale=1.6, title='Stage', title_fontsize=6.5, fontsize=6.5); lab(ax, 'B')
ax = fig.add_subplot(gs[1, 0])
H = B.set_index('variable')[[c for c in B.columns if c.startswith('dAIC_step_at_')]]; H.columns = ['\u2265F1', '\u2265F2', '\u2265F3', '\u2265F4']
H = H.loc[list(MARKERS)]; H.index = EN
im = ax.imshow(H.values.T, cmap='RdBu_r', vmin=-15, vmax=15, aspect='auto')
ax.set_yticks(range(4)); ax.set_yticklabels(H.columns); ax.set_xticks(range(len(H))); ax.set_xticklabels(H.index, rotation=35, ha='right', fontsize=6.5); ax.set_ylabel('Step at stage')
for i in range(H.shape[0]):
    for j in range(H.shape[1]): ax.text(i, j, f'{H.values[i, j]:.0f}', ha='center', va='center', fontsize=6, color='white' if abs(H.values[i, j]) > 9 else 'k')
cb = plt.colorbar(im, ax=ax, fraction=.035, pad=.02); cb.set_label('\u0394AIC (step \u2212 linear)', fontsize=6.5); cb.ax.tick_params(labelsize=6)
ax.spines[['left', 'bottom']].set_visible(False); ax.tick_params(length=0); lab(ax, 'C', -0.22)
ax = fig.add_subplot(gs[1, 1])
f3 = C[C.estadio == 'F3'].Hepatocito.values; f4 = C[C.estadio == 'F4'].Hepatocito.values
ax.hist(f3, bins=18, color=pal['F3'], alpha=.75, lw=0, label='F3 (n=%d)' % len(f3), density=True)
ax.hist(f4, bins=12, histtype='step', color=pal['F4'], lw=1.3, label='F4 (n=%d)' % len(f4), density=True)
g2 = GaussianMixture(2, random_state=0, n_init=5).fit(f3.reshape(-1, 1)); xs = np.linspace(f3.min() - .3, f3.max() + .3, 200)
for k in range(2):
    ax.plot(xs, g2.weights_[k] * np.exp(-(xs - g2.means_[k, 0]) ** 2 / (2 * g2.covariances_[k, 0, 0])) / np.sqrt(2 * np.pi * g2.covariances_[k, 0, 0]), '--', color='k', lw=.8)
ax.set_xlabel('Hepatocyte score'); ax.set_ylabel('Density'); ax.legend(fontsize=6.5)
ax.text(.03, .97, '\u0394BIC = 12.8 for two components\nminor component %.0f%%' % (100 * g2.weights_.min()), transform=ax.transAxes, va='top', fontsize=6.5); lab(ax, 'D')
save(fig, 'JHEP_Fig1_composition_thresholds')

# ---------------- Figure 2: mechanotransduction switch
PROG = {'matrix_crosslinking': 'Matrix crosslinking', 'YAP_TAZ_targets': 'YAP/TAZ targets', 'mechanosensing': 'Mechanosensing machinery', 'TGFb_SMAD': 'TGF-\u03b2/SMAD', 'senescence_SASP': 'Senescence/SASP'}
fig = plt.figure(figsize=(W, W * 0.42)); gs = fig.add_gridspec(1, 3, wspace=0.75, width_ratios=[1.2, 1, 1.1])
ax = fig.add_subplot(gs[0, 0]); cols = [OI['grey'], OI['red'], OI['orange'], OI['blue'], OI['green']]
for (k, en), c in zip(PROG.items(), cols):
    g = S.groupby('estadio')[k].agg(['mean', 'sem']).reindex(ST); ax.errorbar(range(6), g['mean'], yerr=g['sem'], fmt='o-', color=c, ms=3, lw=1.1, capsize=1.5, elinewidth=.7, label=en)
ax.axhline(0, color='0.75', lw=.6); ax.set_xticks(range(6)); ax.set_xticklabels(ST); ax.set_xlabel('Fibrosis stage'); ax.set_ylabel('Programme score (marker z)')
ax.legend(fontsize=6.3, loc='upper left', handlelength=1.1); lab(ax, 'A', -0.22)
ax = fig.add_subplot(gs[0, 1])
x = np.arange(len(PROG)); raw = [Mx.set_index('score').loc[k, 't_step_F4'] for k in PROG]; adj = [Mx.set_index('score').loc[k, 't_step_F4_comp_adj'] for k in PROG]
ax.barh(x + .19, raw, .36, color='0.65', lw=0, label='Unadjusted'); ax.barh(x - .19, adj, .36, color=OI['red'], lw=0, label='Adjusted for composition')
ax.axvline(2, color='0.4', ls='--', lw=.7); ax.axvline(0, color='k', lw=.6)
ax.set_yticks(x); ax.set_yticklabels([PROG[k] for k in PROG], fontsize=6.5); ax.set_xlabel('t of a step at F4'); ax.set_xlim(0, 7.6); ax.legend(fontsize=6.5, loc='upper right'); lab(ax, 'B', -0.95)
ax = fig.add_subplot(gs[0, 2])
for s in ST:
    g = S[S.estadio == s]; ax.scatter(g.HSC, g.YAP_TAZ_targets, s=8, color=pal[s], lw=0, alpha=.85, label=s)
b, a = np.polyfit(S[S.estadio != 'F4'].HSC, S[S.estadio != 'F4'].YAP_TAZ_targets, 1); xs = np.linspace(S.HSC.min(), S.HSC.max(), 50)
ax.plot(xs, a + b * xs, color='k', lw=.9, ls='--', label='fit, Normal\u2013F3')
ax.set_xlabel('Stellate score'); ax.set_ylabel('YAP/TAZ target score'); ax.set_ylim(-1.8, 3.4); ax.legend(fontsize=6, ncol=3, loc='upper left', columnspacing=.6, handletextpad=.3); lab(ax, 'C', -0.28)
save(fig, 'JHEP_Fig2_mechanotransduction_switch')

# ---------------- Figure 6: GWAS loci
G = pd.read_csv(tab('Table_S16c_GWAS_catalog_enrichment_by_trait.csv'))
order = ['liver fibrosis/cirrhosis', 'MASLD/NAFLD', 'ALT (liver enzyme)', 'height (control)', 'educational attainment (control)']
nice = {'liver fibrosis/cirrhosis': 'Liver fibrosis\n/cirrhosis', 'MASLD/NAFLD': 'MASLD', 'ALT (liver enzyme)': 'ALT', 'height (control)': 'Height\n(control)', 'educational attainment (control)': 'Education\n(control)'}
G = G.set_index('trait_set').reindex(order)
fig, axs = plt.subplots(1, 2, figsize=(W * 0.78, W * 0.3)); plt.subplots_adjust(wspace=.42)
ax = axs[0]; colr = [OI['red'], OI['orange'], OI['sky'], '0.7', '0.7']
ax.bar(range(5), 100 * G.frac_in_coupled_pair, color=colr, lw=0, width=.65)
ax.axhline(100 * G.background.iloc[0] if 'background' in G else 21.5, color='k', ls='--', lw=.8, label='All genes')
ax.set_xticks(range(5)); ax.set_xticklabels([nice[t] for t in order], fontsize=6.3); ax.set_ylabel('% of locus genes in a\ncoupled neighbour pair'); ax.legend(fontsize=6.5); lab(ax, 'A', -0.25)
for i, (p, v) in enumerate(zip(G.P_coupled, G.frac_in_coupled_pair)):
    if p < 0.06: ax.text(i, 100 * v + .6, 'P=%.3f' % p, ha='center', fontsize=5.8)
ax = axs[1]
ax.bar(range(5), G.acf_windows, color=colr, lw=0, width=.65)
ax.errorbar(range(5), G.acf_random, yerr=2 * G.acf_random_sd, fmt='o', color='k', ms=3, capsize=2.5, lw=.8, label='Random windows \u00b1 2 SD')
ax.set_xticks(range(5)); ax.set_xticklabels([nice[t] for t in order], fontsize=6.3); ax.set_ylabel('Spatial autocorrelation of the\nstage effect around loci'); ax.legend(fontsize=6.5); lab(ax, 'B', -0.25)
for i, (p, v) in enumerate(zip(G.P_acf, G.acf_windows)):
    if p < 0.06: ax.text(i, v + .004, 'P=%.3f' % p, ha='center', fontsize=5.8)
save(fig, 'JHEP_Fig6_GWAS_neighbourhoods')
print('JHEP figures written to', FIG)
