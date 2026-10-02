"""16 — Journal of Hepatology figure set (final numbering). Run after scripts 01–15.
Fig 1  The positional framework (spectrum, GTEx tissue specificity, invariance)
Fig 2  Linear cell replacement with thresholds at F3–F4
Fig 3  Mechanotransduction switch at F4
Fig 4  Cis-coupled gene neighbourhoods within each cell type
Fig 5  Genome-wide circos map of the fibrotic response
Fig 6  Two layers (GSEA, TF activity) and lineage-intrinsic targets
Supplementary: MASLD and cirrhosis GWAS loci in coupled neighbourhoods
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
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial', 'Helvetica'],
                     'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 6.5,
                     'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
                     'xtick.major.size': 2.4, 'ytick.major.size': 2.4, 'lines.linewidth': 1.0,
                     'xtick.top': False, 'ytick.right': False, 'xtick.minor.visible': False, 'ytick.minor.visible': False,
                     'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
                     'legend.handlelength': 1.3, 'legend.borderaxespad': 0.3, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                     'savefig.transparent': False, 'figure.dpi': 200})
MM = 1 / 25.4; W2 = 183 * MM; FO = os.path.join(FIG, 'jhep'); os.makedirs(FO, exist_ok=True)
ST = ['Normal', 'F0', 'F1', 'F2', 'F3', 'F4']; pal = dict(zip(ST, sns.color_palette('viridis', 6)))
OI = {'blue': '#0072B2', 'orange': '#E69F00', 'green': '#009E73', 'red': '#D55E00', 'purple': '#CC79A7', 'sky': '#56B4E9', 'grey': '#7F7F7F', 'yellow': '#F0E442'}
def save(fig, name):
    fig.savefig(f'{FO}/{name}.pdf', bbox_inches='tight'); fig.savefig(f'{FO}/{name}.png', dpi=400, bbox_inches='tight')
    fig.savefig(f'{FO}/{name}.tiff', dpi=300, bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'}); plt.close(fig)
def lab(ax, s, dx=-0.2, dy=1.15): ax.text(dx, dy, s, transform=ax.transAxes, fontsize=8, fontweight='bold', va='top', ha='left')
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

fig = plt.figure(figsize=(W2, W2 * 0.84)); gs = fig.add_gridspec(3, 6, height_ratios=[1.0, 0.92, 0.92], hspace=0.88, wspace=2.1)
# a — genome-wide consensus spectrum
ax = fig.add_subplot(gs[0, :])
m = U.set_index(['chr', 'k']).mean_log_w; frac = U.set_index(['chr', 'k']).frac_samples_w_gt3
idx = [(c, k) for c, k in zip(U.chr, U.k)]; x = np.arange(len(U)); chr_of = U.chr.values; u = U.universal.values
for i, c in enumerate(CHR):
    sel = chr_of == c; ax.plot(x[sel], np.exp(U.mean_log_w.values[sel] + EULER), color=OI['grey'] if i % 2 else '#4d4d4d', lw=0.35)
ax.axhline(1, color='k', lw=0.5, ls=':'); ax.scatter(x[u], np.exp(U.mean_log_w.values[u] + EULER), s=3.5, color=OI['red'], lw=0, zorder=5,
                                                     label=f'Universal peaks: >3× background in ≥90% of 437 biopsies (n={u.sum()})')
for c, k, nm, dy in [('7', 131, 'chr7 k=131\nperiod ≈ 7 genes', 12), ('1', 117, 'chr1 k=117\n≈18 genes', 12), ('19', 565, 'chr19 k=565\n≈2.6', 12), ('12', 201, 'chr12 k=201\n≈5', 12), ('2', 168, 'chr2 k=168\n≈7.5', 12)]:
    i = idx.index((c, k)); ax.annotate(nm, (x[i], np.exp(U.mean_log_w.values[i] + EULER)), fontsize=6, xytext=(0, dy), textcoords='offset points',
                                       ha='center', arrowprops=dict(arrowstyle='-', lw=0.4, color='0.4'))
xt = [np.where(chr_of == c)[0].mean() for c in CHR]
ax.set_xticks(xt); ax.set_xticklabels([c if (c == 'X' or int(c) % 2 == 1) else '' for c in CHR]); ax.set_yscale('log'); ax.set_ylim(0.3, 110)
ax.set_ylabel('Power / 1/f background'); ax.set_xlabel('Chromosome (spatial frequency increases within each; odd chromosomes labelled)', labelpad=2)
ax.legend(loc='lower left', bbox_to_anchor=(0, 1.01)); lab(ax, 'A', -0.12, 1.24)
# B — GTEx tissue specificity (full-width second row)
ax = fig.add_subplot(gs[1, 0:3]); e = S14e.sort_values('odds_vs_background', ascending=False)
ax.bar(range(len(e)), np.log2(e.odds_vs_background), color=[OI['red'] if t == 'liver' else '0.6' for t in e.tissue], lw=0, width=0.68)
ax.axhline(0, color='k', lw=0.5); ax.set_xticks(range(len(e))); ax.set_xticklabels([pretty[t] for t in e.tissue], fontsize=6, rotation=35, ha='right')
ax.set_ylabel('Enrichment of liver-biopsy\npeaks (log$_2$ odds)')
ax.set_xlabel('GTEx tissue', labelpad=1); lab(ax, 'B', -0.24, 1.18)
# C — period density
ax = fig.add_subplot(gs[1, 3:6]); per = U.period.values
ax.hexbin(np.log10(per), np.log10(np.exp(U.mean_log_w.values + EULER)), gridsize=38, cmap='Greys', mincnt=1, linewidths=0)
ax.scatter(np.log10(per[u]), np.log10(np.exp(U.mean_log_w.values[u] + EULER)), s=3.5, color=OI['red'], lw=0, zorder=5)
ax.set_xticks([0.3, 1, 2, 3]); ax.set_xticklabels(['2', '10', '100', '1000']); ax.set_yticks([0, 0.5, 1]); ax.set_yticklabels(['1', '3', '10'])
ax.set_xlabel('Period (genes)'); ax.set_ylabel('Mean power / background'); lab(ax, 'C', -0.20, 1.18)
# D — reproducibility
ax = fig.add_subplot(gs[2, 0:2]); ax.hist(U.frac_samples_w_gt3, bins=50, color=OI['grey'], lw=0); ax.set_yscale('log')
ax.axvline(0.9, color=OI['red'], ls='--', lw=0.8); ax.text(0.88, ax.get_ylim()[1] * 0.5, f'n = {u.sum()}', ha='right', fontsize=6, color=OI['red'])
ax.set_xlabel('Fraction of biopsies with power >3×'); ax.set_ylabel('Number of frequencies'); lab(ax, 'D', -0.36, 1.18)
# E — the architecture is the same in normal liver and cirrhosis (all 10,007 frequencies)
ax = fig.add_subplot(gs[2, 2:4]); lw_ = np.log(Wl); nm_ = lw_[M.index[M.estadio == 'Normal']].mean(1); f4_ = lw_[M.index[M.estadio == 'F4']].mean(1)
xa = np.log10(np.exp(nm_ + EULER)); ya = np.log10(np.exp(f4_ + EULER)); uu = U.set_index(['chr', 'k']).universal.reindex(lw_.index).fillna(False).values
ax.hexbin(xa, ya, gridsize=38, cmap='Greys', mincnt=1, linewidths=0, bins='log')
ax.scatter(xa[uu], ya[uu], s=3, color=OI['red'], lw=0, zorder=5); lo_, hi_ = -1.3, 1.35; ax.plot([lo_, hi_], [lo_, hi_], color='k', lw=0.5, ls=':')
ax.set_xlim(lo_, hi_); ax.set_ylim(lo_, hi_); ax.set_xticks([-1, 0, 1]); ax.set_xticklabels(['0.1', '1', '10']); ax.set_yticks([-1, 0, 1]); ax.set_yticklabels(['0.1', '1', '10'])
ax.set_xlabel('Power / background, Normal'); ax.set_ylabel('Power / background, F4'); ax.text(0.04, 0.96, f'r = {np.corrcoef(nm_, f4_)[0, 1]:.2f}', transform=ax.transAxes, va='top', fontsize=6.5)
lab(ax, 'E', -0.34, 1.18)
# F — universal peaks keep their amplitude at every stage
ax = fig.add_subplot(gs[2, 4:6]); uidx = pd.MultiIndex.from_frame(U[U.universal][['chr', 'k']].astype({'chr': str}))
upw = np.exp(np.log(Wl.loc[uidx]).mean())
bp = ax.boxplot([upw[M.index[M.estadio == st]].values for st in ST], tick_labels=ST, patch_artist=True, showfliers=False, widths=0.62,
                medianprops=dict(color='white', lw=1.1), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6), boxprops=dict(lw=0))
for p_, st in zip(bp['boxes'], ST): p_.set_facecolor(pal[st])
from scipy.stats import spearmanr as _sr
_s = M.index[M.estadio != 'Control']; rho_, pv_ = _sr(M.loc[_s, 'orden'], upw[_s])
ax.set_ylabel('Mean power at the 186\nuniversal peaks (× background)'); ax.set_xlabel('Fibrosis stage'); ax.tick_params(axis='x', labelsize=6.5)
ax.set_ylim(3.0, 5.6); ax.axhline(1, color='0.6', lw=0.5, ls=':')
ax.text(0.04, 0.96, f'Stage: ρ = {rho_:.2f}, P = {pv_:.2f}', transform=ax.transAxes, ha='left', va='top', fontsize=6.5); lab(ax, 'F', -0.36, 1.18)
save(fig, 'Fig1_design_and_positional_framework')
# ================= Fig 2 — fibrosis as linear cell replacement, with thresholds; composition does not define cirrhosis
from scipy.stats import spearmanr as _spr
T6 = pd.read_csv(tab('Table_S6c_stage_effect_per_gene_with_without_composition.csv'))
F4L = pd.read_csv(tab('Table_S7c_F4_likeness_per_sample.csv'), index_col=0)
fig = plt.figure(figsize=(W2, W2 * 0.9)); gs2 = fig.add_gridspec(3, 3, height_ratios=[1, 0.62, 1], hspace=0.7, wspace=0.6)
# A — composition trajectories with model fits (linear dotted; trend + step at F4 solid where ΔAIC > 2)
ax = fig.add_subplot(gs2[0, 0]); xs6 = np.arange(6); stepF4 = B.set_index('variable')['dAIC_step_at_F4']
for ct in CT_COL:
    g = C.groupby('orden')[ct].agg(['mean', 'sem']).reindex(range(6)); ax.errorbar(xs6, g['mean'], yerr=g['sem'], fmt='o', color=CT_COL[ct], ms=2.4, capsize=1.2, elinewidth=0.6, label=EN[ct])
    d_ = C[[ct, 'orden']].dropna(); o = d_.orden.values.astype(float); y = d_[ct].values
    if stepF4.get(ct, 0) > 2:
        Xs = np.column_stack([np.ones_like(o), o, (o >= 5).astype(float)]); bs = np.linalg.lstsq(Xs, y, rcond=None)[0]
        ax.plot(xs6, bs[0] + bs[1] * xs6 + bs[2] * (xs6 >= 5), color=CT_COL[ct], lw=1.0)
    else:
        b1 = np.polyfit(o, y, 1); ax.plot(xs6, np.polyval(b1, xs6), color=CT_COL[ct], lw=0.8, ls=':')
ax.set_xticks(xs6); ax.set_xticklabels(ST, fontsize=6.5); ax.axhline(0, color='0.8', lw=0.4); ax.set_ylabel('Composition score (marker z)'); ax.set_xlabel('Fibrosis stage'); ax.set_ylim(-1.05, 1.9)
h_, l_ = ax.get_legend_handles_labels(); l_ = [x.replace('Stellate / myofibroblast', 'Stellate') for x in l_]
leg1 = ax.legend(h_, l_, ncol=2, fontsize=6, loc='upper left', handlelength=1.0, columnspacing=0.6, labelspacing=0.22, handletextpad=0.3); ax.add_artist(leg1)
ax.legend([Line2D([], [], color='k', lw=0.8, ls=':'), Line2D([], [], color='k', lw=1.0)], ['Linear trend', 'Trend + step at F4'],
          fontsize=6, loc='lower right', handlelength=1.6, labelspacing=0.22, handletextpad=0.4); lab(ax, 'A', -0.3, 1.16)
# B — per-gene stage effect before vs after composition adjustment
ax = fig.add_subplot(gs2[0, 1]); xg, yg = T6.t_stage.values, T6.t_stage_comp_adjusted.values
ax.hexbin(xg, yg, gridsize=45, cmap='Greys', mincnt=1, linewidths=0, bins='log', extent=(-20, 20, -12, 12))
ax.plot([-20, 20], [-20, 20], color='0.5', lw=0.5, ls=':'); ax.axhline(0, color='0.7', lw=0.4); ax.axvline(0, color='0.7', lw=0.4); ax.set_xlim(-20, 20); ax.set_ylim(-12, 12)
lbl_pos = {'THBS2': (12, 10.6), 'IL32': (-3, 10.6), 'TIMP1': (17, 5.5), 'CYP2C19': (-15, 7.5), 'KRT7': (14, -10), 'COL1A2': (16, -5.5)}
for gname, (tx, ty) in lbl_pos.items():
    r_ = T6[T6.gene_name == gname]
    if len(r_): ax.annotate(gname, (r_.t_stage.iloc[0], r_.t_stage_comp_adjusted.iloc[0]), xytext=(tx, ty), fontsize=6, style='italic', ha='center', arrowprops=dict(arrowstyle='-', lw=0.35, color='0.4'))
ve = 100 * (1 - np.var(yg) / np.var(xg))
ax.text(0.03, 0.04, f'SD of t {np.std(xg):.1f} → {np.std(yg):.1f}\ncomposition explains {ve:.0f}%', transform=ax.transAxes, ha='left', fontsize=6.3)
ax.set_xlabel('Stage effect per gene (t)'); ax.set_ylabel('Composition-adjusted (t)'); lab(ax, 'B', -0.3, 1.16)
# C — global variance follows hepatocyte loss
ax = fig.add_subplot(gs2[0, 2])
for st in ST: g = C[C.estadio == st]; ax.scatter(g.Hepatocito, g.offset, s=4, color=pal[st], alpha=.85, lw=0, label=st)
ax.set_xlabel('Hepatocyte score'); ax.set_ylabel('Global transcriptome variance'); ax.set_xlim(-2.3, 3.0); ax.set_ylim(3.09, 3.40)
ax.text(0.98, 0.04, 'Stage effect\nt = −8.0 → +1.4\nafter adjustment', transform=ax.transAxes, fontsize=6.2, ha='right', va='bottom')
ax.legend(loc='upper right', fontsize=6, markerscale=1.5, ncol=1, handletextpad=0.2, labelspacing=0.3, title='Stage', title_fontsize=6); lab(ax, 'C', -0.3, 1.16)

# D — threshold map (ΔAIC)
ax = fig.add_subplot(gs2[1, :]); H = B.set_index('variable')[[c for c in B.columns if c.startswith('dAIC_step_at_')]]; H.columns = ['≥F1', '≥F2', '≥F3', '≥F4']
H = H.loc[['Hepatocito', 'HSC', 'Colangiocito', 'Macrofago', 'Linfocito', 'Endotelio']]; H.index = [EN[i] for i in H.index]; H = H.T
im = ax.imshow(H.values, cmap='RdBu_r', vmin=-15, vmax=15, aspect='auto')
ax.set_yticks(range(4)); ax.set_yticklabels(['Step from ' + c for c in H.index], fontsize=6.5); ax.set_xticks(range(H.shape[1])); ax.set_xticklabels([c.replace(' / myofibroblast', '') for c in H.columns], fontsize=6.5)
for i in range(H.shape[0]):
    for j in range(H.shape[1]): ax.text(j, i, f'{H.values[i, j]:.0f}', ha='center', va='center', fontsize=6.5, color='white' if abs(H.values[i, j]) > 9 else 'k')
cb = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.012); cb.set_label('ΔAIC, step − linear\n(>2 favours a threshold)', fontsize=6.3); cb.ax.tick_params(labelsize=6)
ax.spines[['left', 'bottom']].set_visible(False); ax.tick_params(length=0); lab(ax, 'D', -0.085, 1.22)

# F — composition tracks disease activity as much as fibrosis
ax = fig.add_subplot(gs2[2, 0]); cts = list(CT_COL); rs = [_spr(C[c], C.orden, nan_policy='omit')[0] for c in cts]; rn = [_spr(C[c], C.nas, nan_policy='omit')[0] for c in cts]
yy = np.arange(len(cts)); ax.barh(yy + 0.19, rs, 0.38, color='0.35', lw=0, label='Fibrosis stage'); ax.barh(yy - 0.19, rn, 0.38, color=OI['sky'], lw=0, label='Activity (NAS)')
ax.set_yticks(yy); ax.set_yticklabels([EN[c].replace(' / myofibroblast', '') for c in cts], fontsize=6.5); ax.axvline(0, color='k', lw=0.5); ax.invert_yaxis()
ax.set_xlabel('Spearman ρ with composition score'); ax.legend(fontsize=6, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, columnspacing=1.0); lab(ax, 'E', -0.45, 1.16)
# G — composition does not separate F3 from F4
ax = fig.add_subplot(gs2[2, 1]); F4L_ = F4L[F4L.estadio.isin(ST)]
bp = ax.boxplot([F4L_[F4L_.estadio == st].F4_likeness.values for st in ST], tick_labels=ST, patch_artist=True, showfliers=False, widths=0.6,
                medianprops=dict(color='white', lw=1.1), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6), boxprops=dict(lw=0))
for p_, st in zip(bp['boxes'], ST): p_.set_facecolor(pal[st])
thr = F4L_[F4L_.estadio == 'F4'].F4_likeness.quantile(0.25); ax.axhline(thr, color=OI['red'], ls='--', lw=0.7)
frac3 = (F4L_[F4L_.estadio == 'F3'].F4_likeness > thr).mean()
ax.text(0.03, 0.96, f'{100*frac3:.0f}% of F3 above the\nF4 lower quartile', transform=ax.transAxes, fontsize=6.2, va='top', color=OI['red'])
ax.set_ylabel('F4-likeness of composition'); ax.set_xlabel('Fibrosis stage'); ax.tick_params(axis='x', labelsize=6.5); lab(ax, 'F', -0.32, 1.16)
# H — …but only F4 has switched on the mechanical programme
ax = fig.add_subplot(gs2[2, 2]); f3_ = F4L_[F4L_.estadio == 'F3']; grp = {'F3, other': f3_[f3_.F4_likeness <= thr], 'F3, F4-like\ncomposition': f3_[f3_.F4_likeness > thr], 'F4': F4L_[F4L_.estadio == 'F4']}
cols_h = [pal['F3'], OI['red'], pal['F4']]; rng_ = np.random.default_rng(1)
for i, (k_, d_) in enumerate(grp.items()):
    v = d_.YAP_TAZ_targets.values; ax.scatter(i + rng_.uniform(-0.18, 0.18, len(v)), v, s=4, color=cols_h[i], alpha=0.8, lw=0)
    ax.errorbar(i, v.mean(), yerr=v.std() / np.sqrt(len(v)), fmt='_', color='k', ms=14, mew=1.2, capsize=0, elinewidth=1.0)
ax.set_xticks(range(3)); ax.set_xticklabels(list(grp), fontsize=6.3); ax.set_ylabel('YAP/TAZ target score'); ax.axhline(0, color='0.75', lw=0.5)
from scipy.stats import mannwhitneyu as _mw
p1 = _mw(grp['F3, F4-like\ncomposition'].YAP_TAZ_targets, grp['F4'].YAP_TAZ_targets).pvalue
ax.plot([1, 1, 2, 2], [2.75, 2.85, 2.85, 2.75], color='k', lw=0.6); ax.text(1.5, 2.9, f'P = {p1:.1g}', ha='center', va='bottom', fontsize=6.3)
ax.set_ylim(-1.1, 3.3); lab(ax, 'G', -0.3, 1.16)
save(fig, 'Fig2_linear_replacement_thresholds')
# ================= Fig 3 — the F3→F4 mechanical switch (programmes, genes, pathways, drivers, positional scale)
PROG = {'matrix_crosslinking': ('Matrix crosslinking', OI['grey']), 'YAP_TAZ_targets': ('YAP/TAZ targets', OI['red']), 'mechanosensing': ('Mechanosensing', OI['orange']), 'TGFb_SMAD': ('TGF-β / SMAD', OI['blue']), 'senescence_SASP': ('Senescence / SASP', OI['purple'])}
R = pd.read_csv(tab('Table_S15_mechanics_scores_thresholds.csv')).set_index('score')
SWg = pd.read_csv(tab('Table_S17a_F4_switch_per_gene.csv'), index_col=0); GSw = pd.read_csv(tab('Table_S17b_GSEA_F4_switch.csv'))
TFw = pd.read_csv(tab('Table_S17c_TF_activity_F4_switch.csv')); BPw = pd.read_csv(tab('Table_S17e_F4_switch_power_by_band.csv'))
fig = plt.figure(figsize=(W2, W2 * 1.02)); gs3 = fig.add_gridspec(3, 3, height_ratios=[1, 1.05, 1], hspace=0.78, wspace=0.62)
# A — programme scores by stage
ax = fig.add_subplot(gs3[0, 0])
for k, (nm, col) in PROG.items():
    g = S.groupby('estadio')[k].agg(['mean', 'sem']).reindex(ST); ax.errorbar(range(6), g['mean'], yerr=g['sem'], fmt='o-', color=col, ms=2.5, lw=0.9, capsize=1.5, elinewidth=0.6, label=nm)
ax.set_xticks(range(6)); ax.set_xticklabels(ST, fontsize=6.5); ax.axhline(0, color='0.75', lw=0.5); ax.set_ylabel('Programme score (marker z)'); ax.set_xlabel('Fibrosis stage')
ax.legend(fontsize=6, loc='upper left', handlelength=1.1, labelspacing=0.25); lab(ax, 'A', -0.3, 1.16)
# B — step at F4, before / after composition
ax = fig.add_subplot(gs3[0, 1]); x = np.arange(len(PROG))
ax.bar(x - 0.2, [R.loc[k, 't_step_F4'] for k in PROG], 0.4, color='0.62', lw=0, label='Unadjusted')
ax.bar(x + 0.2, [R.loc[k, 't_step_F4_comp_adj'] for k in PROG], 0.4, color=OI['orange'], lw=0, label='Composition-adjusted')
ax.axhline(2, color='0.4', ls='--', lw=0.6); ax.set_xticks(x); ax.set_xticklabels(['Crosslinking', 'YAP/TAZ', 'Mechanosensing', 'TGF-β/SMAD', 'SASP'], fontsize=6.5, rotation=35, ha='right')
ax.set_ylabel('t of step at F4'); ax.set_ylim(0, 8.2); ax.legend(fontsize=6, loc='upper right'); lab(ax, 'B', -0.3, 1.16)
# C — YAP/TAZ vs stellate abundance
ax = fig.add_subplot(gs3[0, 2])
for st in ST: g = S[S.estadio == st]; ax.scatter(g.HSC, g.YAP_TAZ_targets, s=5, color=pal[st], alpha=.85, lw=0, label=st)
sub = S[S.estadio != 'F4']; b = np.polyfit(sub.HSC, sub.YAP_TAZ_targets, 1); xs = np.linspace(S.HSC.min(), S.HSC.max(), 50)
ax.plot(xs, np.polyval(b, xs), color='k', lw=0.8, ls='--', label='Fit, Normal–F3')
ax.set_xlabel('Stellate / myofibroblast score'); ax.set_ylabel('YAP/TAZ target score'); ax.legend(fontsize=6, loc='upper left', ncol=2, markerscale=1.3, handletextpad=0.2, columnspacing=0.6); lab(ax, 'C', -0.3, 1.16)

# D — switch genes: step at F4 before vs after composition adjustment
ax = fig.add_subplot(gs3[1, 0]); d = SWg.dropna(subset=['t_step_F4_comp_adj', 't_step_F4']).copy()
mech_cols = ['YAP/TAZ targets', 'Mechanosensing', 'MRTF–SRF targets']; d['mech'] = d[mech_cols].any(axis=1)
ax.scatter(d.t_step_F4, d.t_step_F4_comp_adj, s=1.5, color='0.83', lw=0, rasterized=True)
up = d[d.switch_gene]; dn = d[d.switch_gene_down]
ax.scatter(up.t_step_F4, up.t_step_F4_comp_adj, s=2.5, color=OI['red'], lw=0, rasterized=True, label=f'Up at F4 (n={len(up)})')
ax.scatter(dn.t_step_F4, dn.t_step_F4_comp_adj, s=2.5, color=OI['blue'], lw=0, rasterized=True, label=f'Down at F4 (n={len(dn)})')
hl = d[d.switch_gene & d.mech]; ax.scatter(hl.t_step_F4, hl.t_step_F4_comp_adj, s=13, facecolor='none', edgecolor='k', lw=0.5, zorder=4)
pos = {'CCN1': (-4.6, 7.9), 'CCN2': (-6.6, 6.6), 'AMOTL2': (-4.6, 5.3), 'THBS1': (1.2, 8.1), 'GADD45B': (7.0, 2.9), 'TPM1': (7.0, 1.6), 'TAGLN': (7.0, 0.3)}
texts = []
for _, r in hl.sort_values('t_step_F4_comp_adj', ascending=False).iterrows():
    if r.gene_name in pos:
        tx, ty = pos[r.gene_name]
        ax.annotate(r.gene_name, (r.t_step_F4, r.t_step_F4_comp_adj), xytext=(tx, ty), fontsize=6, style='italic', ha='left' if tx > r.t_step_F4 else 'right', va='center',
                    arrowprops=dict(arrowstyle='-', color='0.45', lw=0.35, shrinkA=0, shrinkB=2))
lim = 8.5; ax.plot([-lim, lim], [-lim, lim], color='0.6', lw=0.5, ls=':'); ax.axhline(0, color='0.75', lw=0.4); ax.axvline(0, color='0.75', lw=0.4)
ax.set_xlim(-lim, 9.8); ax.set_ylim(-lim, 9.0)
ax.set_xlabel('Step at F4, unadjusted (t)'); ax.set_ylabel('Step at F4, composition-adjusted (t)')
ax.legend(fontsize=6, loc='lower right', markerscale=2.2, handletextpad=0.2, borderaxespad=0.2); lab(ax, 'D', -0.3, 1.14)

# E — heatmap: top switch genes by stage (mean z of batch-corrected expression)
ax = fig.add_subplot(gs3[1, 1:]); A_ = pd.read_pickle(inter('expr_adj.pkl'))
topg = SWg[SWg.switch_gene].sort_values('t_step_F4_comp_adj', ascending=False).head(28)
Z = A_.loc[topg.index]; Z = Z.sub(Z.mean(1), axis=0).div(Z.std(1), axis=0)
Hm = pd.DataFrame({st: Z[M.index[M.estadio == st]].mean(axis=1) for st in ST}); Hm.index = topg.gene_name.values
vm = float(np.nanpercentile(np.abs(Hm.values), 97)); im = ax.imshow(Hm.T.values, cmap='RdBu_r', vmin=-vm, vmax=vm, aspect='auto', interpolation='nearest')
ax.set_yticks(range(len(ST))); ax.set_yticklabels(ST); ax.set_xticks(range(len(Hm))); ax.set_xticklabels(Hm.index, rotation=60, ha='right', fontsize=6, style='italic')
mech_all = set(sum([[g for g in v] for v in [SWg[SWg[c]].gene_name for c in mech_cols]], []))
for t_, gname in zip(ax.get_xticklabels(), Hm.index):
    if gname in mech_all: t_.set_color(OI['red']); t_.set_fontweight('bold')
ax.axhline(4.5, color='k', lw=0.8); ax.tick_params(length=0); ax.spines[['left', 'bottom']].set_visible(False)
cb = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.01); cb.set_label('Mean z', fontsize=6.5); cb.ax.tick_params(labelsize=6)
ax.set_title('Top switch genes (red: mechanotransduction sets)', fontsize=6.5, loc='left', pad=3); lab(ax, 'E', -0.07, 1.14)

# F — pathways switching at F4
ax = fig.add_subplot(gs3[2, 0]); keepterms = ['TNFA_SIGNALING_VIA_NFKB', 'HYPOXIA', 'MECH: MRTF–SRF targets', 'MECH: YAP/TAZ targets', 'MECH: Mechanosensing', 'EPITHELIAL_MESENCHYMAL_TRANSITION', 'TGF_BETA_SIGNALING', 'P53_PATHWAY', 'MECH: Senescence/SASP', 'ANGIOGENESIS']
g_ = GSw.set_index('Term').reindex([t for t in keepterms if t in set(GSw.Term)]).iloc[::-1]
nmf = lambda t: t.replace('MECH: ', '').replace('_', ' ').title().replace('Tnfa Signaling Via Nfkb', 'TNFα / NF-κB').replace('Epithelial Mesenchymal Transition', 'EMT').replace('Tgf Beta Signaling', 'TGF-β').replace('P53 Pathway', 'p53').replace('Mrtf–Srf Targets', 'MRTF–SRF targets').replace('Yap/Taz Targets', 'YAP/TAZ targets').replace('Senescence/Sasp', 'Senescence/SASP')
cols_ = [OI['red'] if t.startswith('MECH') else '0.55' for t in g_.index]
ax.barh(range(len(g_)), g_.NES, color=cols_, lw=0, height=0.68); ax.set_yticks(range(len(g_))); ax.set_yticklabels([nmf(t) for t in g_.index], fontsize=6.5)
ax.axvline(0, color='k', lw=0.5); ax.set_xlabel('NES, step at F4 (adjusted)'); ax.set_xlim(0, 3.8)
for i, (v, q) in enumerate(zip(g_.NES, g_['FDR q-val'])): ax.text(v + 0.05, i, '***' if q < 0.001 else ('**' if q < 0.01 else ('*' if q < 0.05 else '')), va='center', fontsize=6)
lab(ax, 'F', -0.62, 1.14)
# G — transcription-factor drivers
ax = fig.add_subplot(gs3[2, 1]); show = ['SMAD4', 'RELA', 'STAT3', 'SMAD3', 'HIF1A', 'TP53', 'NFKB1', 'EGR1', 'KLF4', 'JUN', 'SRF', 'TEAD4', 'TEAD1']
tq = TFw.set_index('tf').reindex([t for t in show if t in set(TFw.tf)]).iloc[::-1]
colt = [OI['orange'] if q < 0.05 else '0.75' for q in tq.q_adj]
ax.barh(range(len(tq)), tq.t_step_F4_comp_adj, color=colt, lw=0, height=0.68); ax.set_yticks(range(len(tq))); ax.set_yticklabels(tq.index, style='italic', fontsize=6.5)
ax.axvline(0, color='k', lw=0.5); ax.set_xlabel('TF activity, step at F4 (adjusted t)')
ax.legend(handles=[Patch(fc=OI['orange'], label='FDR < 0.05'), Patch(fc='0.75', label='n.s.')], fontsize=6, loc='lower right'); lab(ax, 'G', -0.34, 1.14)
# H — positional scale of the switch
ax = fig.add_subplot(gs3[2, 2]); BPw = BPw.iloc[::-1].reset_index(drop=True)
labs_b = [b.replace('(', '').replace(']', '').replace(', ', '–').replace('.0', '') for b in BPw.band]
ax.bar(range(len(BPw)), BPw.z, color=[OI['red'] if z > 2 else (OI['blue'] if z < -2 else '0.7') for z in BPw.z], lw=0, width=0.7)
ax.axhline(0, color='k', lw=0.5); ax.axhline(2, color='0.5', ls=':', lw=0.6); ax.axhline(-2, color='0.5', ls=':', lw=0.6)
ax.set_xticks(range(len(BPw))); ax.set_xticklabels(labs_b, rotation=50, ha='right', fontsize=6); ax.set_xlabel('Period (genes)'); ax.set_ylabel('Power of switch along genome\n(z vs gene-order null)')
lab(ax, 'H', -0.32, 1.14)
save(fig, 'Fig3_mechanical_switch')
# ======================================================= FIGURE 4
A = pd.read_csv(tab('Table_S5d_acf_with_without_composition.csv')); D = pd.read_csv(tab('Table_S5g_concordance_by_distance.csv'), index_col=0)
SC = pd.read_csv(tab('Table_S8c_within_type_spatial_autocorrelation.csv')); CO = pd.read_csv(tab('Table_S8b_within_type_neighbour_coexpression.csv'))
EQ = pd.read_csv(tab('Table_S13e_shared_GTEx_eQTL_coupled_vs_uncoupled.csv'))
en = {'Fagocito_mononuclear': 'Mono. phagocyte', 'Endotelio': 'Endothelium', 'Mesenquima_HSC': 'Mesenchyme/HSC', 'Colangiocito': 'Cholangiocyte',
      'Plasma': 'Plasma', 'T_NK': 'T/NK', 'B': 'B', 'pDC': 'pDC', 'Hepatocito': 'Hepatocyte'}
SN_AC = pd.read_csv(tab('Table_S18f_snRNA_neighbour_coupling.csv')); SN_CO = pd.read_csv(tab('Table_S18h_snRNA_neighbour_coexpression.csv'))
fig = plt.figure(figsize=(W2, W2 * 0.66)); g4 = fig.add_gridspec(2, 6, height_ratios=[1, 1.05], hspace=0.78, wspace=2.2)
axs = [fig.add_subplot(g4[0, 0:3]), fig.add_subplot(g4[0, 3:6]), fig.add_subplot(g4[1, 0:2]), fig.add_subplot(g4[1, 2:4]), fig.add_subplot(g4[1, 4:6])]
ax = axs[0]; ax.axhspan(-A['null_p97.5'].max(), A['null_p97.5'].max(), color='0.9', lw=0, label='95% permutation null')
ax.plot(A.lag, A.acf_unadjusted, 'o-', color='k', ms=3, lw=0.9, label='Stage effect per gene')
ax.plot(A.lag, A.acf_composition_adjusted, 's--', color=OI['orange'], ms=3, lw=0.9, label='Adjusted for cell composition')
ax.set_xscale('log'); ax.set_xlabel('Distance (genes)'); ax.set_ylabel('Spatial autocorrelation'); ax.legend(); lab(ax, 'A', -0.17, 1.16)
ax = axs[1]; xx = range(len(D))
ax.plot(xx, D.r_t, 'o-', color='k', ms=3, lw=0.9, label='Concordance of stage effect')
ax.plot(xx, D.r_t_comp_adj, 'o--', color=OI['orange'], ms=3, lw=0.9, label='…adjusted for composition')
ax.plot(xx, D.coexpr, 's-', color=OI['blue'], ms=3, lw=0.9, label='Co-expression across biopsies')
ax.set_xticks(list(xx)); ax.set_xticklabels(['overlap', '<1', '1–10', '10–50', '50–100', '100–500', '>500'], rotation=45, ha='right')
ax.set_xlabel('Intergenic distance (kb)'); ax.set_ylabel('Correlation'); ax.axhline(0, color='0.6', lw=0.5); ax.set_ylim(-0.08, 0.52); ax.legend(loc='upper right'); lab(ax, 'B', -0.17, 1.16)
# C — scRNA-seq (dissociated cells): cirrhosis vs healthy, non-parenchymal populations
ax = axs[2]; SCs = SC[SC.type != 'Hepatocito'].sort_values('acf_lag1'); y = np.arange(len(SCs))
ax.barh(y, SCs.acf_lag1, color=OI['red'], lw=0, height=0.66)
ax.errorbar([0] * len(SCs), y, xerr=SCs['null_p97.5'], fmt='none', color='k', capsize=2, lw=0.6)
ax.set_yticks(y); ax.set_yticklabels([en[t] for t in SCs.type]); ax.set_xlabel('Neighbour autocorrelation\n(cirrhosis vs healthy)'); ax.set_title('scRNA-seq, 10 donors', fontsize=7, pad=3)
lab(ax, 'C', -0.62, 1.2)
# D — snRNA-seq (nuclei): stage effect across 47 donors, all populations including hepatocytes
ax = axs[3]; pops = ['Hepatocyte', 'Endothelium', 'HSC', 'Cholangiocyte', 'B/Plasma', 'T/NK', 'Macrophage']
sa = SN_AC.set_index('population').reindex(pops).iloc[::-1]; y = np.arange(len(sa))
ax.barh(y, sa.acf_lag1, color=[OI['green'] if p == 'Hepatocyte' else OI['red'] for p in sa.index], lw=0, height=0.66)
ax.errorbar([0] * len(sa), y, xerr=sa['null_p97.5'], fmt='none', color='k', capsize=2, lw=0.6, label='95% null')
ax.set_yticks(y); ax.set_yticklabels(sa.index); ax.set_xlabel('Neighbour autocorrelation\n(stage effect, F0–F4)'); ax.set_title('snRNA-seq, 47 donors', fontsize=7, pad=3); ax.legend(loc='lower right')
lab(ax, 'D', -0.55, 1.2)
# E — snRNA-seq: co-expression of neighbouring genes across single nuclei
ax = axs[4]; sc_ = SN_CO[SN_CO.donors >= 8].sort_values('r_neighbours'); y = np.arange(len(sc_))
ax.barh(y + 0.19, sc_.r_neighbours, 0.38, color=[OI['green'] if t == 'Hepatocyte' else OI['blue'] for t in sc_.type], lw=0, label='Neighbouring genes')
ax.barh(y - 0.19, sc_.r_random, 0.38, color='0.65', lw=0, label='Random pairs')
ax.set_yticks(y); ax.set_yticklabels([f'{t} ({int(d)})' for t, d in zip(sc_.type, sc_.donors)]); ax.set_xlabel('Co-expression across nuclei\n(snRNA-seq; donors in brackets)')
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=6, columnspacing=0.8, handlelength=1.2); lab(ax, 'E', -0.62, 1.2)
r = EQ.iloc[1]
fig.text(0.99, -0.02, f'Coupled neighbours share a same-direction liver eQTL more often than uncoupled pairs at matched distance: {100*r.frac_coupled:.0f}% vs {100*r.frac_uncoupled:.0f}%, odds ratio {r.MH_OR_distance_stratified:.2f} (95% CI {r.CI_low:.2f}–{r.CI_high:.2f}), P = {r.MH_p:.3f}',
         ha='right', va='top', fontsize=6)
save(fig, 'Fig4_neighbourhood_coupling')
# ================= Fig 6 (GWAS)
Gw = pd.read_csv(tab('Table_S16c_GWAS_catalog_enrichment_by_trait.csv')); order = ['liver fibrosis/cirrhosis', 'MASLD/NAFLD', 'ALT (liver enzyme)', 'height (control)', 'educational attainment (control)']; Gw = Gw.set_index('trait_set').loc[order]
lbl = ['Fibrosis/\ncirrhosis', 'MASLD', 'ALT', 'Height', 'Education']; cols = [OI['red'], OI['orange'], OI['sky'], '0.6', '0.6']
fig, axs = plt.subplots(1, 2, figsize=(W2 * 0.85, W2 * 0.36)); plt.subplots_adjust(wspace=0.45)
ax = axs[0]; ax.bar(range(5), 100 * Gw.frac_in_coupled_pair, color=cols, lw=0, width=0.65); ax.axhline(100 * Gw.bg_coupled.iloc[0], color='k', ls='--', lw=0.7, label='All genes')
for i, (f, p) in enumerate(zip(Gw.frac_in_coupled_pair, Gw.P_coupled)): ax.text(i, 100 * f + 0.8, f'P={p:.3f}' if p < 0.1 else 'n.s.', ha='center', fontsize=6)
ax.set_xticks(range(5)); ax.set_xticklabels(lbl, fontsize=6.5); ax.set_xlabel('GWAS trait (grey = non-liver controls)'); ax.set_ylabel('% of locus genes in a\ncis-coupled neighbour pair'); ax.set_ylim(0, 40); ax.legend(loc='upper right'); lab(ax, 'A')
ax = axs[1]; ax.bar(range(5), Gw.acf_windows, color=cols, lw=0, width=0.65); ax.errorbar(range(5), Gw.acf_random, yerr=2 * Gw.acf_random_sd, fmt='o', color='k', ms=3, capsize=3, lw=0.6, label='Random windows ± 2 s.d.')
for i, p in enumerate(Gw.P_acf): ax.text(i, Gw.acf_windows.iloc[i] + 0.008, f'P={p:.3f}' if p < 0.1 else 'n.s.', ha='center', fontsize=6)
ax.set_xticks(range(5)); ax.set_xticklabels(lbl, fontsize=6.5); ax.set_xlabel('GWAS trait'); ax.set_ylabel('Spatial autocorrelation of stage\neffect, ±5 genes around loci'); ax.legend(loc='upper right'); lab(ax, 'B')
save(fig, 'SuppFig_GWAS_loci_neighbourhoods')
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
tot = sum(N.values()); hubsize = int(0.115 * tot)
nfkb = set(Dor[Dor.tf.isin(['NFKB1', 'RELA'])].target); cand = set(Cn.gene_name)
circos = Circos({**{c: N[c] for c in CHR}, **{short[h]: hubsize for h in hubs}}, space=1.8, start=-28, end=314)
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


for rr, lbl in [(92, '1'), (77, '2'), (67, '3')]:
    circos.text(lbl, r=rr, deg=323, size=7.5, fontweight='bold', color='0.15', va='center', ha='center')
circos.text('Rings', r=108, deg=323, size=6.5, color='0.35', va='center', ha='center')
fig = plt.figure(figsize=(W2, W2 * 0.93)); axc = fig.add_axes([0.11, 0.165, 0.76, 0.82], polar=True); circos.plotfig(ax=axc)
from matplotlib.legend_handler import HandlerTuple
hdc = [(Patch(fc=OI['red']), Patch(fc=OI['blue'])), (Patch(fc='#b05a2c'), Patch(fc='#2c5f8c')),
       (Patch(fc=OI['red'], alpha=0.8), Patch(fc=OI['blue'], alpha=0.8)), Line2D([], [], color=OI['green'], lw=1.6),
       Line2D([], [], marker='o', ls='', color=OI['orange'], ms=3.5), Patch(fc=plt.cm.tab10(0), alpha=0.35),
       Line2D([], [], color='0.45', lw=1.0)]
lbc = ['Ring 1 — stage effect per gene (up / down)', 'Ring 2 — same, adjusted for cell composition',
       'Ring 3 — DE neighbourhood, ≥3 genes (up / down)', 'Ring 3 — lineage-intrinsic target gene',
       'Ring 3 — NF-κB target gene', 'Hallmark pathway sector (GSEA NES)', 'Link: neighbourhood → pathway']
fig.legend(hdc, lbc, handler_map={tuple: HandlerTuple(ndivide=None, pad=0.15)}, loc='lower center', bbox_to_anchor=(0.5, 0.005),
           ncol=4, fontsize=6.5, handlelength=2.0, columnspacing=1.0, labelspacing=0.6, handletextpad=0.5, frameon=False)
save(fig, 'Fig5_genome_circos')

# ---- Fig 6: two layers (A, B) on top, lineage-intrinsic targets (C) full width below
fig = plt.figure(figsize=(W2, W2 * 1.27)); gs = fig.add_gridspec(4, 2, height_ratios=[0.92, 0.95, 0.80, 0.92], hspace=0.60, wspace=0.62)
# A — GSEA two layers
ax = fig.add_subplot(gs[0, 0])
paths = ['TNFA_SIGNALING_VIA_NFKB', 'EPITHELIAL_MESENCHYMAL_TRANSITION', 'INFLAMMATORY_RESPONSE', 'ANGIOGENESIS', 'P53_PATHWAY', 'ADIPOGENESIS', 'FATTY_ACID_METABOLISM', 'OXIDATIVE_PHOSPHORYLATION', 'BILE_ACID_METABOLISM']
paths = [p for p in paths if p in G1.index and p in G2.index]
nm = lambda p: p.replace('_', ' ').title().replace('Tnfa Signaling Via Nfkb', 'TNFα / NF-κB').replace('Epithelial Mesenchymal Transition', 'EMT').replace('P53', 'p53').replace('Oxidative Phosphorylation', 'OXPHOS')
y = np.arange(len(paths))
for i, p in enumerate(paths): ax.plot([G1.loc[p, 'NES'], G2.loc[p, 'NES']], [i, i], color='0.75', lw=1)
ax.scatter([G1.loc[p, 'NES'] for p in paths], y, s=16, color='0.35', zorder=3, label='Raw stage statistic')
ax.scatter([G2.loc[p, 'NES'] for p in paths], y, s=16, color=OI['orange'], zorder=3, label='Composition-adjusted')
ax.axvline(0, color='k', lw=0.5); ax.set_yticks(y); ax.set_yticklabels([nm(p) for p in paths], fontsize=6); ax.invert_yaxis()
ax.set_xlabel('GSEA normalised enrichment score'); ax.set_xlim(-3.6, 4.6); ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.01), ncol=2, fontsize=6.5, handletextpad=0.4, columnspacing=1.6); lab(ax, 'A', -0.40, 1.20)
# C — TF activity
ax = fig.add_subplot(gs[0, 1]); tfs = [t for t in ['HNF4A', 'HNF1A', 'FOXA1', 'PPARA', 'NFKB1', 'RELA', 'TP53', 'SMAD3', 'SMAD4', 'KLF4', 'FOXM1'] if t in TF.index]
y = np.arange(len(tfs))
for i, t in enumerate(tfs): ax.plot([TF.loc[t, 't_stage'], TF.loc[t, 't_stage_comp_adjusted']], [i, i], color='0.75', lw=1)
ax.scatter([TF.loc[t, 't_stage'] for t in tfs], y, s=16, color='0.35', zorder=3)
ax.scatter([TF.loc[t, 't_stage_comp_adjusted'] for t in tfs], y, s=16, color=OI['orange'], zorder=3)
ax.axvline(0, color='k', lw=0.5); ax.set_yticks(y); ax.set_yticklabels(tfs, style='italic', fontsize=6); ax.invert_yaxis()
ax.axhline(3.5, color='0.8', lw=0.5, ls=':'); ax.set_xlabel('t of stage on inferred TF activity')
ax.set_xlim(-9.5, 21.5)
for tick, t in zip(ax.get_yticklabels(), tfs): tick.set_color(OI['green'] if t in ('HNF4A', 'HNF1A', 'FOXA1', 'PPARA') else OI['red'])
hd3 = [Line2D([], [], color=OI['green'], lw=2.4, label='Hepatocyte identity'), Line2D([], [], color=OI['red'], lw=2.4, label='Intrinsic stress response')]
ax.legend(handles=hd3, loc='lower center', bbox_to_anchor=(0.5, 1.01), ncol=2, fontsize=6.5, handlelength=1.1, columnspacing=1.4)
lab(ax, 'B', -0.20, 1.20)
# ---- target timing data
TGt = pd.read_csv(tab('Table_S17h_target_timing.csv'), index_col=0)
LIN = {'Mesenquima_HSC': ('Stellate', OI['red']), 'Fagocito_mononuclear': ('Macrophage', OI['orange']), 'Endotelio': ('Endothelium', OI['blue']), 'Colangiocito': ('Cholangiocyte', OI['purple'])}
DEV = {'clinical, liver': ('Clinical, liver', OI['red']), 'clinical, inflammasome': ('Clinical, other', OI['orange']), 'clinical, other fibrosis': ('Clinical, other', OI['orange']),
       'clinical, other': ('Clinical, other', OI['orange']), 'approved, other': ('Clinical, other', OI['orange']), 'preclinical': ('Preclinical', OI['sky']), 'tool compounds': ('Preclinical', OI['sky']),
       'biomarker': ('Biomarker', OI['green']), 'none': ('No modulator', '0.62'), 'not curated': ('Not curated', '0.88')}

# C — when does each target rise? trend versus switch
ax = fig.add_subplot(gs[1, 0])
ax.axvspan(2, 9, ymin=0, ymax=1, color=OI['sky'], alpha=0.07, lw=0); ax.axhline(0, color='0.75', lw=0.4); ax.axvline(0, color='0.75', lw=0.4)
for lk, (ln, lc) in LIN.items():
    d_ = TGt[TGt.lineage == lk]; ax.scatter(d_.t_trend_comp_adj, d_.t_step_F4_comp_adj, s=12, color=lc, lw=0.3, edgecolor='white', alpha=0.9, label=ln, zorder=3)
HEP = pd.read_csv(tab('Table_S9c_hepatocyte_intrinsic_genes.csv'), index_col=0); SWh = pd.read_csv(tab('Table_S17a_F4_switch_per_gene.csv'), index_col=0)
HEP = HEP.join(SWh[['t_trend_comp_adj', 't_step_F4_comp_adj']], how='left')
ax.scatter(HEP.t_trend_comp_adj, HEP.t_step_F4_comp_adj, s=13, marker='D', color=OI['green'], lw=0.3, edgecolor='white', alpha=0.9, label='Hepatocyte (bulk-defined)', zorder=3)
sw_ = TGt[TGt.switch_gene.fillna(False)]; ax.scatter(sw_.t_trend_comp_adj, sw_.t_step_F4_comp_adj, s=26, facecolor='none', edgecolor='k', lw=0.6, zorder=4)
lbl6 = {'THBS2': (7.6, 3.0), 'TREM2': (8.3, -3.9), 'IL32': (7.6, -1.4), 'COL1A1': (4.6, -5.6), 'LGALS3': (1.3, -2.6), 'CCN2': (-3.2, 6.6), 'CCN1': (0.6, 7.4), 'TIMP1': (4.9, 6.2), 'IER3': (-3.2, 4.4), 'VWF': (4.6, 4.4), 'MAT1A': (-6.0, 1.6), 'FGF21': (1.9, -5.1)}
for gname, (tx, ty) in lbl6.items():
    r_ = TGt[TGt.gene_name == gname]
    if not len(r_): r_ = HEP[HEP.gene_name == gname]
    if len(r_): ax.annotate(gname, (r_.t_trend_comp_adj.iloc[0], r_.t_step_F4_comp_adj.iloc[0]), xytext=(tx, ty), fontsize=6, style='italic', ha='center', va='center',
                            arrowprops=dict(arrowstyle='-', lw=0.35, color='0.45', shrinkA=0, shrinkB=2), zorder=5)
ax.set_xlim(-8.6, 9.0); ax.set_ylim(-6.2, 8.0); ax.axvspan(-8.6, -2, color=OI['green'], alpha=0.05, lw=0)
ax.text(8.8, 7.7, 'rising,\nlinear', ha='right', va='top', fontsize=6.3, color=OI['blue']); ax.text(-8.4, 7.7, 'falling,\nlinear', ha='left', va='top', fontsize=6.3, color=OI['green']); ax.set_xlabel('Trend with stage, adjusted (t)'); ax.set_ylabel('Step at F4, adjusted (t)')
hc, lc_ = ax.get_legend_handles_labels(); hc.append(Line2D([], [], marker='o', ls='', markerfacecolor='none', markeredgecolor='k', ms=4.5)); lc_.append('Switch at F4')
ax.legend(hc, lc_, fontsize=6, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=3, handletextpad=0.15, columnspacing=0.7, labelspacing=0.25, borderaxespad=0.2); lab(ax, 'C', -0.24, 1.28)

# D — trajectories of representative early and switch targets
ax = fig.add_subplot(gs[1, 1]); A6 = pd.read_pickle(inter('expr_adj.pkl')); nm6 = keep.set_index('gene_id').gene_name.astype(str)
traj = [('THBS2', OI['red'], '-'), ('TREM2', OI['orange'], '-'), ('IL32', OI['purple'], '-'), ('CCN2', OI['blue'], '--'), ('CCN1', OI['sky'], '--'), ('TIMP1', OI['green'], '--')]
for gname, col, ls_ in traj:
    gid = nm6.index[nm6 == gname]
    if not len(gid): continue
    v = A6.loc[gid[0]]; z = (v - v.mean()) / v.std(); g = z.groupby(M.estadio).agg(['mean', 'sem']).reindex(ST)
    ax.errorbar(range(6), g['mean'], yerr=g['sem'], color=col, ls=ls_, lw=1.0, marker='o', ms=2.4, capsize=1.2, elinewidth=0.5, label=gname)
ax.set_xticks(range(6)); ax.set_xticklabels(ST, fontsize=6.5); ax.axhline(0, color='0.8', lw=0.4); ax.set_xlabel('Fibrosis stage'); ax.set_ylabel('Expression (z)')
h_, l_ = ax.get_legend_handles_labels()
lg1 = ax.legend(h_[:3], l_[:3], title='Early, linear', title_fontsize=6, fontsize=6, loc='upper left', handlelength=1.6, labelspacing=0.25, borderaxespad=0.3); ax.add_artist(lg1)
ax.legend(h_[3:], l_[3:], title='Switch at F4', title_fontsize=6, fontsize=6, loc='upper left', bbox_to_anchor=(0.30, 1.0), handlelength=1.6, labelspacing=0.25, borderaxespad=0.3)
for t_ in ax.get_legend().get_texts() + lg1.get_texts(): t_.set_fontstyle('italic')
ax.set_ylim(-1.5, 2.3); lab(ax, 'D', -0.2, 1.16)

# E — target landscape: within-lineage induction, timing and pharmacology
ax = fig.add_subplot(gs[2, :]); _p = ax.get_position(); ax.set_position([_p.x0, _p.y0, _p.width * 0.83, _p.height])
early = TGt[TGt.timing == 'early, linear'].sort_values('score', ascending=False).head(15)
swt = TGt[TGt.timing == 'switch at F4'].sort_values('t_step_F4_comp_adj', ascending=False).head(9)
sel = pd.concat([early, swt]); rows_ = ['Mesenquima_HSC', 'Fagocito_mononuclear', 'Endotelio', 'Colangiocito', 'Hepatocito_sn']
SNt = pd.read_csv(tab('Table_S18b_snRNA_within_population_stage_t.csv'), index_col=0)['Hepatocyte']; sel = sel.assign(t_sc_Hepatocito_sn=SNt.reindex(sel.gene_name.values).values)
for j, (_, r) in enumerate(sel.iterrows()):
    for i, rk in enumerate(rows_):
        t_ = r[f't_sc_{rk}']
        if pd.notna(t_) and t_ > 0:
            ax.scatter(j, i, s=min(t_, 5) ** 2 * 5.5, color=plt.cm.Reds(0.3 + 0.7 * min(t_, 5) / 5), edgecolor='0.3', lw=0.25, zorder=3)
    dv = DEV.get(r.development_stage, ('Not curated', '0.88'))
    ax.add_patch(plt.Rectangle((j - 0.42, 5.05), 0.84, 0.5, color=dv[1], lw=0, clip_on=False))
ax.axvline(len(early) - 0.5, color='k', lw=0.6, ls=':')
ax.text((len(early) - 1) / 2, -1.25, 'Early, linear targets', ha='center', fontsize=6.5, color=OI['blue'], fontweight='bold')
ax.text(len(early) + (len(swt) - 1) / 2, -1.25, 'Switch-at-F4 targets', ha='center', fontsize=6.5, fontweight='bold')
ax.set_xlim(-0.6, len(sel) - 0.4); ax.set_ylim(5.75, -0.6); ax.set_yticks(list(range(5)) + [5.3]); ax.set_yticklabels([LIN[r][0] for r in rows_[:4]] + ['Hepatocyte*', 'Drug status'], fontsize=6.5)
ax.axhline(3.5, color='0.8', lw=0.5)
ax.set_xticks(range(len(sel))); ax.set_xticklabels(sel.gene_name, rotation=60, ha='right', fontsize=6.3, style='italic')
ax.spines[['left', 'bottom']].set_visible(False); ax.tick_params(length=0); ax.grid(False)
sz = [Line2D([], [], marker='o', ls='', color=plt.cm.Reds(0.3 + 0.7 * v / 5), markeredgecolor='0.3', markeredgewidth=0.25, ms=np.sqrt(v ** 2 * 5.5), label=f't = {v}') for v in (2, 3, 4)]
dvh = [Patch(fc=c, label=l) for l, c in dict(v for v in DEV.values()).items()]
lg = ax.legend(handles=sz, title='Within-lineage\ninduction (t)\nscRNA-seq;\n*snRNA-seq, F0–F4', title_fontsize=6, fontsize=6, loc='upper left', bbox_to_anchor=(1.02, 1.12), labelspacing=0.75, borderaxespad=0)
ax.add_artist(lg); ax.legend(handles=dvh, title='Drug status', title_fontsize=6, fontsize=6, loc='lower left', bbox_to_anchor=(1.02, -0.40), labelspacing=0.25, borderaxespad=0, handlelength=1.0)
lab(ax, 'E', -0.085, 1.22)
# F — receptors of hepatocyte-directed MASH drugs: bulk change and the part left after removing hepatocyte loss
ax = fig.add_subplot(gs[3, 0]); RF = pd.read_csv(tab('Table_S9d_hepatocyte_drug_targets_reference.csv'), index_col=0)
drug = {'THRB': 'resmetirom', 'PPARA': 'lanifibranor', 'PPARD': 'lanifibranor', 'PPARG': 'lanifibranor', 'NR1H4': 'obeticholic acid', 'FASN': 'denifanstat', 'ACACA': 'firsocostat',
        'ACACB': 'firsocostat', 'SCD': 'aramchol', 'HSD17B13': 'HSD17B13 siRNA', 'PNPLA3': 'PNPLA3 ASO', 'DGAT2': 'ION224', 'FGFR1': 'FGF21 analogues', 'KLB': 'FGF21 analogues'}
RF = RF[RF.gene_name.isin(drug)].sort_values('t_bulk_stage'); yy = np.arange(len(RF))
for i, (_, r) in enumerate(RF.iterrows()): ax.plot([r.t_bulk_stage, r.t_bulk_comp_adjusted], [i, i], color='0.78', lw=1.0, zorder=1)
ax.scatter(RF.t_bulk_stage, yy, s=15, color='0.35', zorder=3, label='Bulk stage effect')
ax.scatter(RF.t_bulk_comp_adjusted, yy, s=15, color=OI['green'], zorder=3, label='After removing composition')
RFs = pd.read_csv(tab('Table_S18e_snRNA_hepatocyte_drug_receptors.csv'), index_col=0).set_index('gene_name').t_snRNA_hepatocyte
ax.scatter(RFs.reindex(RF.gene_name).values, yy, s=16, marker='D', facecolor='white', edgecolor=OI['green'], lw=0.8, zorder=4, label='Hepatocyte nuclei (snRNA-seq)')
ax.axvline(0, color='k', lw=0.5); ax.axvline(-2, color='0.6', lw=0.5, ls=':'); ax.set_yticks(yy)
ax.set_yticklabels([f'{g} · {drug[g]}' for g in RF.gene_name], fontsize=6.2)
for t_ in ax.get_yticklabels():
    g = t_.get_text().split(' · ')[0]
    if g == 'THRB': t_.set_fontweight('bold')
for i, (_, r) in enumerate(RF.iterrows()):
    if pd.notna(r.get('t_step_F4_comp_adj')) and r.t_step_F4_comp_adj < -2: ax.text(np.nanmax([r.t_bulk_comp_adjusted, RFs.get(r.gene_name, np.nan)]) + 0.7, i, '▼F4', fontsize=5.6, va='center', color=OI['red'])
ax.set_xlabel('t of stage on expression'); ax.set_xlim(-10.5, 8.5)
ax.legend(fontsize=6, loc='lower center', bbox_to_anchor=(0.45, 1.0), ncol=2, handletextpad=0.3, columnspacing=1.0); lab(ax, 'F', -0.62, 1.28)
# G — hepatocyte-intrinsic genes (bulk-defined): strongest falling and rising
ax = fig.add_subplot(gs[3, 1]); HP = pd.read_csv(tab('Table_S9c_hepatocyte_intrinsic_genes.csv'), index_col=0)
hp = pd.concat([HP[HP.direction == 'down'].sort_values('t_bulk_comp_adjusted').head(9), HP[HP.direction == 'up'].sort_values('t_bulk_comp_adjusted', ascending=False).head(9).iloc[::-1]])
yy = np.arange(len(hp)); ax.barh(yy, hp.t_bulk_comp_adjusted, color=[OI['green'] if v < 0 else OI['red'] for v in hp.t_bulk_comp_adjusted], lw=0, height=0.68)
ax.set_yticks(yy); ax.set_yticklabels(hp.gene_name, style='italic', fontsize=6.2); ax.axvline(0, color='k', lw=0.5); ax.invert_yaxis()
HS = pd.read_csv(tab('Table_S18d_snRNA_hepatocyte_class.csv'), index_col=0).set_index('gene_name').t_snRNA_hepatocyte
ax.scatter(HS.reindex(hp.gene_name).values, yy, s=14, marker='D', facecolor='white', edgecolor='k', lw=0.7, zorder=4, label='Hepatocyte nuclei (snRNA-seq)')
ax.legend(fontsize=6, loc='lower center', bbox_to_anchor=(0.5, 1.0), handletextpad=0.3)
ax.set_xlabel('Stage effect (t): bars, bulk after removing composition')
ax.text(0.02, 0.03, f'{(HP.direction == "down").sum()} falling, {(HP.direction == "up").sum()} rising\nhepatocyte-specific genes', transform=ax.transAxes, ha='left', va='bottom', fontsize=6.2)
lab(ax, 'G', -0.2, 1.28)
save(fig, 'Fig6_two_layers_and_targets')
T1 = Cn.head(25)[['gene_name', 'cell_type_up_sc', 't_bulk_stage', 't_bulk_comp_adjusted', 't_sc_max', 'coupled_neighbour', 'target_class', 'drug_landscape', 'development_stage']].copy()
T1.columns = ['Gene', 'Cell type (scRNA-seq)', 'Bulk t (stage)', 'Bulk t (composition-adjusted)', 'Max within-cell-type t', 'Cis-coupled neighbour', 'Target class', 'Agents / status', 'Development stage']
T1['Cell type (scRNA-seq)'] = T1['Cell type (scRNA-seq)'].str.replace('Mesenquima_HSC', 'stellate').str.replace('Fagocito_mononuclear', 'macrophage').str.replace('Endotelio', 'endothelium').str.replace('Colangiocito', 'cholangiocyte')
T1.round(1).to_csv(os.path.join(FO, 'Table1_targets.csv'), index=False)
print('JHEP figure set written to', FO)
