"""17 — Supplementary Figures S1–S9 for the Journal of Hepatology submission (numbered in order of citation).
New here: S1 (positional spectral method and controls), S3 (stage effects on the spectrum, gene-order control),
S4 (orientation, distance and TADs), S5 (shared liver eQTL variants). S2, S6, S7, S8 and S9 are produced by scripts
13, 11, 16, 10 and 18 and copied here under their supplementary numbers.
Outputs: results/figures/jhep/supplementary/SuppFig_S<n>_<name>.{pdf,png}
"""
import os, shutil, numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, scienceplots, seaborn as sns
from scipy import stats
from common import *

plt.style.use(['science', 'nature', 'no-latex'])
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial', 'Helvetica'], 'font.size': 7, 'axes.labelsize': 7,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 6.5, 'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
                     'xtick.top': False, 'ytick.right': False, 'xtick.minor.visible': False, 'ytick.minor.visible': False,
                     'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False, 'pdf.fonttype': 42})
W2 = 183 / 25.4; SO = os.path.join(FIG, 'jhep', 'supplementary'); os.makedirs(SO, exist_ok=True)
ST = ['Normal', 'F0', 'F1', 'F2', 'F3', 'F4']; pal = dict(zip(ST, sns.color_palette('viridis', 6)))
OI = {'blue': '#0072B2', 'orange': '#E69F00', 'green': '#009E73', 'red': '#D55E00', 'purple': '#CC79A7', 'sky': '#56B4E9', 'grey': '#7F7F7F'}
def lab(ax, s, dx=-0.24, dy=1.16): ax.text(dx, dy, s, transform=ax.transAxes, fontsize=8, fontweight='bold', va='top')
def save(fig, name):
    fig.savefig(os.path.join(SO, name + '.pdf'), bbox_inches='tight'); fig.savefig(os.path.join(SO, name + '.png'), dpi=300, bbox_inches='tight'); plt.close(fig)
blab = lambda cs: [str(c).replace('(', '').replace(']', '').replace(', ', '–').replace('.0', '') for c in cs]

X, keep = pd.read_pickle(inter('expr.pkl')); A = pd.read_pickle(inter('expr_adj.pkl')); M = pd.read_pickle(inter('meta.pkl'))
Wz, Wl = pd.read_pickle(inter('W_adj.pkl')); Wmt = pd.read_pickle(inter('W_mt.pkl')); grid = load_grid(); N = chrom_lengths(grid)
U = pd.read_csv(tab('Table_S2a_consensus_spectrum_universal_peaks.csv')); U['chr'] = U.chr.astype(str)

# ======================= S1 — the positional spectral method and its controls
fig = plt.figure(figsize=(W2, W2 * 0.62)); gs = fig.add_gridspec(2, 3, hspace=0.75, wspace=0.5)
c0 = '7'; smp = M.index[M.estadio == 'F2'][0]; g = keep[keep.chr == c0].sort_values('grid_index'); n = N[c0]
v = A.loc[g.gene_id, smp].values; S = np.zeros(n); S[g.grid_index.values - 1] = v - v.mean()
ax = fig.add_subplot(gs[0, 0]); ax.vlines(np.arange(1, n + 1), 0, S, color='0.35', lw=0.3); ax.axhline(0, color='k', lw=0.4)
ax.set_xlabel(f'Gene position on chromosome {c0} (grid index)'); ax.set_ylabel('Centred log expression'); ax.text(0.02, 0.97, f'one F2 biopsy; {len(g)} genes, {n} slots', transform=ax.transAxes, va='top', fontsize=6); lab(ax, 'A')
P = np.abs(np.fft.rfft(S)[1:n // 2 + 1]) ** 2; k = np.arange(1, len(P) + 1); lx = np.log10(k / n); b = np.polyfit(lx, np.log10(P + 1e-12), 1)
ax = fig.add_subplot(gs[0, 1]); ax.loglog(k / n, P, color='0.55', lw=0.4); ax.loglog(k / n, 10 ** (b[0] * lx + b[1]), color=OI['red'], lw=1.0, label=f'1/f background fit (slope {b[0]:.2f})')
ax.set_xlabel('Spatial frequency (cycles per gene slot)'); ax.set_ylabel('Periodogram power'); ax.legend(loc='lower left'); lab(ax, 'B')
w = P / 10 ** (b[0] * lx + b[1]); w = w / w.mean()
ax = fig.add_subplot(gs[0, 2]); ax.plot(k, w, color='0.45', lw=0.4); ax.axhline(3, color=OI['red'], ls='--', lw=0.6, label='3 × background')
up = U[(U.chr == c0) & U.universal].k.values; ax.scatter(up, w[up - 1], s=6, color=OI['orange'], zorder=4, label='Universal peaks (all biopsies)')
ax.set_yscale('log'); ax.set_ylim(1e-3, 60); ax.set_xlabel('Frequency index k'); ax.set_ylabel('Whitened power'); ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=6, columnspacing=0.8); lab(ax, 'C', -0.24, 1.26)
# D null calibration: whitened power of non-universal frequencies against the exponential law of a periodogram ordinate
nonu = ~U.set_index(['chr', 'k']).universal.reindex(Wl.index).fillna(False).values; rng = np.random.default_rng(0)
vals = Wl.values[nonu][:, rng.choice(Wl.shape[1], 60, replace=False)].ravel(); vals = vals[np.isfinite(vals)]
ax = fig.add_subplot(gs[1, 0]); bins = np.linspace(0, 10, 51); ax.hist(vals, bins=bins, density=True, color='0.6', lw=0, label='Non-universal frequencies')
xs = np.linspace(0, 10, 200); ax.plot(xs, np.exp(-xs), color=OI['red'], lw=1.0, label='Exp(1)')
ax.set_yscale('log'); ax.set_xlabel('Whitened power'); ax.set_ylabel('Density'); ax.legend(loc='upper right', fontsize=6)
ax.text(0.03, 0.05, f'P(w > 3): observed {np.mean(vals > 3):.3f}, Exp(1) {np.exp(-3):.3f}', transform=ax.transAxes, fontsize=6); lab(ax, 'D')
# E grid-mask control
ax = fig.add_subplot(gs[1, 1]); mk = (U.w_mask > 0) & np.isfinite(U.w_mask)
ax.scatter(np.log10(U.w_mask[mk]), np.log10(U.mean_power_over_null[mk]), s=1, color='0.6', lw=0, rasterized=True)
ax.scatter(np.log10(U.w_mask[mk & U.universal]), np.log10(U.mean_power_over_null[mk & U.universal]), s=4, color=OI['orange'], lw=0)
r_m = np.corrcoef(np.log(U.w_mask[mk]), np.log(U.mean_power_over_null[mk]))[0, 1]
ax.set_xlabel('log$_{10}$ power of the empty-slot mask'); ax.set_ylabel('log$_{10}$ consensus power'); ax.text(0.03, 0.95, f'r = {r_m:.3f}', transform=ax.transAxes, va='top', fontsize=6.5); lab(ax, 'E')
# F replication of the universal peaks in each cohort separately
ax = fig.add_subplot(gs[1, 2]); cohs = ['GSE130970', 'GSE135251', 'GSE162694']
data = [U[U.universal][f'frac_w_gt3_{c}'].dropna().values for c in cohs] + [U[~U.universal.astype(bool)][f'frac_w_gt3_{c}'].dropna().values for c in cohs]
pos = [0, 1, 2, 3.6, 4.6, 5.6]; bp = ax.boxplot(data, positions=pos, widths=0.6, patch_artist=True, showfliers=False, medianprops=dict(color='k', lw=0.8), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6), boxprops=dict(lw=0))
for p_, c_ in zip(bp['boxes'], [OI['orange']] * 3 + ['0.7'] * 3): p_.set_facecolor(c_)
ax.set_xticks(pos); ax.set_xticklabels([c.replace('GSE', '') for c in cohs] * 2, rotation=45, ha='right', fontsize=6); ax.axhline(0.9, color=OI['red'], ls='--', lw=0.6)
ax.text(1, 1.06, 'universal peaks', ha='center', fontsize=6, color=OI['orange']); ax.text(4.6, 1.06, 'other frequencies', ha='center', fontsize=6, color='0.4')
ax.set_ylim(0, 1.12); ax.set_ylabel('Fraction of biopsies\n> 3 × background'); ax.set_xlabel('Cohort (GSE)'); lab(ax, 'F', -0.3)
save(fig, 'SuppFig_S1_spectral_method_and_controls')

# ======================= S3 — stage effects on the spectrum and the gene-order control
fig = plt.figure(figsize=(W2, W2 * 0.6)); gs = fig.add_gridspec(2, 3, hspace=0.8, wspace=0.5)
Pb = pd.read_csv(tab('Table_S2e_band_profile_by_condition_log.csv'), index_col=0).reindex(ST)
ax = fig.add_subplot(gs[0, 0])
for s_ in ST: ax.plot(range(Pb.shape[1]), np.exp(Pb.loc[s_].values.astype(float) + EULER), 'o-', color=pal[s_], ms=2.2, lw=0.8, label=s_)
ax.axhline(1, color='0.6', ls=':', lw=0.6); ax.set_xticks(range(Pb.shape[1])); ax.set_xticklabels(blab(Pb.columns), rotation=45, ha='right', fontsize=6)
ax.set_xlabel('Period (genes)'); ax.set_ylabel('Power / null'); ax.legend(ncol=2, fontsize=6, columnspacing=0.6, handlelength=1.0); lab(ax, 'A')
Bd = pd.read_csv(tab('Table_S2d_stage_effect_per_band_log.csv'))
ax = fig.add_subplot(gs[0, 1]); ax.bar(range(len(Bd)), Bd.t_stage, color=[OI['red'] if t > 2 else (OI['blue'] if t < -2 else '0.7') for t in Bd.t_stage], lw=0, width=0.7)
ax.axhline(2, color='0.5', ls=':', lw=0.6); ax.axhline(-2, color='0.5', ls=':', lw=0.6); ax.axhline(0, color='k', lw=0.4)
ax.set_xticks(range(len(Bd))); ax.set_xticklabels(blab(Bd.iloc[:, 0]), rotation=45, ha='right', fontsize=6); ax.set_xlabel('Period (genes)'); ax.set_ylabel('t of stage on band power'); lab(ax, 'B')
C6 = pd.read_csv(tab('Table_S6b_composition_scores_per_sample.csv'), index_col=0); C6 = C6[C6.estadio.isin(ST)]
ax = fig.add_subplot(gs[0, 2]); bp = ax.boxplot([C6[C6.estadio == s_].offset.values for s_ in ST], tick_labels=ST, patch_artist=True, showfliers=False, widths=0.6,
                                               medianprops=dict(color='white', lw=1), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6), boxprops=dict(lw=0))
for p_, s_ in zip(bp['boxes'], ST): p_.set_facecolor(pal[s_])
SP = pd.read_csv(tab('Table_S3c_specparam_vs_stage.csv')).set_index('parameter')
ax.set_ylabel('Aperiodic offset'); ax.set_xlabel('Fibrosis stage'); ax.tick_params(axis='x', labelsize=6.5)
ax.text(0.97, 0.97, f"stage t = {SP.loc['offset', 't']:.1f}", transform=ax.transAxes, ha='right', va='top', fontsize=6.3); lab(ax, 'C')
E3 = pd.read_csv(tab('Table_S3d_exponent_per_chromosome_vs_stage.csv')); tcol = [c for c in E3.columns if c.startswith('t')][0]; ccol = E3.columns[0]
ax = fig.add_subplot(gs[1, 0:2]); ax.bar(range(len(E3)), E3[tcol], color=[OI['red'] if t > 2 else (OI['blue'] if t < -2 else '0.75') for t in E3[tcol]], lw=0, width=0.7)
ax.axhline(0, color='k', lw=0.4); ax.set_xticks(range(len(E3))); ax.set_xticklabels(E3[ccol].astype(str), fontsize=6.5); ax.set_xlabel('Chromosome'); ax.set_ylabel('t of stage on aperiodic exponent')
ax.text(0.01, 0.95, f"genome-wide exponent: P = {SP.loc['exponent', 'p']:.2f}", transform=ax.transAxes, va='top', fontsize=6.3); lab(ax, 'D', -0.1)
F11 = pd.read_csv(tab('Table_S11_spectral_fingerprint_real_vs_shuffled.csv')); tasks = ['early_vs_advanced', 'normal_vs_advanced', 'ordinal']
feats = [('spectrum', 'Gene order', OI['blue']), ('spectrum_shuffled_order', 'Shuffled order', OI['sky']), ('expression_2000', 'Expression', '0.6')]
ax = fig.add_subplot(gs[1, 2])
for j, (f_, nm, col) in enumerate(feats):
    vv = [F11[(F11.features == f_) & (F11.task == t_)].value.iloc[0] for t_ in tasks]; ax.bar(np.arange(3) + (j - 1) * 0.26, vv, 0.26, color=col, lw=0, label=nm)
ax.set_xticks(range(3)); ax.set_xticklabels(['F0–F1 vs\nF3–F4 (AUC)', 'Normal vs\nF3–F4 (AUC)', 'Stage\n(Spearman)'], fontsize=6); ax.set_ylim(0, 1.18); ax.axhline(0.5, color='0.5', ls=':', lw=0.6)
ax.set_ylabel('Leave-one-cohort-out'); ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.22), ncol=3, fontsize=5.8, columnspacing=0.6, handlelength=1.0); lab(ax, 'E', -0.3, 1.28)
save(fig, 'SuppFig_S3_stage_effects_on_spectrum')

# ======================= S4 — orientation, distance and TADs
fig, axs = plt.subplots(1, 4, figsize=(W2, W2 * 0.3)); plt.subplots_adjust(wspace=0.65)
O5 = pd.read_csv(tab('Table_S5h_concordance_by_orientation.csv'))
ax = axs[0]; x_ = np.arange(len(O5)); ax.bar(x_ - 0.18, O5.r_t, 0.36, color='k', lw=0, label='Stage effect'); ax.bar(x_ + 0.18, O5.coexpr, 0.36, color=OI['blue'], lw=0, label='Co-expression')
ax.set_xticks(x_); ax.set_xticklabels(O5.orientation, rotation=30, ha='right'); ax.set_ylabel('Correlation, adjacent pairs'); ax.set_ylim(0, 0.42); ax.legend(loc='upper right', fontsize=6); lab(ax, 'A', -0.42)
Pp = pd.read_pickle(inter('pairs.pkl')).rename(columns={'intergenic_bp': 'dist_bp', 'same_tad': 'mismo_tad'}); q2 = Pp.dropna(subset=['stat1', 'stat2']).copy(); bins = [0, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, 5e6]; q2['db'] = pd.cut(q2.dist_bp, bins)
ax = axs[1]
for same, col, nm in [(True, OI['red'], 'Same TAD'), (False, OI['blue'], 'Different TADs')]:
    g_ = q2[q2.mismo_tad == same].groupby('db', observed=True); r_ = g_.apply(lambda h: np.corrcoef(h.stat1, h.stat2)[0, 1] if len(h) >= 60 else np.nan)
    ax.plot([b_.right / 1e3 for b_ in r_.index], r_.values, 'o-', color=col, ms=2.6, lw=0.9, label=nm)
ax.set_xscale('log'); ax.axhline(0, color='0.6', lw=0.4); ax.set_xlabel('Distance between neighbours (kb)'); ax.set_ylabel('Concordance of stage effect'); ax.legend(loc='upper right', fontsize=6); lab(ax, 'B', -0.42)
T5 = pd.read_csv(tab('Table_S5i_TAD_test_stratified.csv')).iloc[0]
ax = axs[2]; ax.errorbar([0], [T5.null_mean], yerr=[2 * T5.null_sd], fmt='o', color='0.5', capsize=4, label='Permutation null ± 2 SD'); ax.scatter([0], [T5.observed_delta_r], s=30, color=OI['red'], zorder=4, label='Observed')
ax.set_xlim(-1, 1); ax.set_xticks([]); ax.set_ylim(0, 0.3); ax.set_ylabel('Same-TAD minus different-TAD\nconcordance, distance-matched'); ax.text(0.5, 0.04, f'P = {T5.p:.2f}', transform=ax.transAxes, ha='center', fontsize=6.5); ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), fontsize=6); lab(ax, 'C', -0.42, 1.3)
K5 = pd.read_csv(tab('Table_S5k_neighbourhoods_within_one_TAD.csv'))
ax = axs[3]; ax.bar(range(len(K5)), 100 * K5.observed_within_one_TAD, color='0.55', lw=0, width=0.6, label='Observed')
ax.errorbar(range(len(K5)), 100 * K5.chance_mean, yerr=200 * K5.chance_sd, fmt='o', color='k', ms=3, capsize=3, lw=0.7, label='Random windows ± 2 SD')
ax.set_xticks(range(len(K5))); ax.set_xticklabels([f'≥{m}' for m in K5.min_genes]); ax.set_xlabel('Genes per DE neighbourhood'); ax.set_ylabel('% within one TAD'); ax.set_ylim(0, 100); ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), fontsize=6); lab(ax, 'D', -0.42, 1.3)
save(fig, 'SuppFig_S4_orientation_distance_TADs')

# ======================= S5 — shared liver eQTL variants
E13 = pd.read_csv(tab('Table_S13e_shared_GTEx_eQTL_coupled_vs_uncoupled.csv')); F13 = pd.read_csv(tab('Table_S13f_concordance_by_shared_eQTL_distance.csv')); G13 = pd.read_csv(tab('Table_S13g_models_concordance_on_shared_eQTL.csv'))
E13all = E13.copy(); E13 = E13[E13.coupling == E13.coupling.iloc[0]].reset_index(drop=True)   # stage coupling (|t| > 3), the primary definition
fig, axs = plt.subplots(1, 4, figsize=(W2, W2 * 0.34), gridspec_kw={'width_ratios': [1.0, 1.1, 1.0, 1.1]}); plt.subplots_adjust(wspace=0.95)
ax = axs[0]; x_ = np.arange(len(E13))
ax.bar(x_ - 0.18, 100 * E13.frac_coupled, 0.36, color=OI['red'], lw=0, label=f'Coupled (n = {int(E13.n_coupled.iloc[0])})')
ax.bar(x_ + 0.18, 100 * E13.frac_uncoupled, 0.36, color='0.6', lw=0, label=f'Uncoupled (n = {int(E13.n_uncoupled.iloc[0])})')
ax.set_xticks(x_); ax.set_xticklabels(['Any', 'Same\ndirection', 'Lead\nvariant'], fontsize=6.3); ax.set_xlabel('Shared significant eQTL variant'); ax.set_ylabel('% of adjacent eGene pairs'); ax.set_ylim(0, 60)
ax.legend(loc='upper right', fontsize=6); lab(ax, 'A', -0.3)
ax = axs[1]; dists = list(dict.fromkeys(F13.distance))
for sh, col, nm in [(True, OI['red'], 'Share same-direction variant'), (False, '0.5', 'No shared variant')]:
    d_ = F13[F13.shares_same_dir == sh].set_index('distance').reindex(dists); ax.plot(range(len(dists)), d_.r_stage_effect, 'o-', color=col, ms=3, lw=0.9, label=nm)
ax.set_xticks(range(len(dists))); ax.set_xticklabels([d.replace('kb', '') for d in dists], rotation=50, ha='right', fontsize=6); ax.set_xlabel('Intergenic distance (kb)'); ax.set_ylabel('Concordance of stage effect'); ax.axhline(0, color='0.7', lw=0.4); ax.legend(fontsize=6); lab(ax, 'B', -0.3)
ax = axs[2]; nm_ = {'shares': 'Any shared variant', 'shares_same_dir': 'Same-direction variant', 'lead_shared': 'Shared lead variant'}
y_ = np.arange(len(G13)); ax.barh(y_, G13.delta_coexpr, color=[OI['red'] if p < 0.05 else '0.65' for p in G13.p_coexpr], lw=0, height=0.6)
for i, r in G13.iterrows(): ax.text(r.delta_coexpr + 0.003, i, f'P = {r.p_coexpr:.1g}', va='center', fontsize=6)
ax.set_yticks(y_); ax.set_yticklabels([nm_.get(p, p).replace(' variant', '\nvariant') for p in G13.predictor], fontsize=6.3); ax.invert_yaxis(); ax.set_xlabel('Increase in co-expression across biopsies\n(distance-adjusted)'); ax.set_xlim(0, 0.13); lab(ax, 'C', -0.62)
ax = axs[3]; F = E13all.reset_index(drop=True); y_ = np.arange(len(F))[::-1]
short_o = {'share >=1 significant variant': 'any', 'share variant with same-direction effect': 'same direction', 'lead variant of one is significant for the other': 'lead'}
for yy, (_, r) in zip(y_, F.iterrows()):
    col = OI['red'] if r.coupling.startswith('coupled (stage') else OI['blue']
    ax.errorbar(r.MH_OR_distance_stratified, yy, xerr=[[r.MH_OR_distance_stratified - r.CI_low], [r.CI_high - r.MH_OR_distance_stratified]], fmt='o', color=col, ms=3, capsize=2, lw=0.8)
ax.axvline(1, color='k', lw=0.5, ls=':'); ax.set_xscale('log'); ax.set_yticks(y_); ax.set_yticklabels([short_o.get(o, o) for o in F.outcome], fontsize=6)
ax.set_xlabel('Distance-stratified odds ratio'); ax.text(1.02, 0.78, 'stage\ncoupling\n(118 pairs)', transform=ax.transAxes, color=OI['red'], fontsize=6, va='center')
ax.text(1.02, 0.25, 'after\ncomposition\n(52 pairs)', transform=ax.transAxes, color=OI['blue'], fontsize=6, va='center'); lab(ax, 'D', -0.5)
save(fig, 'SuppFig_S5_shared_eQTL_variants')

# ======================= copy figures produced elsewhere under their supplementary numbers
for src, dst in [(os.path.join(FIG, 'genome_research', 'ExtendedData_Fig2_gtex_tissues'), 'SuppFig_S2_GTEx_tissue_specificity'),
                 (os.path.join(FIG, 'genome_research', 'ExtendedData_Fig1_mouse_validation'), 'SuppFig_S6_mouse_replication'),
                 (os.path.join(FIG, 'jhep', 'SuppFig_GWAS_loci_neighbourhoods'), 'SuppFig_S7_GWAS_loci'),
                 (os.path.join(FIG, 'genome_research', 'Fig7_circos_TF_Hallmark'), 'SuppFig_S8_TF_Hallmark_map'),
                 (os.path.join(FIG, 'jhep', 'SuppFig_snRNAseq_validation'), 'SuppFig_S9_snRNAseq_validation')]:
    for ext in ('pdf', 'png'):
        if os.path.exists(src + '.' + ext): shutil.copy(src + '.' + ext, os.path.join(SO, dst + '.' + ext))
print('supplementary figures:', sorted(f for f in os.listdir(SO) if f.endswith('.pdf')))
print(f'mask r = {r_m:.3f}')
