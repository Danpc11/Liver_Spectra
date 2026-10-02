"""11 — Validation in own mouse fatty-liver RNA-seq (CTL vs UTP, n = 3 + 3; ARC-UTP samples excluded).
Tests in a second species and an independent perturbation: (i) invariance of the positional architecture,
(ii) spatial coupling of the transcriptional change, (iii) distance decay of neighbour concordance,
(iv) whether human fibrosis-coupled neighbour pairs are specifically coupled in mouse (syntenic pairs).
Input: data/raw/own/mcounts.tsv. Outputs: Table_S12_*.csv, ExtendedData_Fig1_mouse_validation.{pdf,png}
"""
import os, re, numpy as np, pandas as pd, pyannotables as pa
from scipy import stats
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from common import *

SAMPLES = ['CTL1', 'CTL2', 'CTL3', 'UTP1', 'UTP2', 'UTP3']
MCHR = [str(i) for i in range(1, 20)] + ['X']
m = pd.read_csv(os.path.join(RAW, 'own', 'mcounts.tsv'), sep='\t', index_col=0)
m.index = m.index.str.split('.').str[0]; m = m.groupby(level=0).sum()[SAMPLES]
ann = pa.tables()['mus_musculus-GRCm38-ensembl100']; ann = ann[~ann.index.duplicated()]
ann = ann[ann.Chromosome.astype(str).isin(MCHR)][['Chromosome', 'Start', 'End', 'Strand', 'gene_name']].copy(); ann['chr'] = ann.Chromosome.astype(str)
keep = m.index[m.median(axis=1) >= 10].intersection(ann.index)
A = ann.loc[keep].sort_values(['chr', 'Start']).copy(); A['pos'] = A.groupby('chr').cumcount() + 1
print('mouse genes (expressed, with coordinates):', len(A))

# ---- DE
md = pd.DataFrame({'g': ['CTL'] * 3 + ['UTP'] * 3}, index=SAMPLES)
dds = DeseqDataSet(counts=m.loc[A.index].T.astype(int), metadata=md, design='~ g', quiet=True); dds.deseq2()
st = DeseqStats(dds, contrast=['g', 'UTP', 'CTL'], quiet=True); st.summary(); de = st.results_df; de['gene_name'] = A.gene_name.reindex(de.index).values
de.to_csv(tab('Table_S12a_mouse_DESeq2_UTP_vs_CTL.csv')); print('DE padj<0.05:', int((de.padj < .05).sum()))

# ---- invariance of positional architecture
lc = np.log2(m.loc[A.index] / m.loc[A.index].sum() * 1e6 + 1); rows, names = [], []
for c in MCHR:
    g = A[A.chr == c]; v = lc.loc[g.index].values; v = v - v.mean(0); n = len(g)
    P = np.abs(np.fft.rfft(v, axis=0)[1:n // 2 + 1]) ** 2; rows.append(whiten_matrix(P, n)); names += [(c, k) for k in range(1, P.shape[0] + 1)]
W = pd.DataFrame(np.vstack(rows), index=pd.MultiIndex.from_tuples(names, names=['chr', 'k']), columns=SAMPLES)
u = (W > 3).all(axis=1)
S = pd.DataFrame({'period': [len(A[A.chr == c]) / k for c, k in W.index], 'mean_power_CTL': np.exp(np.log(W[SAMPLES[:3]]).mean(1) + EULER), 'mean_power_UTP': np.exp(np.log(W[SAMPLES[3:]]).mean(1) + EULER), 'universal_all6': u}, index=W.index)
S.reset_index().round(4).to_csv(tab('Table_S12b_mouse_spectrum_universal_peaks.csv'), index=False)
tU = stats.ttest_rel(np.log(S.loc[u, 'mean_power_UTP']), np.log(S.loc[u, 'mean_power_CTL']))
print(f'universal peaks (>3x in all 6): {u.sum()} of {len(W)}; median power CTL {S.loc[u, "mean_power_CTL"].median():.2f}x, UTP {S.loc[u, "mean_power_UTP"].median():.2f}x; paired t {tU.statistic:.2f}, P={tU.pvalue:.2g}')

# ---- spatial coupling of the DE statistic
s = de.stat.reindex(A.index).fillna(0); byc = [s[A.chr == c].values for c in MCHR]; rng = np.random.default_rng(0); lags = (1, 2, 3, 5, 10)
obs = acf_by_chr(byc, lags); null = np.array([acf_by_chr([rng.permutation(v) for v in byc], lags) for _ in range(200)])
ACF = pd.DataFrame({'lag': lags, 'acf_observed': obs, 'null_p2.5': np.quantile(null, .025, 0), 'null_p97.5': np.quantile(null, .975, 0), 'z': (obs - null.mean(0)) / null.std(0)})
ACF.round(4).to_csv(tab('Table_S12c_mouse_spatial_autocorrelation.csv'), index=False); print(ACF.round(3).to_string(index=False))
sig = np.sign(de.stat) * (de.padj < 0.05).astype(int); sig = sig.reindex(A.index).fillna(0)
def adj(sv): return sum(int(np.sum((sv[A.chr.values == c][:-1] == sv[A.chr.values == c][1:]) & (sv[A.chr.values == c][:-1] != 0))) for c in MCHR)
o = adj(sig.values); nl = [adj(np.concatenate([rng.permutation(sig.values[A.chr.values == c]) for c in MCHR])) for _ in range(200)]
# note: permuted vector concatenated in chromosome order matches A ordering (A sorted by chr, Start)
print(f'adjacent concordant DE pairs: {o} vs null {np.mean(nl):.0f} ± {np.std(nl):.0f} (z={(o - np.mean(nl)) / np.std(nl):.1f})')

# ---- distance decay (non-paralog neighbours)
rows = []
for c, g in A.groupby('chr'):
    g = g.sort_values('pos'); a = g.iloc[:-1]; b = g.iloc[1:]
    ibp = np.maximum(b.Start.values - a.End.values, a.Start.values - b.End.values)
    fam = [re.sub(r'[0-9]+[A-Za-z]?$', '', str(x)) == re.sub(r'[0-9]+[A-Za-z]?$', '', str(y)) and len(re.sub(r'[0-9]+[A-Za-z]?$', '', str(x))) >= 3 for x, y in zip(a.gene_name, b.gene_name)]
    rows.append(pd.DataFrame({'g1': a.index, 'g2': b.index, 'intergenic_bp': ibp, 'paralog': fam, 's1': s.reindex(a.index).values, 's2': s.reindex(b.index).values}))
PM = pd.concat(rows); PM = PM[~PM.paralog]; PM['dbin'] = pd.cut(PM.intergenic_bp, [-1e9, 0, 1e3, 1e4, 5e4, 1e5, 5e5, 1e9], labels=['overlap', '<1 kb', '1-10 kb', '10-50 kb', '50-100 kb', '100-500 kb', '>500 kb'])
DEC = PM.groupby('dbin', observed=True).apply(lambda h: pd.Series({'n': len(h), 'r_DE': np.corrcoef(h.s1, h.s2)[0, 1]})); DEC.round(4).to_csv(tab('Table_S12d_mouse_concordance_by_distance.csv')); print(DEC.round(3).to_string())

# ---- syntenic human pairs: specific or generic coupling?
P = pd.read_csv(tab('Table_S5f_adjacent_pairs_distance_orientation_TAD.csv')); P = P[~P.paralog_family]
mpos = A.assign(sym=A.gene_name.str.upper()).drop_duplicates('sym').set_index('sym')[['chr', 'pos']]
P = P[P.n1.isin(mpos.index) & P.n2.isin(mpos.index)].copy()
P['same'] = (mpos.chr.reindex(P.n1).values == mpos.chr.reindex(P.n2).values) & (np.abs(mpos.pos.reindex(P.n1).values - mpos.pos.reindex(P.n2).values) <= 2)
syn = P[P.same].copy(); ms = de.assign(sym=de.gene_name.str.upper()).dropna(subset=['sym']).drop_duplicates('sym').set_index('sym').stat
syn['m1'] = ms.reindex(syn.n1).values; syn['m2'] = ms.reindex(syn.n2).values; syn = syn.dropna(subset=['m1', 'm2'])
syn['human_coupled'] = (np.sign(syn.t1) == np.sign(syn.t2)) & (syn.t1.abs() > 3) & (syn.t2.abs() > 3)
rnd = [np.corrcoef(rng.choice(ms.dropna().values, len(syn)), rng.choice(ms.dropna().values, len(syn)))[0, 1] for _ in range(300)]
def rboot(df, B=500):
    r = [np.corrcoef(*df.sample(len(df), replace=True, random_state=i)[['m1', 'm2']].values.T)[0, 1] for i in range(B)]; return np.percentile(r, [2.5, 97.5])
c_ = syn[syn.human_coupled]; n_ = syn[~syn.human_coupled]
SY = pd.DataFrame([{'set': 'all syntenic neighbour pairs', 'n': len(syn), 'r_mouse': np.corrcoef(syn.m1, syn.m2)[0, 1], 'ci': rboot(syn)},
                   {'set': 'coupled in human fibrosis', 'n': len(c_), 'r_mouse': np.corrcoef(c_.m1, c_.m2)[0, 1], 'ci': rboot(c_)},
                   {'set': 'not coupled in human', 'n': len(n_), 'r_mouse': np.corrcoef(n_.m1, n_.m2)[0, 1], 'ci': rboot(n_)},
                   {'set': 'random gene pairs', 'n': len(syn), 'r_mouse': np.mean(rnd), 'ci': np.percentile(rnd, [2.5, 97.5])}])
SY['ci_low'] = SY.ci.str[0]; SY['ci_high'] = SY.ci.str[1]; SY.drop(columns='ci').round(4).to_csv(tab('Table_S12e_mouse_syntenic_pairs.csv'), index=False); print(SY.drop(columns='ci').round(3).to_string(index=False))

# ---- relation to human signature
H = pd.read_csv(tab('Table_S6c_stage_effect_per_gene_with_without_composition.csv')).drop_duplicates('gene_name').set_index('gene_name')
j = pd.DataFrame({'mouse_stat': ms}).join(H[['t_stage', 't_stage_comp_adjusted']], how='inner')
print('Spearman mouse UTP-vs-CTL vs human stage: %.3f; vs composition-adjusted: %.3f (n=%d orthologs by symbol)' % (stats.spearmanr(j.mouse_stat, j.t_stage, nan_policy='omit')[0], stats.spearmanr(j.mouse_stat, j.t_stage_comp_adjusted, nan_policy='omit')[0], len(j)))
pd.to_pickle((S, ACF, DEC, SY), inter('mouse_validation.pkl'))

# ---- Extended Data Figure 1
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, scienceplots
plt.style.use(['science', 'nature', 'no-latex'])
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial', 'Helvetica'], 'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6,
                     'xtick.top': False, 'ytick.right': False, 'xtick.minor.visible': False, 'ytick.minor.visible': False, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False, 'pdf.fonttype': 42})
OI = {'blue': '#0072B2', 'orange': '#E69F00', 'red': '#D55E00', 'grey': '#7F7F7F', 'green': '#009E73'}
def lab(ax, t, dx=-0.22, dy=1.16): ax.text(dx, dy, t.upper(), transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')
Hacf = pd.read_csv(tab('Table_S5a_spatial_autocorrelation_DE.csv'))
fig, axs = plt.subplots(2, 2, figsize=(183 / 25.4 * 0.8, 183 / 25.4 * 0.7)); axs = axs.ravel(); plt.subplots_adjust(hspace=0.65, wspace=0.45)
ax = axs[0]; x = np.log10(S.mean_power_CTL); y = np.log10(S.mean_power_UTP)
ax.hexbin(x, y, gridsize=45, cmap='Greys', mincnt=1, linewidths=0); ax.scatter(x[u], y[u], s=4, color=OI['red'], lw=0, label=f'universal peaks (n={u.sum()})')
lim = [-1.2, 1.3]; ax.plot(lim, lim, color='k', lw=0.5, ls=':'); ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xticks([-1, 0, 1]); ax.set_xticklabels(['0.1', '1', '10']); ax.set_yticks([-1, 0, 1]); ax.set_yticklabels(['0.1', '1', '10'])
ax.set_xlabel('Power / background, CTL'); ax.set_ylabel('Power / background, UTP'); ax.legend(loc='upper left'); ax.text(0.98, 0.04, f'peaks: paired P = {tU.pvalue:.2f}', transform=ax.transAxes, ha='right', fontsize=6); lab(ax, 'a')
ax = axs[1]; ax.axhspan(ACF['null_p2.5'].min(), ACF['null_p97.5'].max(), color='0.9', lw=0, label='95% permutation null')
ax.plot(ACF.lag, ACF.acf_observed, 'o-', color=OI['blue'], ms=3, lw=0.9, label='Mouse, UTP vs CTL'); hh = Hacf[Hacf.lag <= 10]; ax.plot(hh.lag, hh.acf_observed, 's--', color='k', ms=3, lw=0.9, label='Human, fibrosis stage')
ax.set_xscale('log'); ax.set_xlabel('Distance (genes)'); ax.set_ylabel('Spatial autocorrelation\nof the DE statistic'); ax.legend(loc='upper right'); lab(ax, 'b')
ax = axs[2]; ax.plot(range(len(DEC)), DEC.r_DE, 'o-', color=OI['blue'], ms=3, lw=0.9); ax.axhline(0, color='0.6', lw=0.5)
ax.set_xticks(range(len(DEC))); ax.set_xticklabels(['overlap', '<1', '1–10', '10–50', '50–100', '100–500', '>500'], rotation=45, ha='right'); ax.set_xlabel('Intergenic distance (kb), mouse non-paralog neighbours'); ax.set_ylabel('Concordance of UTP effect'); lab(ax, 'c')
ax = axs[3]; cols = [OI['grey'], OI['red'], OI['orange'], '0.75']
ax.bar(range(4), SY.r_mouse, color=cols, lw=0, width=0.65); ax.errorbar(range(4), SY.r_mouse, yerr=[SY.r_mouse - SY.ci_low, SY.ci_high - SY.r_mouse], fmt='none', color='k', lw=0.6, capsize=2)
ax.set_xticks(range(4)); ax.set_xticklabels(['all\nsyntenic', 'coupled\nin human', 'not coupled\nin human', 'random\npairs'], fontsize=6); ax.set_ylabel('Concordance of UTP effect\nbetween neighbours (r)'); ax.axhline(0, color='0.6', lw=0.5); lab(ax, 'd')
fig.savefig(os.path.join(FIG, 'genome_research', 'ExtendedData_Fig1_mouse_validation.pdf'), bbox_inches='tight'); fig.savefig(os.path.join(FIG, 'genome_research', 'ExtendedData_Fig1_mouse_validation.png'), dpi=400, bbox_inches='tight'); print('figure ok')
