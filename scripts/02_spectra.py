"""02 — Positional spectra per sample (log-expression and within-cohort z), whitening, permutation null,
universal peaks, grid-mask control, consensus spectra per condition, stage effects per frequency/band, multitaper.
Outputs: Table_S2_*.csv, intermediate W_adj.pkl, W_mt.pkl
"""
import numpy as np, pandas as pd
from scipy import stats
from scipy.signal.windows import dpss
from statsmodels.stats.multitest import multipletests
from common import *

X, keep = pd.read_pickle(inter('expr.pkl')); A = pd.read_pickle(inter('expr_adj.pkl')); M = pd.read_pickle(inter('meta.pkl'))
grid = load_grid(); N = chrom_lengths(grid)

# within-cohort gene standardisation (deviation spectra)
Zg = pd.concat([((A[M.index[M.cohorte == c]].T - A[M.index[M.cohorte == c]].mean(axis=1)) / A[M.index[M.cohorte == c]].std(axis=1).replace(0, 1)).T for c in COHORTS], axis=1)[M.index]
Wl = spectra(A, keep, N); Wz = spectra(Zg, keep, N)
tapers = {n: dpss(n, 3, 5) for n in set(N.values())}
Wmt = spectra(A, keep, N, tapers=tapers)
pd.to_pickle((Wz, Wl), inter('W_adj.pkl')); pd.to_pickle(Wmt, inter('W_mt.pkl'))
per = pd.Series([N[c] / k for c, k in Wl.index], index=Wl.index)

# permutation null for the whitened periodogram
rng = np.random.default_rng(0); sub = A.iloc[:, rng.choice(A.shape[1], 60, replace=False)]
null = np.concatenate([np.log(spectra(sub, keep, N, rng=rng).values.ravel()) for _ in range(5)])
print('null: mean log w %.3f (theory -0.577), sd %.3f, 99th pct w %.2f' % (null.mean(), null.std(), np.exp(np.quantile(null, .99))))

# grid-mask control
mask_rows = []
for c in CHR:
    g = keep[keep.chr == c]; m = np.zeros(N[c]); m[g.grid_index.values - 1] = 1
    P = np.abs(np.fft.rfft(m - m.mean())[1:N[c] // 2 + 1]) ** 2
    w = whiten_matrix(P[:, None], N[c])[:, 0]; mask_rows += [(c, k + 1, w[k]) for k in range(len(w))]
mask = pd.DataFrame(mask_rows, columns=['chr', 'k', 'w_mask']).set_index(['chr', 'k']).w_mask

# universal peaks
Lw = np.log(Wl); m = Lw.mean(axis=1); frac = (Wl > 3).mean(axis=1)
T = pd.DataFrame({'period': per, 'mean_log_w': m, 'mean_power_over_null': np.exp(m + EULER), 'frac_samples_w_gt3': frac, 'w_mask': mask.reindex(Wl.index)})
for c in COHORTS: T['frac_w_gt3_' + c] = (Wl[M.index[M.cohorte == c]] > 3).mean(axis=1)
T['universal'] = T.frac_samples_w_gt3 >= 0.9
T.reset_index().round(4).to_csv(tab('Table_S2a_consensus_spectrum_universal_peaks.csv'), index=False)
print('universal peaks:', int(T.universal.sum()), '| in every cohort:', int((T[[f'frac_w_gt3_{c}' for c in COHORTS]] >= 0.8).all(axis=1).sum()),
      '| corr(log mean w, log mask):', round(np.corrcoef(np.log(T.mean_power_over_null), np.log(T.w_mask))[0, 1], 3))

# consensus spectra per condition and stage models, for log and z spectra
s = M.index[M.estadio != 'Control']
for lab, W in [('log', Wl), ('z', Wz), ('multitaper', Wmt)]:
    Lw = np.log(W); rows = []
    for st in STAGES + ['Control']:
        ss = M.index[M.estadio == st]; Ls = Lw[ss]; mu = Ls.mean(axis=1); se = Ls.std(axis=1) / np.sqrt(len(ss))
        t = (mu + EULER) / se; p = 2 * stats.t.sf(np.abs(t), len(ss) - 1)
        rows.append(pd.DataFrame({'estadio': st, 'n': len(ss), 'log_w_mean': mu, 'se': se, 't_vs_null': t, 'q': multipletests(p, method='fdr_bh')[1], 'power_over_null': np.exp(mu + EULER)}))
    R = pd.concat(rows).reset_index(); R['period'] = per.loc[pd.MultiIndex.from_frame(R[['chr', 'k']])].values
    R.round(4).to_csv(tab(f'Table_S2b_condition_spectra_{lab}.csv'), index=False)
    D = design(M, s); t, df = ols_t(Lw[s].values.T, D); p = 2 * stats.t.sf(np.abs(t), df); q = multipletests(p, method='fdr_bh')[1]
    E = pd.DataFrame({'slope_t': t, 'p': p, 'q': q, 'period': per.values}, index=Lw.index).reset_index()
    E.round(5).to_csv(tab(f'Table_S2c_stage_effect_per_frequency_{lab}.csv'), index=False)
    band = pd.cut(per, BANDS); B = Lw.groupby(band.values, observed=True).mean().T
    out = []
    for b in B.columns:
        tb, _ = ols_t(B.loc[s, [b]].values, D); out.append({'band': str(b), 't_stage': tb[0]})
    pd.DataFrame(out).to_csv(tab(f'Table_S2d_stage_effect_per_band_{lab}.csv'), index=False)
    Bm = B.groupby(M.estadio).mean().reindex(STAGES + ['Control']); Bm.columns = [str(c) for c in Bm.columns]
    Bm.round(4).to_csv(tab(f'Table_S2e_band_profile_by_condition_{lab}.csv'))
    print(f'[{lab}] stage-associated frequencies q<0.05: {(q < 0.05).sum()}; sd between samples {Lw.std(axis=1).median():.3f}')
