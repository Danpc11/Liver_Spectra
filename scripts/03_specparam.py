"""03 — Spectral parameterisation per sample and chromosome (aperiodic offset/exponent + periodic peaks),
after Donoghue et al. 2020; stage models. Outputs: Table_S3_specparam_*.csv
"""
import numpy as np, pandas as pd, statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from common import *

X, keep = pd.read_pickle(inter('expr.pkl')); A = pd.read_pickle(inter('expr_adj.pkl')); M = pd.read_pickle(inter('meta.pkl'))
grid = load_grid(); N = chrom_lengths(grid)


def robust_aperiodic(lx, L, passes=2, thr=2.0):
    mask = np.ones_like(L, dtype=bool); b = np.zeros(L.shape[1]); a = np.zeros(L.shape[1])
    for _ in range(passes + 1):
        for j in range(L.shape[1]):
            mm = mask[:, j]; b[j], a[j] = np.polyfit(lx[mm], L[mm, j], 1)
        res = L - (np.outer(lx, b) + a); sd = res.std(axis=0); mask = res < thr * sd
    return a, b, res

rows, peaks = [], []
for c in CHR:
    g = keep[keep.chr == c].sort_values('grid_index'); vals = A.loc[g.gene_id].values; n = N[c]
    S = np.zeros((n, vals.shape[1])); S[g.grid_index.values - 1] = vals - vals.mean(axis=0)
    P = np.abs(np.fft.rfft(S, axis=0)[1:n // 2 + 1]) ** 2
    Ps = pd.DataFrame(P).rolling(5, center=True, min_periods=1).mean().values
    f = np.arange(1, P.shape[0] + 1) / n; lx = np.log10(f); L = np.log10(Ps + 1e-12)
    a, b, res = robust_aperiodic(lx, L); sd = res.std(axis=0)
    ispk = (res > 2.5 * sd) & (res >= np.roll(res, 1, axis=0)) & (res >= np.roll(res, -1, axis=0)); ispk[[0, -1]] = False
    rows.append(pd.DataFrame({'sample': A.columns, 'chr': c, 'offset': a, 'exponent': -b, 'sd_resid': sd, 'n_peaks': ispk.sum(axis=0)}))
    ii, jj = np.where(ispk); peaks.append(pd.DataFrame({'sample': A.columns[jj], 'chr': c, 'k': ii + 1, 'period': n / (ii + 1), 'height_log10': res[ii, jj]}))
SP = pd.concat(rows); PK = pd.concat(peaks)
SP.to_csv(tab('Table_S3a_specparam_per_sample_chromosome.csv'), index=False); PK.to_csv(tab('Table_S3b_specparam_peaks_per_sample.csv'), index=False)

s = M.index[M.estadio != 'Control']; D = design(M, s); j = list(D.columns).index('orden'); out = []
for var in ['offset', 'exponent', 'n_peaks', 'sd_resid']:
    y = SP.groupby('sample')[var].mean().loc[s]; m = sm.OLS(y.values, D.values).fit()
    out.append({'parameter': var, 'slope_per_stage': m.params[j], 't': m.tvalues[j], 'p': m.pvalues[j],
                **{st: SP.groupby('sample')[var].mean().loc[M.index[M.estadio == st]].mean() for st in STAGES}})
pd.DataFrame(out).to_csv(tab('Table_S3c_specparam_vs_stage.csv'), index=False); print(pd.DataFrame(out).round(4).to_string(index=False))
rows = []
for c in CHR:
    y = SP[SP.chr == c].set_index('sample').exponent.loc[s]; m = sm.OLS(y.values, D.values).fit(); rows.append({'chr': c, 'slope': m.params[j], 't': m.tvalues[j], 'p': m.pvalues[j]})
EC = pd.DataFrame(rows); EC['q'] = multipletests(EC.p, method='fdr_bh')[1]; EC.to_csv(tab('Table_S3d_exponent_per_chromosome_vs_stage.csv'), index=False)
