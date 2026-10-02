"""05 — Marker-based cell composition; composition-adjusted stage effects; neighbour pairs by genomic distance,
orientation and liver TAD (hg19 / GRCh37); threshold-vs-trend and within-stage bimodality.
Outputs: Table_S6_*, Table_S7_*, Table_S5d/e_*
"""
import os, re, numpy as np, pandas as pd, statsmodels.api as sm, pyannotables as pa
from sklearn.mixture import GaussianMixture
from common import *

X, keep = pd.read_pickle(inter('expr.pkl')); A = pd.read_pickle(inter('expr_adj.pkl')); M = pd.read_pickle(inter('meta.pkl')); de, de2 = pd.read_pickle(inter('de.pkl'))
SP = pd.read_csv(tab('Table_S3a_specparam_per_sample_chromosome.csv'))
names = keep.set_index('gene_id').gene_name.astype(str).reindex(A.index)
Az = ((A.T - A.mean(axis=1)) / A.std(axis=1).replace(0, 1)).T

# ---- composition scores
comp = pd.DataFrame({ct: Az[names.str.match(p).values].mean() for ct, p in MARKERS.items()})
pd.Series({ct: int(names.str.match(p).sum()) for ct, p in MARKERS.items()}).to_csv(tab('Table_S6a_marker_counts.csv'))
C = comp.join(M[['cohorte', 'estadio', 'orden', 'nas', 'sexo', 'edad']]).join(SP.groupby('sample').agg(offset=('offset', 'mean'), n_peaks=('n_peaks', 'mean')))
C.round(4).to_csv(tab('Table_S6b_composition_scores_per_sample.csv')); pd.to_pickle(comp, inter('comp.pkl'))

# ---- stage effect per gene with/without composition
s = M.index[M.estadio != 'Control']; D0 = design(M, s); D1 = design(M, s, extra=comp)
t0, _ = ols_t(A[s].values.T, D0); t1, _ = ols_t(A[s].values.T, D1)
T = pd.DataFrame({'gene_id': A.index, 'gene_name': names.values, 't_stage': t0, 't_stage_comp_adjusted': t1}); T.to_csv(tab('Table_S6c_stage_effect_per_gene_with_without_composition.csv'), index=False)
print('variance of t explained by composition: %.1f%%' % (100 * (1 - np.var(t1) / np.var(t0))))
# offset ~ stage with/without composition
_rows6d = []
for var in ['offset', 'n_peaks']:
    m0 = sm.OLS(C.loc[s, var].values, D0.values).fit(); m1 = sm.OLS(C.loc[s, var].values, D1.values).fit(); j = list(D0.columns).index('orden')
    print(f'{var}: stage t {m0.tvalues[j]:.2f} -> {m1.tvalues[j]:.2f}; R2 {m0.rsquared:.3f} -> {m1.rsquared:.3f}')
    _rows6d.append({'response': var, 'model_without': 'ordinal stage + cohort + sex', 'model_with': '+ six composition scores', 't_stage_without': m0.tvalues[j], 't_stage_with': m1.tvalues[j], 'R2_without': m0.rsquared, 'R2_with': m1.rsquared, 'n': len(s)})
pd.DataFrame(_rows6d).round(4).to_csv(tab('Table_S6d_spectral_parameters_vs_stage_with_composition.csv'), index=False)
# ACF of t with/without composition
Kt = keep.set_index('gene_id').reindex(A.index); order = Kt.sort_values(['chr', 'grid_index']); Ts = T.set_index('gene_id')
rng = np.random.default_rng(0); lags = (1, 2, 3, 5, 10)
byc0 = [Ts.loc[order[order.chr == c].index, 't_stage'].values for c in CHR]; byc1 = [Ts.loc[order[order.chr == c].index, 't_stage_comp_adjusted'].values for c in CHR]
null = np.array([acf_by_chr([rng.permutation(v) for v in byc0], lags) for _ in range(100)])
pd.DataFrame({'lag': lags, 'acf_unadjusted': acf_by_chr(byc0, lags), 'acf_composition_adjusted': acf_by_chr(byc1, lags), 'null_p97.5': np.quantile(null, .975, 0)}).round(4).to_csv(tab('Table_S5d_acf_with_without_composition.csv'), index=False)

# ---- genomic coordinates (GRCh37) and liver TADs
G37 = pa.tables()['homo_sapiens-GRCh37-ensembl100']; G37 = G37[~G37.index.duplicated()][['Chromosome', 'Start', 'End', 'Strand']]; G37.columns = ['chr37', 'start', 'end', 'strand']
K = keep.set_index('gene_id').join(G37).join(de2[['stat', 'padj']]).join(Ts[['t_stage', 't_stage_comp_adjusted']]); K = K[K.chr37.notna()]; K['mid'] = (K.start + K.end) // 2
d = os.path.join(EXT, 'TAD-stability-heritability-master', 'data', '20binsTADlandscape', 'Liver_leung2015') + '/'
b6 = pd.read_csv(d + 'bin_6_Liver_leung2015.bed', sep='\t', header=None, names=['chr', 's', 'e']); b15 = pd.read_csv(d + 'bin_15_Liver_leung2015.bed', sep='\t', header=None, names=['chr', 's', 'e'])
b6['bin'] = b6.e - b6.s + 1; b6['exp_end'] = b6.s + 10 * b6.bin - 1; ends = {c: set(g.e) for c, g in b15.groupby('chr')}
b6['end'] = [next((e for e in range(r.exp_end - 3, r.exp_end + 4) if e in ends.get(r.chr, set())), np.nan) for _, r in b6.iterrows()]
TAD = b6.dropna(subset=['end']).rename(columns={'s': 'start'}); TAD['chr'] = TAD.chr.str.replace('chr', ''); TAD = TAD[TAD.chr.isin(CHR)].sort_values(['chr', 'start']).reset_index(drop=True); TAD['tad_id'] = np.arange(len(TAD))
TAD[['chr', 'start', 'end', 'tad_id']].to_csv(tab('Table_S5e_liver_TADs_reconstructed_hg19.csv'), index=False)
K['tad'] = -1
for c, g in TAD.groupby('chr'):
    idx = K.index[K.chr == c]; mm = K.loc[idx, 'mid'].values; st_ = g.start.values; en = g.end.values.astype(int); ids = g.tad_id.values
    j = np.searchsorted(st_, mm, side='right') - 1; ok = (j >= 0) & (mm <= en[np.clip(j, 0, len(en) - 1)]); K.loc[idx[ok], 'tad'] = ids[j[ok]]
pd.to_pickle(K, inter('K_tad.pkl'))

# ---- adjacent pairs: distance, orientation, family, TAD, co-expression
Am = A.values; gi = {g: i for i, g in enumerate(A.index)}; Ac = (Am - Am.mean(1, keepdims=True)) / (Am.std(1, keepdims=True) + 1e-9)
rows = []
for c, g in K.sort_values(['chr', 'grid_index']).groupby('chr'):
    g = g.reset_index(); a = g.iloc[:-1].reset_index(drop=True); b = g.iloc[1:].reset_index(drop=True)
    inter_bp = np.maximum(b.start.values - a.end.values, a.start.values - b.end.values)
    ori = np.where((a.strand == '-') & (b.strand == '+'), 'divergent', np.where((a.strand == '+') & (b.strand == '-'), 'convergent', 'tandem'))
    fam = [re.sub(r'[0-9]+[A-Z]?$', '', str(x)) == re.sub(r'[0-9]+[A-Z]?$', '', str(y)) and len(re.sub(r'[0-9]+[A-Z]?$', '', str(x))) >= 3 for x, y in zip(a.gene_name, b.gene_name)]
    co = [np.dot(Ac[gi[x]], Ac[gi[y]]) / Ac.shape[1] for x, y in zip(a.gene_id, b.gene_id)]
    rows.append(pd.DataFrame({'chr': c, 'g1': a.gene_id, 'g2': b.gene_id, 'n1': a.gene_name, 'n2': b.gene_name, 'orientation': ori, 'intergenic_bp': inter_bp, 'paralog_family': fam, 'tad1': a.tad, 'tad2': b.tad,
                              'stat1': a.stat, 'stat2': b.stat, 't1': a.t_stage, 't2': b.t_stage, 't1c': a.t_stage_comp_adjusted, 't2c': b.t_stage_comp_adjusted, 'coexpr': co}))
P = pd.concat(rows); P['same_tad'] = (P.tad1 >= 0) & (P.tad1 == P.tad2); P.round(4).to_csv(tab('Table_S5f_adjacent_pairs_distance_orientation_TAD.csv'), index=False); pd.to_pickle(P, inter('pairs.pkl'))
q = P[~P.paralog_family].copy(); q['dbin'] = pd.cut(q.intergenic_bp, [-1e9, 0, 1e3, 1e4, 5e4, 1e5, 5e5, 1e9], labels=['overlap', '<1kb', '1-10kb', '10-50kb', '50-100kb', '100-500kb', '>500kb'])
q.groupby('dbin', observed=True).apply(lambda h: pd.Series({'n': len(h), 'r_t': np.corrcoef(h.t1, h.t2)[0, 1], 'r_t_comp_adj': np.corrcoef(h.t1c, h.t2c)[0, 1], 'coexpr': h.coexpr.mean()})).round(4).to_csv(tab('Table_S5g_concordance_by_distance.csv'))
q[q.intergenic_bp < 5e4].groupby('orientation').apply(lambda h: pd.Series({'n': len(h), 'r_t': np.corrcoef(h.t1, h.t2)[0, 1], 'coexpr': h.coexpr.mean()})).round(4).to_csv(tab('Table_S5h_concordance_by_orientation.csv'))
# TAD test stratified by distance
qq = P[(P.tad1 >= 0) & (P.tad2 >= 0)].dropna(subset=['stat1', 'stat2']).copy(); qq['dbin'] = pd.cut(qq.intergenic_bp.clip(lower=0), [-1, 5e4, 1e5, 2e5, 5e5, 1e6, 5e6, 1e9])
def dif(lab): return np.corrcoef(qq[lab].stat1, qq[lab].stat2)[0, 1] - np.corrcoef(qq[~lab].stat1, qq[~lab].stat2)[0, 1]
obs = dif(qq.same_tad.values); nl = []
for _ in range(500):
    lab = qq.groupby('dbin', observed=True).same_tad.transform(lambda v: rng.permutation(v.values)).astype(bool).values; nl.append(dif(lab))
pd.DataFrame({'observed_delta_r': [obs], 'null_mean': [np.mean(nl)], 'null_sd': [np.std(nl)], 'p': [(np.sum(np.array(nl) >= obs) + 1) / 501]}).to_csv(tab('Table_S5i_TAD_test_stratified.csv'), index=False)

# ---- thresholds and bimodality
Dm = pd.get_dummies(C.loc[s, ['cohorte']], drop_first=True).astype(float); Dm['orden'] = C.loc[s, 'orden'].astype(float); Dm.insert(0, 'const', 1.0); rows = []
for var in list(MARKERS) + ['offset', 'n_peaks']:
    y = C.loc[s, var].values; m0 = sm.OLS(y, Dm.values).fit(); r = {'variable': var, 'AIC_linear': m0.aic}
    for br in (2, 3, 4, 5):
        Db = Dm.copy(); Db['step'] = (C.loc[s, 'orden'] >= br).astype(float).values; m1 = sm.OLS(y, Db.values).fit(); r[f'dAIC_step_at_{STAGES[br]}'] = m0.aic - m1.aic; r[f't_step_{STAGES[br]}'] = m1.tvalues[-1]
    rows.append(r)
pd.DataFrame(rows).round(3).to_csv(tab('Table_S7a_threshold_vs_linear.csv'), index=False)
rows = []
for var in ['HSC', 'Colangiocito', 'Hepatocito']:
    for st in STAGES:
        y = C[C.estadio == st][var].values.reshape(-1, 1)
        if len(y) < 15: continue
        g1 = GaussianMixture(1, random_state=0).fit(y); g2 = GaussianMixture(2, random_state=0, n_init=5).fit(y)
        rows.append({'variable': var, 'estadio': st, 'n': len(y), 'dBIC_1_minus_2': g1.bic(y) - g2.bic(y), 'means_2comp': np.round(np.sort(g2.means_.ravel()), 2).tolist(), 'min_weight': round(g2.weights_.min(), 2)})
pd.DataFrame(rows).to_csv(tab('Table_S7b_within_stage_bimodality.csv'), index=False)
print('done')

# ---- containment of contiguous DE neighbourhoods within a single liver TAD, against random windows of the same size
#      on the same chromosome (Supplementary Table S5k; Supplementary Fig. S4D)
Vk = pd.read_csv(tab('Table_S5c_contiguous_DE_neighbourhoods.csv'))
Kt2 = pd.read_pickle(inter('K_tad.pkl')).sort_values(['chr', 'grid_index']); Kt2['gene_name'] = Kt2.gene_name.astype(str)
tad_of = dict(zip(Kt2.gene_name, Kt2.tad)); byc_tad = {c: g.tad.values for c, g in Kt2.groupby('chr')}
def within_one(tads): tads = [t for t in tads]; return all(pd.notna(t) for t in tads) and len(set(tads)) == 1
obs_rows = []
for _, v in Vk.iterrows():
    gs_ = str(v.genes).split(','); tt = [tad_of.get(g, np.nan) for g in gs_]
    if all(pd.isna(t) for t in tt): continue
    obs_rows.append({'chr': str(v.chr), 'n_genes': len(gs_), 'within': within_one(tt)})
OB = pd.DataFrame(obs_rows); rng_k = np.random.default_rng(0); rows_k = []
for nmin in (3, 4, 5):
    sub = OB[OB.n_genes >= nmin]; perm = []
    for _ in range(500):
        hits = 0; tot = 0
        for _, r in sub.iterrows():
            arr = byc_tad.get(r.chr)
            if arr is None or len(arr) <= r.n_genes: continue
            st_ = rng_k.integers(0, len(arr) - r.n_genes); w_ = arr[st_:st_ + r.n_genes]; tot += 1; hits += within_one(list(w_))
        perm.append(hits / max(tot, 1))
    perm = np.array(perm)
    rows_k.append({'min_genes': nmin, 'n_neighbourhoods': len(sub), 'observed_within_one_TAD': sub.within.mean(), 'chance_mean': perm.mean(), 'chance_sd': perm.std(),
                   'P_greater': (np.sum(perm >= sub.within.mean()) + 1) / (len(perm) + 1)})
pd.DataFrame(rows_k).round(4).to_csv(tab('Table_S5k_neighbourhoods_within_one_TAD.csv'), index=False); print(pd.DataFrame(rows_k).round(3).to_string(index=False))
