"""14 — Shared significant eQTL variants between adjacent genes (GTEx v8 liver, all significant variant–gene pairs).
For adjacent non-paralogous pairs in which both genes are liver eGenes: (i) do they share a significant variant?
(ii) same-direction effect at shared variants (same allele raises both or lowers both)? Compared between pairs coupled
in the fibrosis stage effect and uncoupled pairs, stratified by intergenic distance (Mantel–Haenszel), plus logistic
models. Outputs: Table_S13d–g
"""
import os, numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats
from statsmodels.stats.contingency_tables import StratifiedTable
from common import *

SV = pd.read_csv(os.path.join(EXT, 'Liver.v8.signif_variant_gene_pairs.txt.gz'), sep='\t', usecols=['variant_id', 'gene_id', 'pval_nominal', 'slope'])
SV['gene'] = SV.gene_id.str.split('.').str[0]
EG = pd.read_csv(os.path.join(EXT, 'Liver.v8.egenes.txt.gz'), sep='\t', usecols=['gene_id', 'qval']); EG['gene'] = EG.gene_id.str.split('.').str[0]
egenes = set(EG.gene[EG.qval <= 0.05]); print('liver eGenes (q<=0.05): %d; significant variant-gene pairs: %d' % (len(egenes), len(SV)))
var = SV.groupby('gene').variant_id.apply(set).to_dict(); slope = SV.set_index(['gene', 'variant_id']).slope
lead = SV.loc[SV.groupby('gene').pval_nominal.idxmin()].set_index('gene').variant_id.to_dict()

P = pd.read_csv(tab('Table_S5f_adjacent_pairs_distance_orientation_TAD.csv')); P = P[~P.paralog_family].dropna(subset=['t1', 't2']).copy()
P['both_egenes'] = P.g1.isin(egenes) & P.g2.isin(egenes)
def info(a, b):
    s = var.get(a, set()) & var.get(b, set())
    if not s: return 0, np.nan, False
    sl = np.array([slope.loc[(a, v)] * slope.loc[(b, v)] for v in s]); return len(s), float(np.mean(sl > 0)), (lead.get(a) in var.get(b, set())) or (lead.get(b) in var.get(a, set()))
res = [info(a, b) if e else (0, np.nan, False) for a, b, e in zip(P.g1, P.g2, P.both_egenes)]
P['n_shared_sig_var'] = [r[0] for r in res]; P['frac_same_direction'] = [r[1] for r in res]; P['lead_shared'] = [r[2] for r in res]
P['shares'] = P.n_shared_sig_var > 0; P['shares_same_dir'] = P.shares & (P.frac_same_direction >= 0.5)
P['coupled'] = (np.sign(P.t1) == np.sign(P.t2)) & (P.t1.abs() > 3) & (P.t2.abs() > 3)
P['coupled_comp'] = (np.sign(P.t1c) == np.sign(P.t2c)) & (P.t1c.abs() > 2) & (P.t2c.abs() > 2)
P['uncoupled'] = (P.t1.abs() < 1) | (P.t2.abs() < 1) | (np.sign(P.t1) != np.sign(P.t2))
P['dbin'] = pd.cut(P.intergenic_bp, [-1e9, 0, 1e3, 1e4, 5e4, 1e5, 5e5, 1e9], labels=['overlap', '<1kb', '1-10kb', '10-50kb', '50-100kb', '100-500kb', '>500kb'])
E = P[P.both_egenes].copy(); print('adjacent non-paralog pairs with both genes eGenes: %d (of %d); sharing >=1 significant variant: %d' % (len(E), len(P), int(E.shares.sum())))
P.round(4).to_csv(tab('Table_S13d_adjacent_pairs_shared_GTEx_signif_variants.csv'), index=False)

rows = []
for clab, ccol in [('coupled (stage, |t|>3)', 'coupled'), ('coupled after composition (|t|>2)', 'coupled_comp')]:
    for olab, ocol in [('share >=1 significant variant', 'shares'), ('share variant with same-direction effect', 'shares_same_dir'), ('lead variant of one is significant for the other', 'lead_shared')]:
        sub = E[E[ccol] | E.uncoupled].copy(); sub['unc'] = sub.uncoupled & ~sub[ccol]; tabs = []
        for b, g in sub.groupby('dbin', observed=True):
            t_ = np.array([[(g[ccol] & g[ocol]).sum(), (g[ccol] & ~g[ocol]).sum()], [(g.unc & g[ocol]).sum(), (g.unc & ~g[ocol]).sum()]])
            if (t_.sum(1) > 0).all(): tabs.append(t_)
        st = StratifiedTable(tabs); ci = st.oddsratio_pooled_confint()
        rows.append({'coupling': clab, 'outcome': olab, 'n_coupled': int(sub[ccol].sum()), 'frac_coupled': sub.loc[sub[ccol], ocol].mean(), 'n_uncoupled': int(sub.unc.sum()), 'frac_uncoupled': sub.loc[sub.unc, ocol].mean(),
                     'MH_OR_distance_stratified': st.oddsratio_pooled, 'CI_low': ci[0], 'CI_high': ci[1], 'MH_p': st.test_null_odds().pvalue})
R = pd.DataFrame(rows); R.round(4).to_csv(tab('Table_S13e_shared_GTEx_eQTL_coupled_vs_uncoupled.csv'), index=False); print(R.round(3).to_string(index=False))

# continuous: concordance of stage effect by sharing, within distance bins
rows = []
for b, g in E.groupby('dbin', observed=True):
    for s_, h in g.groupby('shares_same_dir'):
        if len(h) >= 25: rows.append({'distance': str(b), 'shares_same_dir': s_, 'n': len(h), 'r_stage_effect': np.corrcoef(h.t1, h.t2)[0, 1], 'r_comp_adjusted': np.corrcoef(h.t1c, h.t2c)[0, 1], 'coexpr': h.coexpr.mean()})
C = pd.DataFrame(rows); C.round(4).to_csv(tab('Table_S13f_concordance_by_shared_eQTL_distance.csv'), index=False); print(C.round(3).to_string(index=False))
# logistic / linear models across all testable pairs
E['concordant'] = (np.sign(E.t1) == np.sign(E.t2)).astype(float); E['logd'] = np.log10(E.intergenic_bp.clip(lower=0) + 1e3)
out = []
for ocol in ['shares', 'shares_same_dir', 'lead_shared']:
    m = sm.Logit(E.concordant, sm.add_constant(E[['logd']].assign(x=E[ocol].astype(float)))).fit(disp=0)
    ml = sm.OLS(E.coexpr, sm.add_constant(E[['logd']].assign(x=E[ocol].astype(float)))).fit()
    out.append({'predictor': ocol, 'OR_concordant_sign': np.exp(m.params['x']), 'CI_low': np.exp(m.conf_int().loc['x', 0]), 'CI_high': np.exp(m.conf_int().loc['x', 1]), 'p': m.pvalues['x'], 'delta_coexpr': ml.params['x'], 'p_coexpr': ml.pvalues['x']})
M = pd.DataFrame(out); M.round(4).to_csv(tab('Table_S13g_models_concordance_on_shared_eQTL.csv'), index=False); print(M.round(4).to_string(index=False))
