"""12 — Do coupled neighbour pairs share causal eQTL variants? Liver fine-mapped credible sets (eQTL Catalogue,
GTEx liver, QTD000266, SuSiE). For adjacent non-paralogous gene pairs, test whether their credible sets share a variant,
comparing pairs coupled in the fibrosis stage effect with uncoupled pairs at matched intergenic distance.
Outputs: Table_S13a–c
"""
import os, numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.contingency_tables import StratifiedTable
import statsmodels.api as sm
from common import *

CS = pd.read_csv(os.path.join(EXT, 'QTD000266.credible_sets.tsv.gz'), sep='\t')
print('credible sets: %d variants, %d genes, %d sets' % (len(CS), CS.gene_id.nunique(), CS.cs_id.nunique()))
var = CS.groupby('gene_id').variant.apply(set).to_dict()
# variants weighted by PIP: maximum shared PIP product (probability both genes use the same causal variant)
pip = CS.groupby(['gene_id', 'variant']).pip.max()
pipd = {g: d.droplevel(0).to_dict() for g, d in pip.groupby(level=0)}

P = pd.read_csv(tab('Table_S5f_adjacent_pairs_distance_orientation_TAD.csv')); P = P[~P.paralog_family].dropna(subset=['t1', 't2']).copy()
P['both_egenes'] = P.g1.isin(var) & P.g2.isin(var)
def shared(a, b):
    s = var.get(a, set()) & var.get(b, set())
    return len(s), (max(pipd[a][v] * pipd[b][v] for v in s) if s else 0.0)
sh = [shared(a, b) for a, b in zip(P.g1, P.g2)]; P['n_shared_var'] = [x[0] for x in sh]; P['max_shared_pip'] = [x[1] for x in sh]
P['shares_eqtl'] = P.n_shared_var > 0; P['shares_eqtl_pip01'] = P.max_shared_pip >= 0.01
P['coupled'] = (np.sign(P.t1) == np.sign(P.t2)) & (P.t1.abs() > 3) & (P.t2.abs() > 3)
P['coupled_comp'] = (np.sign(P.t1c) == np.sign(P.t2c)) & (P.t1c.abs() > 2) & (P.t2c.abs() > 2)
P['uncoupled'] = (P.t1.abs() < 1) | (P.t2.abs() < 1) | (np.sign(P.t1) != np.sign(P.t2))
P['dbin'] = pd.cut(P.intergenic_bp, [-1e9, 0, 1e3, 1e4, 5e4, 1e5, 5e5, 1e9], labels=['overlap', '<1kb', '1-10kb', '10-50kb', '50-100kb', '100-500kb', '>500kb'])
E = P[P.both_egenes].copy(); print('adjacent non-paralog pairs: %d; both genes with liver eQTL credible sets: %d' % (len(P), len(E)))
P.round(4).to_csv(tab('Table_S13a_adjacent_pairs_shared_eQTL.csv'), index=False)

rows = []
for lab, col in [('coupled (stage, |t|>3)', 'coupled'), ('coupled after composition (|t|>2)', 'coupled_comp')]:
    for outcome in ['shares_eqtl', 'shares_eqtl_pip01']:
        sub = E[E[col] | E.uncoupled]; tabs = []
        for b, g in sub.groupby('dbin', observed=True):
            t_ = np.array([[(g[col] & g[outcome]).sum(), (g[col] & ~g[outcome]).sum()], [(g.uncoupled & ~g[col] & g[outcome]).sum(), (g.uncoupled & ~g[col] & ~g[outcome]).sum()]])
            if (t_.sum(1) > 0).all(): tabs.append(t_)
        st = StratifiedTable(tabs); test = st.test_null_odds()
        rows.append({'coupling_definition': lab, 'eqtl_definition': 'any shared credible variant' if outcome == 'shares_eqtl' else 'shared variant PIP product >= 0.01',
                     'n_coupled': int(sub[col].sum()), 'frac_shared_coupled': sub.loc[sub[col], outcome].mean(), 'n_uncoupled': int((sub.uncoupled & ~sub[col]).sum()), 'frac_shared_uncoupled': sub.loc[sub.uncoupled & ~sub[col], outcome].mean(),
                     'MH_odds_ratio_distance_stratified': st.oddsratio_pooled, 'MH_CI_low': st.oddsratio_pooled_confint()[0], 'MH_CI_high': st.oddsratio_pooled_confint()[1], 'MH_p': test.pvalue})
R = pd.DataFrame(rows); R.round(4).to_csv(tab('Table_S13b_shared_eQTL_coupled_vs_uncoupled.csv'), index=False); print(R.round(3).to_string(index=False))

# continuous: does sharing an eQTL predict concordance of the stage effect? r(t1,t2) by distance bin and eQTL sharing
rows = []
for b, g in E.groupby('dbin', observed=True):
    for s_, h in g.groupby('shares_eqtl'):
        if len(h) >= 20: rows.append({'distance': str(b), 'shares_eqtl': s_, 'n': len(h), 'r_stage_effect': np.corrcoef(h.t1, h.t2)[0, 1], 'r_comp_adjusted': np.corrcoef(h.t1c, h.t2c)[0, 1], 'coexpr': h.coexpr.mean()})
C = pd.DataFrame(rows); C.round(4).to_csv(tab('Table_S13c_concordance_by_shared_eQTL_and_distance.csv'), index=False); print(C.round(3).to_string(index=False))
# logistic model: concordant sign ~ shares_eqtl + log distance
E['concordant'] = (np.sign(E.t1) == np.sign(E.t2)).astype(float); E['logd'] = np.log10(E.intergenic_bp.clip(lower=0) + 1e3)
m = sm.Logit(E.concordant, sm.add_constant(E[['logd']].assign(shares=E.shares_eqtl.astype(float)))).fit(disp=0)
print('logit concordant sign ~ log distance + shared eQTL: OR(shared) = %.2f (95%% CI %.2f-%.2f), P = %.2g' % (np.exp(m.params['shares']), *np.exp(m.conf_int().loc['shares']), m.pvalues['shares']))
