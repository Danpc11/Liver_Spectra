"""07 — GSEA pre-ranked (Hallmark) on ordinal and composition-adjusted statistics; DoRothEA TF activity vs stage;
shared regulators between coupled neighbour pairs; TF x Hallmark overlap. Outputs: Table_S10_*
"""
import os, numpy as np, pandas as pd, gseapy as gp, pyreadr
from scipy import stats
from statsmodels.stats.multitest import multipletests
from collections import Counter
from common import *

X, keep = pd.read_pickle(inter('expr.pkl')); A = pd.read_pickle(inter('expr_adj.pkl')); M = pd.read_pickle(inter('meta.pkl')); de, de2 = pd.read_pickle(inter('de.pkl')); comp = pd.read_pickle(inter('comp.pkl'))
T = pd.read_csv(tab('Table_S6c_stage_effect_per_gene_with_without_composition.csv')).set_index('gene_id'); P = pd.read_pickle(inter('pairs.pkl'))
names = keep.set_index('gene_id').gene_name.astype(str); GMT = os.path.join(EXT, 'hallmark.gmt')

def prerank(stat, out):
    rnk = pd.Series(stat.values, index=names.reindex(stat.index).values).dropna().groupby(level=0).mean().sort_values(ascending=False)
    pre = gp.prerank(rnk=rnk, gene_sets=GMT, permutation_num=1000, seed=0, threads=1, min_size=10, max_size=500, outdir=None, verbose=False)
    G = pre.res2d.copy(); G['Term'] = G.Term.str.replace('HALLMARK_', ''); G = G[['Term', 'NES', 'NOM p-val', 'FDR q-val', 'Lead_genes']].sort_values('NES'); G.to_csv(tab(out), index=False); return G
G = prerank(de2.stat, 'Table_S10a_GSEA_hallmark_ordinal.csv'); G2 = prerank(T.t_stage_comp_adjusted, 'Table_S10b_GSEA_hallmark_composition_adjusted.csv')
print('Hallmark sets FDR<0.05: ordinal', (G['FDR q-val'] < .05).sum(), '| composition-adjusted', (G2['FDR q-val'] < .05).sum())

# DoRothEA (A-C)
D = pyreadr.read_r(os.path.join(EXT, 'dorothea_hs.rda')); D = D[list(D)[0]]; D = D[D.confidence.isin(['A', 'B', 'C'])]; D.to_csv(tab('Table_S10c_dorothea_ABC_regulons_used.csv'), index=False)
Az = ((A.T - A.mean(axis=1)) / A.std(axis=1).replace(0, 1)).T; Az.index = names.reindex(Az.index).values; Az = Az[~Az.index.duplicated()]
act, nt = {}, {}
for tf, g in D.groupby('tf'):
    g = g[g.target.isin(Az.index)]
    if len(g) < 10: continue
    w = g.set_index('target').mor; act[tf] = Az.loc[w.index].mul(w, axis=0).sum() / np.sqrt((w ** 2).sum()); nt[tf] = len(g)
ACT = pd.DataFrame(act); s = M.index[M.estadio != 'Control']
t0, df0 = ols_t(ACT.loc[s].values, design(M, s)); t1, df1 = ols_t(ACT.loc[s].values, design(M, s, extra=comp))
TF = pd.DataFrame({'tf': ACT.columns, 'n_targets': [nt[t] for t in ACT.columns], 't_stage': t0, 't_stage_comp_adjusted': t1})
TF['q'] = multipletests(2 * stats.t.sf(np.abs(t0), df0), method='fdr_bh')[1]; TF['q_adj_comp'] = multipletests(2 * stats.t.sf(np.abs(t1), df1), method='fdr_bh')[1]
TF.sort_values('t_stage').round(4).to_csv(tab('Table_S10d_DoRothEA_TF_activity_vs_stage.csv'), index=False); ACT.to_csv(tab('Table_S10e_TF_activity_per_sample.csv'))
print('TFs q<0.05:', (TF.q < .05).sum(), '| after composition:', (TF.q_adj_comp < .05).sum())

# shared regulators in coupled vs uncoupled neighbour pairs
reg = D.groupby('target').tf.apply(set).to_dict(); Pq = P[~P.paralog_family].copy()
Pq['shared_tf'] = [len(reg.get(a, set()) & reg.get(b, set())) for a, b in zip(Pq.n1, Pq.n2)]
Pq['coupled'] = (np.sign(Pq.t1) == np.sign(Pq.t2)) & (Pq.t1.abs() > 3) & (Pq.t2.abs() > 3); Pq['uncoupled'] = (Pq.t1.abs() < 1) | (Pq.t2.abs() < 1) | (np.sign(Pq.t1) != np.sign(Pq.t2))
a = Pq[Pq.coupled]; u = Pq[Pq.uncoupled]; odds, p = stats.fisher_exact([[(a.shared_tf > 0).sum(), (a.shared_tf == 0).sum()], [(u.shared_tf > 0).sum(), (u.shared_tf == 0).sum()]])
pd.DataFrame({'frac_shared_coupled': [(a.shared_tf > 0).mean()], 'frac_shared_uncoupled': [(u.shared_tf > 0).mean()], 'odds': [odds], 'p': [p]}).to_csv(tab('Table_S10f_shared_TF_coupled_pairs.csv'), index=False)
cu, cd = Counter(), Counter()
for _, r in a.iterrows(): (cu if r.t1 > 0 else cd).update(reg.get(r.n1, set()) & reg.get(r.n2, set()))
pd.DataFrame({'tf': list(cu) + list(cd), 'direction': ['up'] * len(cu) + ['down'] * len(cd), 'n_pairs': list(cu.values()) + list(cd.values())}).to_csv(tab('Table_S10g_shared_TFs_by_direction.csv'), index=False)

# TF x Hallmark overlap matrix (intrinsic TFs x enriched pathways)
sig = TF[TF.q_adj_comp < .05]; tfs = sig.reindex(sig.t_stage_comp_adjusted.abs().sort_values(ascending=False).index).tf.head(16).tolist()
Gs = G[G['FDR q-val'] < .05]; paths = Gs.reindex(Gs.NES.abs().sort_values(ascending=False).index).Term.head(12).tolist()
regs = {t: set(D[D.tf == t].target) for t in tfs}; lead = {r.Term: set(str(r.Lead_genes).split(';')) for r in G.itertuples()}
Mx = pd.DataFrame({t: [len(regs[t] & lead[p]) for p in paths] for t in tfs}, index=paths); Mx = Mx.loc[Mx.sum(1) >= 3, Mx.sum(0) >= 5]; Mx.to_csv(tab('Table_S10h_TF_x_Hallmark_overlap.csv'))
print('done')
