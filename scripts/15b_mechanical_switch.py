"""15b — The F3→F4 mechanical switch, gene by gene.
(1) Per-gene step at F4 on batch-corrected expression (linear trend + step indicator; cohort, sex), with and without
    the six composition scores; DESeq2 F4 vs F3 with and without composition covariates.
(2) GSEA pre-ranked on the composition-adjusted switch statistic: Hallmark + curated mechanotransduction sets.
(3) Drivers: step at F4 in DoRothEA TF activity (TEAD, SRF, NF-kB, SMAD, AP-1…), composition-adjusted.
(4) Positional organisation of the switch: spatial autocorrelation along the genome, power by period band against a
    gene-order permutation null, and enrichment of switch genes in contiguous DE neighbourhoods and coupled pairs.
(5) Cell of origin of switch genes in scRNA-seq (cirrhosis vs healthy t within each population).
Outputs: Table_S15, Table_S15b, Table_S7c, Table_S17a–g
"""
import os, numpy as np, pandas as pd, gseapy as gp
from scipy import stats
from statsmodels.stats.multitest import multipletests
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from common import *

X, keep = pd.read_pickle(inter('expr.pkl')); A = pd.read_pickle(inter('expr_adj.pkl')); M = pd.read_pickle(inter('meta.pkl'))
comp = pd.read_pickle(inter('comp.pkl')); counts = pd.read_pickle(inter('counts.pkl'))
names = keep.set_index('gene_id').gene_name.astype(str)
s = M.index[M.estadio != 'Control']

# ---------- (0) programme scores per biopsy, threshold tests, and F4-likeness of composition (Supplementary Tables S15, S7c)
import statsmodels.api as sm
from sklearn.mixture import GaussianMixture
SETS15 = {'YAP_TAZ_targets': ['CCN1', 'CCN2', 'ANKRD1', 'AMOTL2', 'AXL', 'LATS2', 'TEAD1', 'TEAD4', 'CRIM1', 'F3', 'GADD45B', 'SERPINE1', 'TGFB2', 'THBS1', 'DKK1', 'WWC1', 'AJUBA', 'CYR61', 'CTGF'],
          'matrix_crosslinking': ['LOX', 'LOXL1', 'LOXL2', 'LOXL3', 'LOXL4', 'TGM2', 'PLOD1', 'PLOD2', 'PLOD3', 'P4HA1', 'P4HA2', 'PXDN', 'ELN', 'FBLN5', 'COL1A1', 'COL1A2', 'COL3A1'],
          'mechanosensing': ['PIEZO1', 'PIEZO2', 'ITGB1', 'ITGA5', 'ITGAV', 'ITGB5', 'PTK2', 'RHOA', 'ROCK1', 'ROCK2', 'YAP1', 'WWTR1', 'VCL', 'TLN1', 'FLNA', 'ACTN1', 'MYL9', 'TAGLN', 'CNN1', 'ACTA2'],
          'TGFb_SMAD': ['TGFB1', 'TGFB2', 'TGFB3', 'TGFBR1', 'TGFBR2', 'SMAD2', 'SMAD3', 'SMAD7', 'SERPINE1', 'SKIL', 'ID1', 'JUNB', 'PMEPA1', 'LTBP1', 'LTBP2', 'FN1', 'POSTN'],
          'senescence_SASP': ['CDKN1A', 'CDKN2A', 'CDKN2B', 'SERPINE1', 'IL6', 'IL1A', 'IL1B', 'CXCL8', 'CCL2', 'MMP3', 'IGFBP3', 'IGFBP7', 'GDF15', 'TP53']}
Azs = ((A.T - A.mean(axis=1)) / A.std(axis=1).replace(0, 1)).T; Azs.index = names.reindex(Azs.index).values
SC15 = pd.DataFrame({k: Azs.loc[[g for g in v if g in Azs.index]].mean() for k, v in SETS15.items()}).join(M[['cohorte', 'estadio', 'orden', 'nas']]).join(comp)
SC15 = SC15[SC15.estadio != 'Control']; SC15.round(4).to_csv(tab('Table_S15b_mechanics_scores_per_sample.csv'))
D0 = pd.get_dummies(SC15[['cohorte']], drop_first=True).astype(float); D0['orden'] = SC15.orden.astype(float); D0.insert(0, 'const', 1.0); D1 = D0.join(SC15[list(comp.columns)])
rows15 = []
for k in SETS15:
    y = SC15[k].values; m0 = sm.OLS(y, D0.values).fit(); m1 = sm.OLS(y, D1.values).fit(); j = list(D0.columns).index('orden')
    r = {'score': k, 't_stage': m0.tvalues[j], 't_stage_comp_adj': m1.tvalues[j]}
    for br in (3, 4, 5):
        Db = D0.copy(); Db['step'] = (SC15.orden >= br).astype(float).values; mb = sm.OLS(y, Db.values).fit(); r[f'dAIC_step_{STAGES[br]}'] = m0.aic - mb.aic; r[f't_step_{STAGES[br]}'] = mb.tvalues[-1]
        Dbc = D1.copy(); Dbc['step'] = (SC15.orden >= br).astype(float).values; r[f't_step_{STAGES[br]}_comp_adj'] = sm.OLS(y, Dbc.values).fit().tvalues[-1]
    y3 = SC15[SC15.estadio == 'F3'][k].values.reshape(-1, 1); g1 = GaussianMixture(1, random_state=0).fit(y3); g2 = GaussianMixture(2, random_state=0, n_init=5).fit(y3)
    r['dBIC_F3_bimodal'] = g1.bic(y3) - g2.bic(y3); r['F3_minor_weight'] = g2.weights_.min(); rows15.append(r)
pd.DataFrame(rows15).round(4).to_csv(tab('Table_S15_mechanics_scores_thresholds.csv'), index=False)
cc = list(comp.columns); w = SC15[SC15.estadio == 'F4'][cc].mean() - SC15[SC15.estadio == 'F2'][cc].mean()
Zc = (SC15[cc] - SC15[cc].mean()) / SC15[cc].std(); SC15['F4_likeness'] = Zc @ (w / np.linalg.norm(w))
SC15[['estadio', 'cohorte', 'F4_likeness', 'Hepatocito', 'HSC', 'Colangiocito', 'YAP_TAZ_targets', 'mechanosensing', 'senescence_SASP', 'matrix_crosslinking', 'nas']].round(4).to_csv(tab('Table_S7c_F4_likeness_per_sample.csv'))
thr = SC15[SC15.estadio == 'F4'].F4_likeness.quantile(0.25); f3_ = SC15[SC15.estadio == 'F3']
print('F3 with F4-like composition: %.0f%%; YAP/TAZ F4-like F3 %.2f vs F4 %.2f' % (100 * (f3_.F4_likeness > thr).mean(), f3_[f3_.F4_likeness > thr].YAP_TAZ_targets.mean(), SC15[SC15.estadio == 'F4'].YAP_TAZ_targets.mean()))

MECH = {
    'YAP/TAZ targets': ['CCN1', 'CCN2', 'ANKRD1', 'AMOTL2', 'AXL', 'LATS2', 'TEAD1', 'TEAD4', 'CRIM1', 'F3', 'GADD45B', 'SERPINE1', 'TGFB2', 'THBS1', 'DKK1', 'AJUBA', 'WWC1', 'NUAK2', 'FSTL1', 'CYR61', 'CTGF'],
    'Mechanosensing': ['PIEZO1', 'PIEZO2', 'ITGB1', 'ITGA5', 'ITGAV', 'ITGB5', 'PTK2', 'RHOA', 'ROCK1', 'ROCK2', 'YAP1', 'WWTR1', 'VCL', 'TLN1', 'FLNA', 'ACTN1', 'MYL9', 'TAGLN', 'CNN1', 'ACTA2', 'ZYX', 'TRPV4', 'LMNA'],
    'Matrix crosslinking': ['LOX', 'LOXL1', 'LOXL2', 'LOXL3', 'LOXL4', 'TGM2', 'PLOD1', 'PLOD2', 'PLOD3', 'P4HA1', 'P4HA2', 'PXDN', 'ELN', 'FBLN5', 'COL1A1', 'COL1A2', 'COL3A1'],
    'MRTF–SRF targets': ['ACTA2', 'TAGLN', 'CNN1', 'MYL9', 'VCL', 'CTGF', 'CCN2', 'CYR61', 'CCN1', 'FHL2', 'ACTB', 'ACTG1', 'MYH9', 'TPM1', 'CALD1', 'SRF', 'FOS', 'EGR1'],
    'Senescence/SASP': ['CDKN1A', 'CDKN2A', 'CDKN2B', 'SERPINE1', 'IL6', 'IL1A', 'IL1B', 'CXCL8', 'CCL2', 'MMP3', 'IGFBP3', 'IGFBP7', 'GDF15', 'TP53'],
}

# ---------- (1) per-gene step at F4 (OLS on batch-corrected log expression)
def step_t(Y, extra=None, br=5):
    D = design(M, s, extra=extra); D['step'] = (M.loc[s, 'orden'] >= br).astype(float).values
    Dv = D.values; XtXi = np.linalg.pinv(Dv.T @ Dv); B = XtXi @ Dv.T @ Y; R = Y - Dv @ B
    df = Y.shape[0] - np.linalg.matrix_rank(Dv); se = np.sqrt((R ** 2).sum(0) / df * np.diag(XtXi)[:, None])
    cols = list(D.columns); return B[cols.index('step')] / se[cols.index('step')], B[cols.index('orden')] / se[cols.index('orden')], df
Y = A[s].values.T
t0, tr0, df0 = step_t(Y); t1, tr1, df1 = step_t(Y, extra=comp)
SW = pd.DataFrame({'gene_name': names.reindex(A.index).values, 't_step_F4': t0, 't_trend': tr0, 't_step_F4_comp_adj': t1, 't_trend_comp_adj': tr1}, index=A.index)
SW['p_adj_step'] = 2 * stats.t.sf(np.abs(SW.t_step_F4_comp_adj), df1); SW['q_adj_step'] = multipletests(SW.p_adj_step, method='fdr_bh')[1]

# DESeq2 F4 vs F3, without and with composition covariates
ss = M.index[M.estadio.isin(['F3', 'F4'])]
cnt = pd.concat([counts[c].loc[keep.gene_id] for c in COHORTS], axis=1)[ss].T.astype(int)
md = M.loc[ss, ['cohorte']].copy(); md['cond'] = M.loc[ss, 'estadio'].values
res = {}
for lab_, extra in [('raw', None), ('comp_adj', comp.loc[ss])]:
    mdx = md.copy(); design_ = '~ cohorte + cond'
    if extra is not None:
        for c in extra.columns: mdx[c] = extra[c].values
        design_ = '~ cohorte + ' + ' + '.join(extra.columns) + ' + cond'
    dds = DeseqDataSet(counts=cnt, metadata=mdx, design=design_, quiet=True); dds.deseq2()
    st = DeseqStats(dds, contrast=['cond', 'F4', 'F3'], quiet=True); st.summary(); r = st.results_df; res[lab_] = r
    SW[f'log2FC_F4_vs_F3_{lab_}'] = r.log2FoldChange.reindex(SW.index); SW[f'padj_F4_vs_F3_{lab_}'] = r.padj.reindex(SW.index)
# primary definition: composition-adjusted step at F4 across all 427 biopsies (FDR < 0.05); the DESeq2 F4-vs-F3 contrast
# (104 biopsies, six extra covariates) is underpowered and is reported as supporting evidence only
SW['switch_gene'] = (SW.q_adj_step < 0.05) & (SW.t_step_F4_comp_adj > 0)
SW['switch_gene_down'] = (SW.q_adj_step < 0.05) & (SW.t_step_F4_comp_adj < 0)
SW['switch_supported_DESeq2_raw'] = SW.switch_gene & (SW.padj_F4_vs_F3_raw < 0.05) & (SW.log2FC_F4_vs_F3_raw > 0)
for k, g in MECH.items(): SW[k] = SW.gene_name.isin(g)
SW.round(5).to_csv(tab('Table_S17a_F4_switch_per_gene.csv'))
print('DE F4 vs F3 padj<0.05: raw %d, composition-adjusted %d' % ((res['raw'].padj < .05).sum(), (res['comp_adj'].padj < .05).sum()))
print('switch genes (composition-adjusted step at F4, FDR<0.05): up %d, down %d; up and DE in DESeq2 F4 vs F3: %d' % (SW.switch_gene.sum(), SW.switch_gene_down.sum(), SW.switch_supported_DESeq2_raw.sum()))
for k in MECH: print('  %-22s switch-up %d / %d' % (k, (SW.switch_gene & SW[k]).sum(), SW[k].sum()))

# ---------- (2) GSEA on the composition-adjusted switch statistic
hall = {}
for line in open(os.path.join(EXT, 'hallmark.gmt')):
    p = line.rstrip('\n').split('\t'); hall[p[0].replace('HALLMARK_', '')] = p[2:]
gsets = {**hall, **{('MECH: ' + k): v for k, v in MECH.items()}}
rnk = SW.dropna(subset=['gene_name']).groupby('gene_name').t_step_F4_comp_adj.mean().sort_values(ascending=False)
pre = gp.prerank(rnk=rnk, gene_sets=gsets, permutation_num=1000, seed=0, threads=1, min_size=8, max_size=500, outdir=None, verbose=False)
GS = pre.res2d[['Term', 'NES', 'NOM p-val', 'FDR q-val', 'Lead_genes']].sort_values('NES', ascending=False); GS.to_csv(tab('Table_S17b_GSEA_F4_switch.csv'), index=False)
print('\nGSEA on composition-adjusted F4 switch (FDR<0.05):'); print(GS[GS['FDR q-val'] < 0.05][['Term', 'NES', 'FDR q-val']].round(3).to_string(index=False))

# ---------- (3) drivers: step at F4 in TF activity
ACT = pd.read_csv(tab('Table_S10e_TF_activity_per_sample.csv'), index_col=0)
ta0, _, dfa = step_t(ACT.loc[s].values); ta1, _, dfa1 = step_t(ACT.loc[s].values, extra=comp)
TFS = pd.DataFrame({'tf': ACT.columns, 't_step_F4': ta0, 't_step_F4_comp_adj': ta1})
TFS['q_adj'] = multipletests(2 * stats.t.sf(np.abs(TFS.t_step_F4_comp_adj), dfa1), method='fdr_bh')[1]
TFS = TFS.sort_values('t_step_F4_comp_adj', ascending=False); TFS.round(4).to_csv(tab('Table_S17c_TF_activity_F4_switch.csv'), index=False)
print('\nTF activity switching at F4 (composition-adjusted q<0.05):', (TFS.q_adj < .05).sum()); print(TFS.head(15).round(2).to_string(index=False))
mech_tfs = TFS[TFS.tf.isin(['TEAD1', 'TEAD2', 'TEAD4', 'SRF', 'NFKB1', 'RELA', 'SMAD3', 'SMAD4', 'JUN', 'FOS', 'TP53', 'HIF1A', 'EGR1', 'KLF4'])]
print(mech_tfs.round(2).to_string(index=False))

# ---------- (4) positional organisation of the switch
order = keep.sort_values(['chr', 'grid_index']).copy(); order['t'] = SW.t_step_F4_comp_adj.reindex(order.gene_id).fillna(0).values
byc = [order[order.chr == c].t.values for c in CHR]; rng = np.random.default_rng(0); lags = (1, 2, 3, 5, 10)
obs = acf_by_chr(byc, lags); nul = np.array([acf_by_chr([rng.permutation(v) for v in byc], lags) for _ in range(200)])
ACFd = pd.DataFrame({'lag': lags, 'acf': obs, 'null_p2.5': np.quantile(nul, .025, 0), 'null_p97.5': np.quantile(nul, .975, 0), 'z': (obs - nul.mean(0)) / nul.std(0)})
ACFd.round(4).to_csv(tab('Table_S17d_F4_switch_spatial_autocorrelation.csv'), index=False); print('\nACF of the switch statistic:'); print(ACFd.round(3).to_string(index=False))
N = chrom_lengths(load_grid())
def band_power(vecs):
    out = {}
    for c, v in zip(CHR, vecs):
        g = order[order.chr == c]; n = N[c]; x = np.zeros(n); x[g.grid_index.values - 1] = v - v.mean()
        P = np.abs(np.fft.rfft(x)[1:n // 2 + 1]) ** 2; per = n / np.arange(1, len(P) + 1)
        for b, pw in zip(pd.cut(per, BANDS), P): out[b] = out.get(b, 0) + pw
    tot = sum(out.values()); return pd.Series({str(k): v / tot for k, v in out.items()})
bp_obs = band_power(byc); bp_null = pd.DataFrame([band_power([rng.permutation(v) for v in byc]) for _ in range(200)])
BP = pd.DataFrame({'band': bp_obs.index, 'fraction_power': bp_obs.values, 'null_mean': bp_null.mean().reindex(bp_obs.index).values, 'null_sd': bp_null.std().reindex(bp_obs.index).values})
BP['z'] = (BP.fraction_power - BP.null_mean) / BP.null_sd; BP.round(5).to_csv(tab('Table_S17e_F4_switch_power_by_band.csv'), index=False)
print('\nPower of the switch statistic by period band (z vs gene-order null):'); print(BP.round(3).to_string(index=False))
V = pd.read_csv(tab('Table_S5c_contiguous_DE_neighbourhoods.csv')); nb_genes = set(g for gs_ in V.genes for g in str(gs_).split(','))
P = pd.read_pickle(inter('pairs.pkl')); P['coupled'] = (np.sign(P.t1) == np.sign(P.t2)) & (P.t1.abs() > 3) & (P.t2.abs() > 3)
cpl = set(P[P.coupled].n1) | set(P[P.coupled].n2)
rows = []
for lab_, col in [('switch up', 'switch_gene'), ('switch down', 'switch_gene_down')]:
    g = SW[SW[col]].gene_name; allg = SW.gene_name
    for setname, S_ in [('DE neighbourhood', nb_genes), ('coupled pair', cpl)]:
        a = g.isin(S_).sum(); n = len(g); bg = allg.isin(S_).mean()
        rows.append({'switch_set': lab_, 'n_genes': n, 'in': setname, 'observed': a, 'fraction': a / max(n, 1), 'background': bg, 'odds_P': stats.binomtest(int(a), int(n), bg, alternative='greater').pvalue if n else np.nan})
EN_ = pd.DataFrame(rows); EN_.round(4).to_csv(tab('Table_S17f_F4_switch_in_neighbourhoods.csv'), index=False); print(); print(EN_.round(3).to_string(index=False))

# ---------- (5) cell of origin
TD = pd.read_csv(tab('Table_S8d_within_type_t_cirrhosis_vs_healthy.csv'), index_col=0)
sw = SW[SW.switch_gene]; cols = [c for c in ['Mesenquima_HSC', 'Endotelio', 'Fagocito_mononuclear', 'Colangiocito', 'T_NK', 'B', 'Plasma', 'pDC'] if c in TD.columns]
CO = TD.reindex(sw.index)[cols]; CO['gene_name'] = sw.gene_name; CO['t_step_F4_comp_adj'] = sw.t_step_F4_comp_adj
CO = CO.sort_values('t_step_F4_comp_adj', ascending=False); CO.round(3).to_csv(tab('Table_S17g_F4_switch_genes_cell_of_origin.csv'))
CO['cell_type_max'] = CO[cols].apply(lambda r: r.idxmax() if r.notna().any() else 'not detected', axis=1)
CO.round(3).to_csv(tab('Table_S17g_F4_switch_genes_cell_of_origin.csv'))
print('\ncell type with highest within-type induction among switch genes:'); print(CO.cell_type_max.value_counts().to_string())
print('\ntop switch genes:', ', '.join(CO.gene_name.head(40)))

# ---------- (6) timing of lineage-intrinsic targets: early/linear versus switch at F4 (Supplementary Table S17h)
TG = pd.read_csv(tab('Table_S9a_candidate_targets_lineage_intrinsic.csv'), index_col=0)
TG = TG.join(SW[['t_trend_comp_adj', 't_step_F4_comp_adj', 'q_adj_step', 'switch_gene']], how='left')
lin = ['t_sc_Mesenquima_HSC', 't_sc_Fagocito_mononuclear', 't_sc_Endotelio', 't_sc_Colangiocito']
TG['lineage'] = TG[lin].idxmax(axis=1).str.replace('t_sc_', '')
TG['timing'] = np.where(TG.switch_gene.fillna(False), 'switch at F4', np.where(TG.t_trend_comp_adj > 2, 'early, linear', 'other'))
TG.round(4).to_csv(tab('Table_S17h_target_timing.csv'))
print('\ntarget timing:', TG.timing.value_counts().to_dict()); print(pd.crosstab(TG.lineage, TG.timing).to_string())

# ---------- (7) hepatocyte programme summary, using the hepatocyte-specific gene set of script 08
# (≥2-fold above every other population in the scRNA-seq reference; Supplementary Table S9b column hep_specificity_log2)
G9 = pd.read_csv(tab('Table_S9b_integration_bulk_sc_per_gene.csv'), index_col=0)
EX = pd.read_csv(tab('Table_S8e_mean_expression_by_population.csv'), index_col=0); oth = [c for c in EX.columns if c not in ('gene_name', 'Hepatocito')]
spec = (EX['Hepatocito'] - EX[oth].max(axis=1)).reindex(G9.index); hs = G9[spec > 1]
pd.DataFrame([{'hepatocyte_specific_genes': len(hs), 'down_raw_t_lt_-3': int((hs.t_bulk_stage < -3).sum()), 'down_adjusted_t_lt_-3': int((hs.t_bulk_comp_adjusted < -3).sum()),
               'up_adjusted_t_gt_3': int((hs.t_bulk_comp_adjusted > 3).sum()), 'median_t_raw': hs.t_bulk_stage.median(), 'median_t_adjusted': hs.t_bulk_comp_adjusted.median()}]).round(3).to_csv(tab('Table_S17i_hepatocyte_programme_summary.csv'), index=False)
RF9 = pd.read_csv(tab('Table_S9d_hepatocyte_drug_targets_reference.csv'), index_col=0).join(SW[['t_trend_comp_adj', 't_step_F4_comp_adj']])
RF9.round(4).to_csv(tab('Table_S9d_hepatocyte_drug_targets_reference.csv'))
print('\nhepatocyte programme:', pd.read_csv(tab('Table_S17i_hepatocyte_programme_summary.csv')).to_dict('records'))
print(RF9[['gene_name', 't_bulk_stage', 't_bulk_comp_adjusted', 't_step_F4_comp_adj', 'hep_specificity_log2', 'development_stage']].sort_values('t_bulk_comp_adjusted').round(2).to_string(index=False))
