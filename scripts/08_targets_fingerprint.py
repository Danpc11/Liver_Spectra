"""08 — Lineage-intrinsic target prioritisation (bulk x composition x single cell) and the spectral-fingerprint
control (classification with real vs shuffled gene order). Outputs: Table_S9_*, Table_S11_*
"""
import os, re, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score
from scipy import stats
from common import *

X, keep = pd.read_pickle(inter('expr.pkl')); A = pd.read_pickle(inter('expr_adj.pkl')); M = pd.read_pickle(inter('meta.pkl')); de, de2 = pd.read_pickle(inter('de.pkl'))
T = pd.read_csv(tab('Table_S6c_stage_effect_per_gene_with_without_composition.csv')).set_index('gene_id'); TD = pd.read_csv(tab('Table_S8d_within_type_t_cirrhosis_vs_healthy.csv'), index_col=0)
P = pd.read_pickle(inter('pairs.pkl')); V = pd.read_csv(tab('Table_S5c_contiguous_DE_neighbourhoods.csv'))
DRUG = {'secreted_ligand': r'^(CXCL|CCL|IL[0-9]|TGFB[123]$|PDGF[ABCD]$|THBS|SPP1|TNFSF|CCN[1-6]|LIF$|OSM$|BMP|WNT|FGF[0-9]|VEGF|ANGPT|IGF[12]$|IGFBP|LGALS|S100A|CHI3L1|MMP|TIMP|LOX|SERPINE1|PLAU$|PLAUR$|CTHRC1|POSTN|TNC$|FN1$|APOE$)',
        'membrane_receptor': r'(R$|R[0-9]$|RA$|RB$|RB[0-9]$|RL[0-9]$|^PDGFR|^TGFBR|^ITG|^CD[0-9]+|^TNFRSF|^FZD|^NOTCH|^EGFR|^MET$|^KIT$|^FLT|^KDR|^TEK$|^ACKR|^GPR|^ADGR|^S1PR|^LPAR|^CCR|^CXCR|^TREM|^MARCO)',
        'kinase': r'(^MAP[0-9]?K|^PRK|^CDK|^JAK|^SRC$|^SYK$|^BTK$|^ROCK|^LIMK|^PIK3|^AKT|^MTOR$|^PTK2$|^TYK2$)',
        'metabolic_enzyme': r'(^HSD|^CYP|^ACAC|^SCD$|^DGAT|^FASN$|^PNPLA|^HMGCR$|^NAMPT$|^IDO1$|^ALOX|^PTGS|^MAOA|^MAOB|^DPP4$|^ACE$|^NOX|^CYBB$|^P4HA|^PLOD|^TGM2$|^GLUL$)'}
G = pd.DataFrame({'gene_name': TD.gene_name, 't_bulk_stage': T.t_stage.reindex(TD.index), 't_bulk_comp_adjusted': T.t_stage_comp_adjusted.reindex(TD.index), 'padj_bulk': de2.padj.reindex(TD.index)})
sc_cols = [c for c in ['Mesenquima_HSC', 'Colangiocito', 'Fagocito_mononuclear', 'Endotelio', 'T_NK', 'B', 'Plasma', 'pDC'] if c in TD.columns]
for c in sc_cols: G['t_sc_' + c] = TD[c]
P2 = P[(~P.paralog_family) & (np.sign(P.t1) == np.sign(P.t2)) & (P.t1.abs() > 3) & (P.t2.abs() > 3)]; partner = {}
for _, r in P2.iterrows(): partner.setdefault(r.g1, []).append(r.n2); partner.setdefault(r.g2, []).append(r.n1)
G['coupled_neighbour'] = [','.join(partner.get(g, [])) for g in G.index]; block_genes = set(g for gs in V.genes for g in str(gs).split(',')); G['in_DE_neighbourhood'] = G.gene_name.isin(block_genes)
G['druggable_class'] = G.gene_name.map(lambda n: ';'.join(k for k, p in DRUG.items() if re.search(p, str(n))))
main = ['t_sc_' + c for c in ['Mesenquima_HSC', 'Colangiocito', 'Fagocito_mononuclear', 'Endotelio'] if 't_sc_' + c in G.columns]
G['cell_type_up_sc'] = [','.join(c.replace('t_sc_', '') for c in main if pd.notna(r[c]) and r[c] > 2) for _, r in G.iterrows()]; G['t_sc_max'] = G[main].max(axis=1)
C = G[(G.t_bulk_stage > 4) & (G.t_sc_max > 2) & (G.t_bulk_comp_adjusted > 1.5)].copy()
C['score'] = C.t_bulk_comp_adjusted.clip(0, 10) / 10 + C.t_sc_max.clip(0, 8) / 8 + (C.druggable_class != '').astype(float) + (C.coupled_neighbour != '').astype(float) * .5 + C.in_DE_neighbourhood.astype(float) * .5
DL = pd.read_csv(os.path.join(ROOT, 'data', 'curated', 'drug_landscape.csv')); C.index.name = 'gene_id'; C = C.reset_index().merge(DL, on='gene_name', how='left').set_index('gene_id'); C['development_stage'] = C.development_stage.fillna('not curated')
C.sort_values('score', ascending=False).round(3).to_csv(tab('Table_S9a_candidate_targets_lineage_intrinsic.csv')); G.round(3).to_csv(tab('Table_S9b_integration_bulk_sc_per_gene.csv'))
print('candidates:', len(C))

# ---- spectral fingerprint control (leave-one-cohort-out)
Wz, Wl = pd.read_pickle(inter('W_adj.pkl')); grid = load_grid(); N = chrom_lengths(grid); rng = np.random.default_rng(0)
F = {'spectrum': np.log(Wl).T, 'expression_2000': A.loc[A.var(axis=1).sort_values(ascending=False).index[:2000]].T, 'spectrum_shuffled_order': np.log(spectra(A, keep, N, rng=rng)).T}
per = pd.Series([N[c] / k for c, k in Wl.index], index=Wl.index); band = pd.cut(per, BANDS)
B = np.log(Wl).groupby([Wl.index.get_level_values(0), band.values], observed=True).mean().T; B.columns = [f'{c}_{b}' for c, b in B.columns]; F['bands_x_chromosome'] = B
M['group'] = np.where(M.estadio.isin(['F0', 'F1']), 'early', np.where(M.estadio.isin(['F3', 'F4']), 'advanced', np.where(M.estadio == 'Normal', 'normal', 'other')))
def loco(Xf, task):
    res = []
    for test_c in COHORTS:
        if task == 'early_vs_advanced': idx = M.index[M.group.isin(['early', 'advanced'])]; y = (M.loc[idx, 'group'] == 'advanced').astype(int)
        elif task == 'normal_vs_advanced': idx = M.index[M.group.isin(['normal', 'advanced'])]; y = (M.loc[idx, 'group'] == 'advanced').astype(int)
        else: idx = M.index[M.estadio != 'Control']; y = M.loc[idx, 'orden'].astype(float)
        tr = idx[M.loc[idx, 'cohorte'] != test_c]; te = idx[M.loc[idx, 'cohorte'] == test_c]
        if len(te) < 10 or (task != 'ordinal' and y.loc[te].nunique() < 2): continue
        nc = min(40, len(tr) - 1, Xf.shape[1])
        if task == 'ordinal': mdl = make_pipeline(StandardScaler(), PCA(nc, random_state=0), Ridge(alpha=10.0)).fit(Xf.loc[tr].values, y.loc[tr]); res.append(stats.spearmanr(mdl.predict(Xf.loc[te].values), y.loc[te])[0])
        else: mdl = make_pipeline(StandardScaler(), PCA(nc, random_state=0), LogisticRegression(C=0.5, max_iter=2000)).fit(Xf.loc[tr].values, y.loc[tr]); res.append(roc_auc_score(y.loc[te], mdl.predict_proba(Xf.loc[te].values)[:, 1]))
    return np.mean(res)
rows = [{'features': nm, 'n_features': Xf.shape[1], 'task': task, 'metric': 'Spearman' if task == 'ordinal' else 'AUC', 'value': loco(Xf, task)} for nm, Xf in F.items() for task in ['early_vs_advanced', 'normal_vs_advanced', 'ordinal']]
pd.DataFrame(rows).round(3).to_csv(tab('Table_S11_spectral_fingerprint_real_vs_shuffled.csv'), index=False); print(pd.DataFrame(rows).round(3))
