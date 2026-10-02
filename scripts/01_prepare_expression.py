"""01 — Load counts and GEO metadata, define conditions, normalise, remove cohort effect protecting stage.
Inputs : data/raw/counts/counts_<GSE>.tsv, data/raw/metadata/metadata_cruda_<GSE>.tsv, data/raw/grid/*rejilla_genes.tsv
Outputs: results/tables/Table_S1_samples.csv, results/intermediate/{counts,expr,expr_adj,meta}.pkl
"""
import os, re, glob
import numpy as np, pandas as pd
from common import *

grid = load_grid()
counts, rows = {}, []
for c in COHORTS:
    counts[c] = pd.read_csv(os.path.join(RAW, 'counts', f'counts_{c}.tsv'), sep='\t', index_col=0)
    d = pd.read_csv(os.path.join(RAW, 'metadata', f'metadata_{c}.tsv'), sep='\t')
    ch = [col for col in d.columns if col.startswith('characteristics_ch1')]
    for _, r in d.iterrows():
        kv = {}
        for col in ch:
            v = str(r[col])
            if ':' in v:
                k, val = v.split(':', 1); kv[k.strip().lower()] = val.strip()
        rows.append({'gsm': r.geo_accession, 'cohorte': c, 'source': r.get('source_name_ch1', ''), **kv})
M = pd.DataFrame(rows).set_index('gsm')


def condition(r):
    fs = str(r.get('fibrosis stage'))
    if r.cohorte == 'GSE162694' and fs == 'normal liver histology': return 'Normal'
    if r.cohorte == 'GSE135251' and r.get('group in paper') == 'control': return 'Control'
    if r.cohorte == 'GSE130970' and fs == '0' and float(r.get('nafld activity score', 'nan')) == 0: return 'Normal'
    return 'F' + fs

M['estadio'] = M.apply(condition, axis=1)
M['sexo'] = M.get('sex', pd.Series(index=M.index, dtype=str)).fillna('NA').astype(str).str[0].str.upper()
M['edad'] = pd.to_numeric(M.get('age', pd.Series(index=M.index)).fillna(M.get('age at biopsy', pd.Series(index=M.index))), errors='coerce')
M['nas'] = pd.to_numeric(M.get('nas score', pd.Series(index=M.index)).fillna(M.get('nafld activity score', pd.Series(index=M.index))), errors='coerce')
M['orden'] = M.estadio.map(ORDER_MAP)
all_samples = [s for c in COHORTS for s in counts[c].columns]
M = M.loc[all_samples]
M[['cohorte', 'estadio', 'orden', 'sexo', 'edad', 'nas']].to_csv(tab('Table_S1_samples.csv'))
print(pd.crosstab(M.cohorte, M.estadio))

# genes: on grid and median count >= 10 in every cohort
keep = set(grid.gene_id)
for c in COHORTS:
    keep &= set(counts[c].index[counts[c].median(axis=1) >= 10])
keep = grid[grid.gene_id.isin(keep)].copy()
print('genes retained:', len(keep))

# median-of-ratios normalisation per cohort, log2
L = {}
for c in COHORTS:
    d = counts[c].loc[keep.gene_id]
    lg = np.log(d.replace(0, np.nan)); ref = lg.mean(axis=1)
    sf = np.exp((lg.sub(ref, axis=0)).median())
    L[c] = np.log2(d / sf + 1)
X = pd.concat(L.values(), axis=1)[M.index]

# batch correction protecting condition: expr ~ cohort + condition, subtract cohort terms only
D = pd.get_dummies(M[['cohorte', 'estadio']], drop_first=True).astype(float); D.insert(0, 'const', 1.0)
beta = np.linalg.lstsq(D.values, X.values.T, rcond=None)[0]
coh = [i for i, c in enumerate(D.columns) if c.startswith('cohorte_')]
A = pd.DataFrame(X.values - (D.values[:, coh] @ beta[coh]).T, index=X.index, columns=X.columns)

pd.to_pickle(counts, inter('counts.pkl'))
pd.to_pickle((X, keep), inter('expr.pkl'))
pd.to_pickle(A, inter('expr_adj.pkl'))
M.to_pickle(inter('meta.pkl'))


def r2(Y, F):
    Dm = pd.get_dummies(F).values.astype(float); H = Dm @ np.linalg.pinv(Dm)
    Yc = Y - Y.mean(axis=1, keepdims=True)
    return 1 - np.sum((Yc - (H @ Yc.T).T) ** 2) / np.sum(Yc ** 2)

print('R2 cohort before/after: %.3f / %.3f ; R2 condition before/after: %.3f / %.3f' %
      (r2(X.values, M.cohorte), r2(A.values, M.cohorte), r2(X.values, M.estadio), r2(A.values, M.estadio)))
