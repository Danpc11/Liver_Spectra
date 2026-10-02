"""06 — GSE136103 (Ramachandran 2019): QC, marker annotation, pseudobulk per donor x type, within-type
cirrhosis-vs-healthy statistic, spatial autocorrelation, adjacent concordant pairs, within-type neighbour co-expression.
Streams one sample at a time (fits in ~3 GB RAM). Outputs: Table_S8_*
"""
import os, re, glob, pickle, zlib
import numpy as np, pandas as pd, scipy.io
from scipy import stats
from common import *

X, keep = pd.read_pickle(inter('expr.pkl')); K = pd.read_pickle(inter('K_tad.pkl')); grid_genes = set(keep.gene_id)
order = K.sort_values(['chr', 'grid_index']); pairs = []
for c, g in order.groupby('chr'):
    ids = g.index.values; pairs += list(zip(ids[:-1], ids[1:]))
pairs = pd.DataFrame(pairs, columns=['g1', 'g2'])
MARK = {'Mesenquima_HSC': ['COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM', 'PDGFRB', 'ACTA2', 'RGS5', 'MYH11', 'COL6A3', 'PDGFRA', 'TAGLN'],
        'Endotelio': ['PECAM1', 'VWF', 'CDH5', 'CLEC4G', 'STAB2', 'PLVAP', 'ACKR1', 'FCN3', 'OIT3', 'KDR'],
        'Colangiocito': ['KRT7', 'KRT19', 'EPCAM', 'SOX9', 'CFTR', 'TACSTD2', 'ANXA4'],
        'Hepatocito': ['ALB', 'APOA1', 'HP', 'TTR', 'CYP3A4', 'APOC3', 'TF', 'SERPINA1'],
        'Fagocito_mononuclear': ['CD68', 'CD163', 'LYZ', 'C1QA', 'C1QB', 'MARCO', 'TREM2', 'S100A8', 'FCGR3A', 'CD14', 'VSIG4'],
        'T_NK': ['CD3D', 'CD3E', 'CD2', 'NKG7', 'GNLY', 'IL7R', 'CD8A', 'KLRD1', 'TRAC'],
        'B': ['MS4A1', 'CD79A', 'CD79B', 'CD19'], 'Plasma': ['MZB1', 'JCHAIN', 'IGKC', 'XBP1', 'IGHG1'], 'pDC': ['LILRA4', 'IRF7', 'CLEC4C', 'GZMB']}
files = [f for f in sorted(glob.glob(os.path.join(EXT, 'GSE136103', '*_matrix.mtx.gz'))) if not any(s in f.lower() for s in ['blood', 'mouse'])]
PB, NC, CELLS, COEX = {}, {}, [], []; rng = np.random.default_rng(0)
for f in files:
    base = f.replace('_matrix.mtx.gz', ''); gsm, name = os.path.basename(base).split('_', 1)
    genes = pd.read_csv(base + '_genes.tsv.gz', sep='\t', header=None, names=['id', 'sym']); Mx = scipy.io.mmread(f).tocsc().astype(np.float32)
    nUMI = np.asarray(Mx.sum(0)).ravel(); ngene = np.asarray((Mx > 0).sum(0)).ravel(); mt = genes.sym.str.startswith('MT-').values; pmt = np.asarray(Mx[mt].sum(0)).ravel() / np.maximum(nUMI, 1)
    ok = (nUMI >= 500) & (ngene >= 300) & (pmt < 0.3); Mx = Mx[:, ok]; nUMI = nUMI[ok]; L = Mx.multiply(1e4 / nUMI[None, :]).tocsr(); L.data = np.log1p(L.data)
    idx = {s: i for i, s in enumerate(genes.sym.values)}
    S = np.vstack([np.asarray(L[[idx[g] for g in gl if g in idx]].mean(0)).ravel() for gl in MARK.values()]); Sz = (S - S.mean(1, keepdims=True)) / (S.std(1, keepdims=True) + 1e-9)
    best = Sz.argmax(0); srt = np.sort(Sz, 0); types = np.array(list(MARK)); ct = types[best]; ct[(srt[-1] < 1.0) | ((srt[-1] - srt[-2]) < 0.75)] = 'unassigned'
    cond = 'healthy' if 'healthy' in name.lower() else 'cirrhotic'; donor = re.match(r'(healthy\d|cirrhotic\d)', name.lower()).group(1)
    gi = np.array([i for i, g in enumerate(genes.id) if g in grid_genes]); gids = genes.id.values[gi]; Lg = L[gi]
    for t in np.unique(ct):
        if t == 'unassigned': continue
        cells = np.where(ct == t)[0]; pb = np.asarray(Mx[gi][:, cells].sum(1)).ravel()
        PB.setdefault((donor, cond, t), np.zeros(len(gi))); PB[(donor, cond, t)] += pb; NC[(donor, cond, t)] = NC.get((donor, cond, t), 0) + len(cells)
        if len(cells) >= 150:
            Am = Lg[:, cells].toarray(); expr = (Am > 0).mean(1) >= 0.10; Am = Am[expr]; sub = gids[expr]; pos = {g: j for j, g in enumerate(sub)}
            if len(sub) < 500: continue
            Am = (Am - Am.mean(1, keepdims=True)) / (Am.std(1, keepdims=True) + 1e-9); pp = pairs[pairs.g1.isin(pos) & pairs.g2.isin(pos)]
            r = np.einsum('ij,ij->i', Am[[pos[g] for g in pp.g1]], Am[[pos[g] for g in pp.g2]]) / Am.shape[1]
            i1 = rng.integers(0, len(sub), len(pp)); i2 = rng.integers(0, len(sub), len(pp)); rn = (np.einsum('ij,ij->i', Am[i1], Am[i2]) / Am.shape[1])[i1 != i2]
            COEX.append({'sample': name, 'donor': donor, 'cond': cond, 'type': t, 'n_cells': len(cells), 'n_genes': len(sub), 'r_neighbours': r.mean(), 'r_random': rn.mean()})
    CELLS.append(pd.DataFrame({'gsm': gsm, 'sample': name, 'donor': donor, 'cond': cond, 'type': ct}))
cells = pd.concat(CELLS); cells.to_csv(tab('Table_S8a_single_cell_annotation.csv'), index=False); pd.DataFrame(COEX).to_csv(tab('Table_S8b_within_type_neighbour_coexpression.csv'), index=False)
pickle.dump((PB, NC, gids), open(inter('pb_sc.pkl'), 'wb'))

# within-type DE and spatial autocorrelation
names = keep.set_index('gene_id').gene_name; order = K.reindex(gids).dropna(subset=['grid_index']).sort_values(['chr', 'grid_index'])
types = sorted(set(t for _, _, t in PB)); res, TT = [], {}
MIN_CELLS = {'Hepatocito': 40}
for t in types:
    # hepatocytes are poorly captured by tissue dissociation in this dataset (≈2% of cells), so they are admitted with ≥40 cells per donor;
    # all other populations need ≥100 cells per donor. Hepatocyte estimates are therefore noisier and are interpreted with that caveat.
    min_cells = MIN_CELLS.get(t, 100)
    keys = [k for k in PB if k[2] == t and NC[k] >= min_cells]; h = [k for k in keys if k[1] == 'healthy']; c_ = [k for k in keys if k[1] == 'cirrhotic']
    if len(h) < 3 or len(c_) < 3: continue
    Mx = np.column_stack([PB[k] for k in h + c_]); cpm = np.log2(Mx / Mx.sum(0) * 1e6 + 1); expr = (Mx >= 5).mean(1) >= 0.5
    tstat, p = stats.ttest_ind(cpm[:, len(h):], cpm[:, :len(h)], axis=1, equal_var=False); T = pd.Series(tstat, index=gids); T[~expr] = np.nan
    Tg = T.reindex(order.index).dropna(); byc = [Tg.loc[g.index.intersection(Tg.index)].values for _, g in order.groupby('chr')]; byc = [v for v in byc if len(v) > 20]
    a1 = acf_by_chr(byc, (1,))[0]; rng_t = np.random.default_rng(zlib.crc32(t.encode()))  # per-population seed: results do not depend on which populations are tested
    nul = np.array([acf_by_chr([rng_t.permutation(v) for v in byc], (1,))[0] for _ in range(200)]); TT[t] = T
    res.append({'type': t, 'healthy_donors': len(h), 'cirrhotic_donors': len(c_), 'cells': sum(NC[k] for k in keys), 'genes': int(Tg.notna().sum()), 'acf_lag1': a1, 'null_p97.5': np.quantile(nul, .975), 'z': (a1 - nul.mean()) / nul.std()})
pd.DataFrame(res).round(4).to_csv(tab('Table_S8c_within_type_spatial_autocorrelation.csv'), index=False); print(pd.DataFrame(res).round(3))
TD = pd.DataFrame(TT); TD['gene_name'] = names.reindex(TD.index).values; TD.to_csv(tab('Table_S8d_within_type_t_cirrhosis_vs_healthy.csv'))
# mean expression per population (log2 CPM over all donors) -> lineage specificity of each gene
EX = {}
for t in types:
    keys = [k for k in PB if k[2] == t and NC[k] >= MIN_CELLS.get(t, 100)]
    if len(keys) < 4: continue
    Mx = np.column_stack([PB[k] for k in keys]); EX[t] = np.log2(Mx / Mx.sum(0) * 1e6 + 1).mean(1)
EXd = pd.DataFrame(EX, index=gids); EXd['gene_name'] = names.reindex(EXd.index).values; EXd.round(4).to_csv(tab('Table_S8e_mean_expression_by_population.csv'))
print('populations tested:', [r['type'] for r in res])
