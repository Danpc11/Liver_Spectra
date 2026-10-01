"""04 — DESeq2 (pyDESeq2) ordinal-stage and F4-vs-Normal models; programme enrichment; spatial autocorrelation
of the DE statistic; adjacent DE pairs and contiguous neighbourhoods. Outputs: Table_S4_*, Table_S5_*
"""
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from common import *

counts = pd.read_pickle(inter('counts.pkl')); X, keep = pd.read_pickle(inter('expr.pkl')); M = pd.read_pickle(inter('meta.pkl'))
names = keep.set_index('gene_id').gene_name.astype(str)
cnt = pd.concat([counts[c].loc[keep.gene_id] for c in COHORTS], axis=1)[M.index].T

# F4 vs Normal (cohorts containing both)
s = M.index[(M.estadio.isin(['Normal', 'F4'])) & (M.cohorte.isin(FULL))]
md = M.loc[s, ['cohorte', 'estadio']].rename(columns={'estadio': 'cond'})
dds = DeseqDataSet(counts=cnt.loc[s].astype(int), metadata=md, design='~ cohorte + cond', quiet=True); dds.deseq2()
st = DeseqStats(dds, contrast=['cond', 'F4', 'Normal'], quiet=True); st.summary(); de = st.results_df; de['gene_name'] = names.reindex(de.index).values
de.to_csv(tab('Table_S4a_DESeq2_F4_vs_Normal.csv'))
# ordinal
s2 = M.index[M.estadio != 'Control']; md2 = M.loc[s2, ['cohorte', 'orden']].copy(); md2['orden'] = md2.orden.astype(float)
dds2 = DeseqDataSet(counts=cnt.loc[s2].astype(int), metadata=md2, design='~ cohorte + orden', quiet=True); dds2.deseq2()
st2 = DeseqStats(dds2, contrast=np.array([0, 0, 0, 1.0]), quiet=True); st2.summary(); de2 = st2.results_df; de2['gene_name'] = names.reindex(de2.index).values
de2.to_csv(tab('Table_S4b_DESeq2_ordinal_stage.csv')); pd.to_pickle((de, de2), inter('de.pkl'))
print('DE padj<0.05: F4 vs Normal', (de.padj < .05).sum(), '| ordinal', (de2.padj < .05).sum())

# programme enrichment (ordinal, padj<0.01)
D = keep.set_index('gene_id').join(de2[['stat', 'padj']]); D['stat'] = D.stat.fillna(0)
D['de_up'] = (D.padj < .01) & (D.stat > 0); D['de_dn'] = (D.padj < .01) & (D.stat < 0)
nm = D.gene_name.astype(str); rows = []
for name, pat in PROGRAMMES.items():
    m = nm.str.match(pat).values
    for lab, sel in [('up', D.de_up.values), ('down', D.de_dn.values)]:
        a = (m & sel).sum(); odds, p = fisher_exact([[a, (m & ~sel).sum()], [(~m & sel).sum(), (~m & ~sel).sum()]], alternative='greater')
        rows.append({'programme': name, 'direction': lab, 'n': int(a), 'expected': round(m.sum() * sel.mean(), 1), 'odds': round(odds, 2), 'p': p})
pd.DataFrame(rows).to_csv(tab('Table_S4c_programme_enrichment_DE.csv'), index=False)

# spatial autocorrelation of the DE statistic (compact gene order)
order = D.sort_values(['chr', 'grid_index']); byc = [order[order.chr == c].stat.values for c in CHR]
rng = np.random.default_rng(1); lags = (1, 2, 3, 5, 10, 20, 50)
obs = acf_by_chr(byc, lags); null = np.array([acf_by_chr([rng.permutation(v) for v in byc], lags) for _ in range(200)])
pd.DataFrame({'lag': lags, 'acf_observed': obs, 'null_p2.5': np.quantile(null, .025, 0), 'null_p97.5': np.quantile(null, .975, 0), 'z': (obs - null.mean(0)) / null.std(0)}).round(4).to_csv(tab('Table_S5a_spatial_autocorrelation_DE.csv'), index=False)

# adjacent concordant DE pairs and runs >=3
def adj_pairs(col, Dx):
    return sum(np.sum((Dx[Dx.chr == c].sort_values('grid_index')[col].values.astype(int))[:-1] * (Dx[Dx.chr == c].sort_values('grid_index')[col].values.astype(int))[1:]) for c in CHR)
rows = []
for col in ['de_up', 'de_dn']:
    o = adj_pairs(col, D); nl = []
    for _ in range(200):
        Dp = D.copy()
        for c in CHR:
            idx = Dp.index[Dp.chr == c]; Dp.loc[idx, col] = rng.permutation(Dp.loc[idx, col].values)
        nl.append(adj_pairs(col, Dp))
    rows.append({'direction': col, 'observed': o, 'null_mean': np.mean(nl), 'null_sd': np.std(nl), 'z': (o - np.mean(nl)) / np.std(nl)})
pd.DataFrame(rows).to_csv(tab('Table_S5b_adjacent_DE_pairs.csv'), index=False); print(pd.DataFrame(rows))
D['sig'] = np.sign(D.stat) * (D.padj < .01).astype(int); runs = []
for c in CHR:
    g = D[D.chr == c].sort_values('grid_index').reset_index(); i = 0
    while i < len(g):
        sgn = g.sig[i]
        if sgn == 0: i += 1; continue
        j = i
        while j + 1 < len(g) and g.sig[j + 1] == sgn: j += 1
        if j - i + 1 >= 3:
            runs.append({'chr': c, 'start_grid': g.grid_index[i], 'end_grid': g.grid_index[j], 'n_genes': j - i + 1, 'direction': 'up' if sgn > 0 else 'down', 'mean_stat': g.stat[i:j + 1].mean(), 'genes': ','.join(g.gene_name[i:j + 1].astype(str))})
        i = j + 1
R = pd.DataFrame(runs).sort_values('n_genes', ascending=False); R.to_csv(tab('Table_S5c_contiguous_DE_neighbourhoods.csv'), index=False)
print('runs >=3:', len(R))
