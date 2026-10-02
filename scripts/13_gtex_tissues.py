"""13 — Is the positional architecture hepatic or genomic? Per-sample spectra in 11 GTEx tissues on the same
positional grid; universal peaks per tissue; overlap with the liver-biopsy universal peaks; similarity of consensus
spectra between tissues. Two gene sets: (a) genes expressed in each tissue (median >= 10 counts), (b) genes expressed
in all tissues (common set). Up to 120 randomly chosen samples per tissue. Outputs: Table_S14a–d, ExtendedData_Fig2
"""
import os, glob, gzip, re, numpy as np, pandas as pd
from common import *

rng = np.random.default_rng(0); grid = load_grid(); N = chrom_lengths(grid); gidx = grid.set_index('gene_id')
files = sorted(glob.glob(os.path.join(EXT, 'gtex', '*gct.gz')))
def tname(f): return re.sub(r'^gene_reads_(adult_gtex_v11_|v10_)', '', os.path.basename(f)).replace('_gct.gz', '').replace('.gct.gz', '')
cnt = {}
for f in files:
    with gzip.open(f, 'rt') as fh: fh.readline(); fh.readline(); hdr = fh.readline().rstrip('\n').split('\t')
    samp = hdr[2:]; pick = list(rng.choice(samp, min(120, len(samp)), replace=False))
    d = pd.read_csv(f, sep='\t', skiprows=2, usecols=['Name'] + pick, index_col=0)
    d.index = d.index.str.split('.').str[0]; d = d[~d.index.duplicated()]; d = d[d.index.isin(gidx.index)]
    cnt[tname(f)] = d.astype(np.float32); print(tname(f), d.shape)
tissues = list(cnt)
expr = {t: set(d.index[d.median(axis=1) >= 10]) for t, d in cnt.items()}
common = set.intersection(*expr.values()); print('genes expressed in all tissues:', len(common))

def tissue_spectra(t, genes):
    d = cnt[t].loc[sorted(genes)]; lc = np.log2(d / cnt[t].sum() * 1e6 + 1)
    keep = grid[grid.gene_id.isin(lc.index)]
    return spectra(lc, keep, N)
res = {}
for gs_name in ['own', 'common']:
    W = {t: spectra_t for t, spectra_t in ((t, tissue_spectra(t, expr[t] if gs_name == 'own' else common)) for t in tissues)}
    frac = pd.DataFrame({t: (w > 3).mean(axis=1) for t, w in W.items()}); mlog = pd.DataFrame({t: np.log(w).mean(axis=1) for t, w in W.items()})
    res[gs_name] = (frac, mlog)
    frac.round(4).to_csv(tab(f'Table_S14a_gtex_fraction_samples_w_gt3_{gs_name}_genes.csv')); mlog.round(4).to_csv(tab(f'Table_S14b_gtex_mean_log_power_{gs_name}_genes.csv'))

# liver biopsy universal peaks
U = pd.read_csv(tab('Table_S2a_consensus_spectrum_universal_peaks.csv')); U['chr'] = U.chr.astype(str); U = U.set_index(['chr', 'k'])
biopsy_univ = U.index[U.universal]
rows = []
for gs_name, (frac, mlog) in res.items():
    univ = {t: set(frac.index[frac[t] >= 0.9]) for t in tissues}
    for t in tissues:
        ix = frac.index.intersection(biopsy_univ)
        rows.append({'gene_set': gs_name, 'tissue': t, 'n_samples': cnt[t].shape[1], 'n_genes': len(expr[t]) if gs_name == 'own' else len(common), 'universal_peaks_tissue': len(univ[t]),
                     'biopsy_universal_peaks_testable': len(ix), 'frac_biopsy_peaks_universal_in_tissue': np.mean([p in univ[t] for p in ix]),
                     'jaccard_with_GTEx_liver': len(univ[t] & univ['liver']) / max(1, len(univ[t] | univ['liver'])),
                     'r_mean_log_power_with_biopsy_liver': np.corrcoef(mlog[t].reindex(U.index).dropna(), U.mean_log_w.reindex(mlog[t].reindex(U.index).dropna().index))[0, 1] if gs_name == 'own' else np.nan})
S = pd.DataFrame(rows); S.round(4).to_csv(tab('Table_S14c_gtex_universal_peaks_by_tissue.csv'), index=False); print(S.round(3).to_string(index=False))
for gs_name, (frac, mlog) in res.items():
    Cm = mlog.corr(); Cm.round(4).to_csv(tab(f'Table_S14d_gtex_spectrum_similarity_{gs_name}_genes.csv')); print(f'\n[{gs_name}] correlation of consensus spectra between tissues:\n', Cm.round(2).to_string())
pd.to_pickle(res, inter('gtex_res.pkl'))

# ---- enrichment of biopsy peaks in GTEx liver vs other tissues; shared core
from scipy.stats import fisher_exact
frac, mlog = res['own']; ix = frac.index.intersection(biopsy_univ); rest = frac.index.difference(biopsy_univ)
rows = []
for t in tissues:
    u_t = frac[t] >= 0.9; a = int(u_t.loc[ix].sum()); b = len(ix) - a; c = int(u_t.loc[rest].sum()); d = len(rest) - c
    o, p = fisher_exact([[a, b], [c, d]]); rows.append({'tissue': t, 'biopsy_peaks_universal': a, 'of': len(ix), 'odds_vs_background': o, 'p': p})
F = pd.DataFrame(rows); F.round(5).to_csv(tab('Table_S14e_biopsy_peak_enrichment_by_tissue.csv'), index=False); print(F.round(4).to_string(index=False))
nuniv = (frac >= 0.9).sum(axis=1); core = nuniv[nuniv >= 8]; print('frequencies universal in >=8 of 11 tissues (genomic core):', len(core), '| in exactly 1 tissue:', int((nuniv == 1).sum()), '| in liver only:', int(((frac['liver'] >= 0.9) & (nuniv == 1)).sum()))
pd.DataFrame({'n_tissues_universal': nuniv}).reset_index().to_csv(tab('Table_S14f_peak_sharing_across_tissues.csv'), index=False)

# ---- Extended Data Figure 2
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, scienceplots
plt.style.use(['science', 'nature', 'no-latex'])
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial'], 'font.size': 7, 'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 6, 'xtick.top': False, 'ytick.right': False,
                     'xtick.minor.visible': False, 'ytick.minor.visible': False, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False, 'pdf.fonttype': 42})
def lab(ax, t, dx=-0.22, dy=1.14): ax.text(dx, dy, t.upper(), transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')
pretty = {t: t.replace('_', ' ').replace('brain cerebellar hemisphere', 'brain (cerebellum)').replace('adipose visceral omentum', 'adipose (visceral)').replace('adipose subcutaneous', 'adipose (subcut.)').replace('colon transverse', 'colon').replace('kidney cortex', 'kidney').replace('muscle skeletal', 'muscle').replace('artery aorta', 'aorta').replace('whole blood', 'blood') for t in tissues}
fig = plt.figure(figsize=(183 / 25.4, 183 / 25.4 * 0.62)); gs_ = fig.add_gridspec(2, 3, hspace=0.75, wspace=0.6, width_ratios=[1, 1, 1.25])
ax = fig.add_subplot(gs_[0, 0]); s1 = S[S.gene_set == 'own'].sort_values('universal_peaks_tissue')
ax.barh(range(len(s1)), s1.universal_peaks_tissue, color=['#D55E00' if t == 'liver' else '0.55' for t in s1.tissue], lw=0); ax.set_yticks(range(len(s1))); ax.set_yticklabels([pretty[t] for t in s1.tissue]); ax.set_xlabel('Universal peaks\n(>3× in ≥90% of samples)'); lab(ax, 'a', -0.55)
ax = fig.add_subplot(gs_[0, 1]); f1 = F.sort_values('odds_vs_background')
ax.barh(range(len(f1)), np.log2(f1.odds_vs_background), color=['#D55E00' if t == 'liver' else '0.55' for t in f1.tissue], lw=0); ax.set_yticks(range(len(f1))); ax.set_yticklabels([pretty[t] for t in f1.tissue]); ax.axvline(0, color='k', lw=0.5)
ax.set_xlabel('log2 odds: biopsy-liver peaks\nuniversal in GTEx tissue'); lab(ax, 'b', -0.55)
ax = fig.add_subplot(gs_[:, 2]); Cm = res['common'][1].corr(); order = ['liver', 'kidney_cortex', 'pancreas', 'colon_transverse', 'lung', 'adipose_subcutaneous', 'adipose_visceral_omentum', 'artery_aorta', 'muscle_skeletal', 'brain_cerebellar_hemisphere', 'whole_blood']
Cm = Cm.loc[order, order]; im = ax.imshow(Cm.values, cmap='viridis', vmin=0, vmax=1); ax.set_xticks(range(len(order))); ax.set_xticklabels([pretty[t] for t in order], rotation=60, ha='right'); ax.set_yticks(range(len(order))); ax.set_yticklabels([pretty[t] for t in order])
cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03); cb.set_label('r, consensus spectra (common genes)', fontsize=6); ax.tick_params(length=0); ax.spines[['left', 'bottom']].set_visible(False); lab(ax, 'd', -0.45, 1.06)
ax = fig.add_subplot(gs_[1, :2]); h = nuniv[nuniv > 0].value_counts().sort_index()
ax.bar(h.index, h.values, color='0.55', lw=0, width=0.7); ax.set_xticks(range(1, 12)); ax.set_xlabel('Number of tissues in which a frequency is a universal peak'); ax.set_ylabel('Frequencies'); ax.set_yscale('log'); lab(ax, 'c', -0.1)
fig.savefig(os.path.join(FIG, 'genome_research', 'ExtendedData_Fig2_gtex_tissues.pdf'), bbox_inches='tight'); fig.savefig(os.path.join(FIG, 'genome_research', 'ExtendedData_Fig2_gtex_tissues.png'), dpi=400, bbox_inches='tight'); print('figure ok')
