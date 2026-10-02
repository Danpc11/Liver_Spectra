"""18 — Validation in single-nucleus RNA-seq across the MASLD spectrum (GSE202379; Gribben et al., Nature 2024).
47 donors, 59 samples (healthy and end-stage explants contribute three lobes each), SAF fibrosis stage per biopsy.
Step A (per sample, resumable): marker-based annotation of nuclei, pseudobulk counts per sample × population, programme
scores per population. Step B: per-donor analyses — composition by stage, within-population stage effects, concordance
with the bulk composition-adjusted effect, hepatocyte class and drug receptors, neighbour coupling, the F4 switch.
Stage coding: healthy = Normal; NAFLD/NASH biopsies by SAF F0–F4; end-stage explants analysed as F4 (and excluded in a
sensitivity analysis). Outputs: Table_S18a–h
"""
import os, glob, gzip, sys, warnings, numpy as np, pandas as pd
warnings.filterwarnings('ignore')
from scipy import stats
from common import *

DD = os.path.join(EXT, 'GSE202379'); CACHE = os.path.join(INTER, 'gse202379'); os.makedirs(CACHE, exist_ok=True)
MARK = {'Hepatocyte': ['ALB', 'APOA1', 'APOC3', 'HP', 'TTR', 'CYP3A4', 'CYP2E1', 'SERPINA1', 'TF', 'APOB', 'PCK1', 'ASGR1'],
        'Cholangiocyte': ['KRT7', 'KRT19', 'EPCAM', 'SOX9', 'CFTR', 'ANXA4', 'BICC1', 'PKHD1', 'DCDC2'],
        'HSC': ['COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM', 'PDGFRB', 'RELN', 'LRAT', 'ACTA2', 'RGS5', 'CARMN', 'COLEC11'],
        'Endothelium': ['PECAM1', 'VWF', 'CDH5', 'STAB2', 'CLEC4G', 'FLT1', 'PTPRB', 'KDR', 'LDB2', 'EMCN', 'MMRN1'],
        'Macrophage': ['CD163', 'MARCO', 'VSIG4', 'C1QA', 'C1QB', 'CD68', 'MS4A7', 'F13A1', 'MRC1', 'TREM2', 'CSF1R'],
        'T/NK': ['PTPRC', 'CD247', 'IL7R', 'CD3E', 'NKG7', 'GNLY', 'SKAP1', 'THEMIS', 'CD96'],
        'B/Plasma': ['MS4A1', 'CD79A', 'BANK1', 'IGKC', 'JCHAIN', 'MZB1', 'IGHM']}
PROG = {'YAP_TAZ_targets': ['CCN1', 'CCN2', 'ANKRD1', 'AMOTL2', 'AXL', 'LATS2', 'TEAD1', 'TEAD4', 'CRIM1', 'F3', 'GADD45B', 'SERPINE1', 'TGFB2', 'THBS1', 'DKK1', 'WWC1', 'AJUBA'],
        'mechanosensing': ['PIEZO1', 'PIEZO2', 'ITGB1', 'ITGA5', 'ITGAV', 'ITGB5', 'PTK2', 'RHOA', 'ROCK1', 'ROCK2', 'YAP1', 'WWTR1', 'VCL', 'TLN1', 'FLNA', 'ACTN1', 'MYL9', 'TAGLN', 'CNN1', 'ACTA2'],
        'matrix_crosslinking': ['LOX', 'LOXL1', 'LOXL2', 'LOXL3', 'LOXL4', 'TGM2', 'PLOD1', 'PLOD2', 'PLOD3', 'P4HA1', 'P4HA2', 'PXDN', 'ELN', 'FBLN5', 'COL1A1', 'COL1A2', 'COL3A1']}

def meta():
    rows = {}
    for line in gzip.open(os.path.join(DD, 'GSE202379_series_matrix_txt.gz'), 'rt'):
        if line.startswith('!Sample_'):
            p = line.rstrip('\n').split('\t'); rows.setdefault(p[0], []).append([x.strip('"') for x in p[1:]])
    D = pd.DataFrame({'gsm': rows['!Sample_geo_accession'][0], 'title': rows['!Sample_title'][0]})
    for v in rows['!Sample_characteristics_ch1']: D[v[0].split(':')[0]] = [x.split(': ', 1)[1] if ': ' in x else x for x in v]
    D['F'] = D['saf score'].str.extract(r'F(\d)').astype(float)
    D['stage'] = np.where(D['disease status'] == 'Healthy control', 'Normal', np.where(D['disease status'] == 'end stage', 'F4', 'F' + D.F.fillna(-1).astype(int).astype(str)))
    D['end_stage'] = D['disease status'] == 'end stage'; D['orden'] = D.stage.map({'Normal': 0, 'F0': 1, 'F1': 2, 'F2': 3, 'F3': 4, 'F4': 5})
    return D.set_index('gsm')

def _read(f):
    X = pd.read_csv(f, index_col=0, dtype={0: str}).astype(np.float32)          # genes × nuclei
    genes = X.index.values; lib = X.sum(axis=0).values; L = np.log1p(X.values / lib * 1e4)
    return X, genes, L, {g: i for i, g in enumerate(genes)}

def step_a():
    """Pass 1: raw marker-set scores per nucleus (mean log-normalised expression), saved per sample."""
    for f in sorted(glob.glob(os.path.join(DD, 'GSM*raw_counts_csv.gz'))):
        gsm = os.path.basename(f).split('_')[0]; out = os.path.join(CACHE, gsm + '_scores.npz')
        if os.path.exists(out): continue
        X, genes, L, gi = _read(f)
        S = np.vstack([L[[gi[g] for g in v if g in gi]].mean(0) for v in MARK.values()])
        np.savez_compressed(out, S=S); print(gsm, S.shape[1], flush=True)

def step_b():
    """Pass 2: scale each marker score by its 99th percentile over all nuclei (z-scoring penalises the majority population,
    here hepatocytes), annotate by the top-scoring set, pseudobulk per sample × population.
    Validation: this reproduces the published composition (hepatocytes 68% vs 70%, cholangiocytes 5.3% vs 5.4%)."""
    allS = np.hstack([np.load(f)['S'] for f in sorted(glob.glob(os.path.join(CACHE, '*_scores.npz')))])
    q99 = np.percentile(allS, 99, axis=1, keepdims=True); types = np.array(list(MARK))
    for f in sorted(glob.glob(os.path.join(DD, 'GSM*raw_counts_csv.gz'))):
        gsm = os.path.basename(f).split('_')[0]; out = os.path.join(CACHE, gsm + '.npz')
        if os.path.exists(out): continue
        X, genes, L, gi = _read(f); Z = np.load(os.path.join(CACHE, gsm + '_scores.npz'))['S'] / q99
        srt = np.sort(Z, 0); ct = types[Z.argmax(0)]; ct[(srt[-1] < 0.3) | ((srt[-1] - srt[-2]) < 0.1)] = 'Unassigned'
        P = {k: L[[gi[g] for g in v if g in gi]].mean(0) for k, v in PROG.items()}
        pb, nc, pr = {}, {}, {}
        for t in list(MARK) + ['Unassigned']:
            m = ct == t; nc[t] = int(m.sum())
            if m.sum(): pb[t] = X.values[:, m].sum(1); pr[t] = [float(P[k][m].mean()) for k in PROG]
        np.savez_compressed(out, genes=genes, types=np.array(list(pb)), pb=np.vstack([pb[t] for t in pb]), nc=np.array([nc[t] for t in pb]), prog=np.array([pr[t] for t in pb]))
        print(gsm, X.shape[1], {k: v for k, v in nc.items() if v}, flush=True)

def load_pseudobulk():
    D = meta(); recs = {}
    for f in sorted(glob.glob(os.path.join(CACHE, 'GSM*.npz'))):
        if f.endswith('_scores.npz'): continue
        gsm = os.path.basename(f)[:-4]; z = np.load(f, allow_pickle=True)
        for t, v, n, pr in zip(z['types'], z['pb'], z['nc'], z['prog']): recs[(gsm, str(t))] = (v, int(n), pr)
        genes = z['genes']
    return D, recs, genes

def donor_matrix(D, recs, genes, ctype, min_nuclei=30, include_end_stage=False):
    """Sum lobes per donor; return donor table, counts (genes × donors) and nuclei-weighted programme scores."""
    rows = []; cols = {}; progs = {}
    for donor, g in D.groupby('patient id'):
        if g.end_stage.iloc[0] and not include_end_stage: continue
        v = sum(recs[(s, ctype)][0] for s in g.index if (s, ctype) in recs) if any((s, ctype) in recs for s in g.index) else None
        n = sum(recs[(s, ctype)][1] for s in g.index if (s, ctype) in recs)
        if v is None or n < min_nuclei: continue
        pr = sum(recs[(s, ctype)][2] * recs[(s, ctype)][1] for s in g.index if (s, ctype) in recs) / n
        cols[donor] = v; progs[donor] = pr
        rows.append({'donor': donor, 'stage': g.stage.iloc[0], 'orden': g.orden.iloc[0], 'sex': g.gender.iloc[0], 'end_stage': bool(g.end_stage.iloc[0]), 'nuclei': n})
    T = pd.DataFrame(rows).set_index('donor'); C = pd.DataFrame(cols, index=genes)[T.index]
    P = pd.DataFrame(progs, index=list(PROG)).T.loc[T.index]
    return T, C, P

def stage_t(T, C, min_counts=10):
    cpm = np.log2(C / C.sum(0) * 1e6 + 1); keep = (C >= min_counts).mean(1) >= 0.5; cpm = cpm[keep]
    X = np.column_stack([np.ones(len(T)), T.orden.values.astype(float), (T.sex == 'M').values.astype(float)])
    XtXi = np.linalg.pinv(X.T @ X); B = XtXi @ X.T @ cpm.values.T; R = cpm.values.T - X @ B; df = len(T) - 3
    se = np.sqrt((R ** 2).sum(0) / df * XtXi[1, 1]); return pd.Series(B[1] / se, index=cpm.index), df

def step_c():
    D, recs, genes = load_pseudobulk()
    types = list(MARK)
    # ---- C1 composition by stage (per donor, all lobes pooled)
    comp = []
    for donor, g in D.groupby('patient id'):
        n = {t: sum(recs[(s, t)][1] for s in g.index if (s, t) in recs) for t in types + ['Unassigned']}; tot = sum(n.values())
        comp.append({'donor': donor, 'stage': g.stage.iloc[0], 'orden': g.orden.iloc[0], 'end_stage': bool(g.end_stage.iloc[0]), 'nuclei': tot, **{t: n[t] / tot for t in types}})
    CP = pd.DataFrame(comp); CP.round(4).to_csv(tab('Table_S18a_snRNA_composition_per_donor.csv'), index=False)
    bio = CP[~CP.end_stage]
    print('composition vs stage (biopsies, Spearman):', {t: round(stats.spearmanr(bio.orden, bio[t])[0], 2) for t in types})
    print(CP.groupby('stage')[types].mean().reindex(['Normal', 'F0', 'F1', 'F2', 'F3', 'F4']).round(3).to_string())
    # ---- bulk references
    T6 = pd.read_csv(tab('Table_S6c_stage_effect_per_gene_with_without_composition.csv')).drop_duplicates('gene_name').set_index('gene_name')
    EX = pd.read_csv(tab('Table_S8e_mean_expression_by_population.csv'), index_col=0); oth = [c for c in EX.columns if c not in ('gene_name', 'Hepatocito')]
    hepspec = set(EX[(EX['Hepatocito'] - EX[oth].max(axis=1)) > 1].gene_name.astype(str))
    # ---- C2/C3 within-population stage effects and concordance with bulk
    TT = {}; conc = []
    for t in types:
        for inc in (False, True):
            Td, Cd, Pd = donor_matrix(D, recs, genes, t, include_end_stage=inc)
            if len(Td) < 12: continue
            ts, df = stage_t(Td, Cd); key = f'{t}' + (' (+end stage)' if inc else ''); TT[key] = ts
            j = pd.DataFrame({'sn': ts}).join(T6[['t_stage', 't_stage_comp_adjusted']], how='inner')
            jh = j[j.index.isin(hepspec)]
            conc.append({'population': key, 'donors': len(Td), 'genes': len(ts), 'rho_all_vs_bulk_adjusted': stats.spearmanr(j.sn, j.t_stage_comp_adjusted)[0],
                         'rho_all_vs_bulk_raw': stats.spearmanr(j.sn, j.t_stage)[0], 'n_hepspec': len(jh),
                         'rho_hepspec_vs_bulk_adjusted': stats.spearmanr(jh.sn, jh.t_stage_comp_adjusted)[0] if len(jh) > 20 else np.nan})
    pd.DataFrame(TT).round(4).to_csv(tab('Table_S18b_snRNA_within_population_stage_t.csv'))
    CO = pd.DataFrame(conc); CO.round(4).to_csv(tab('Table_S18c_snRNA_concordance_with_bulk.csv'), index=False); print(); print(CO.round(3).to_string(index=False))
    # ---- C4 hepatocyte class and drug receptors in hepatocyte nuclei
    H9 = pd.read_csv(tab('Table_S9c_hepatocyte_intrinsic_genes.csv'), index_col=0); R9 = pd.read_csv(tab('Table_S9d_hepatocyte_drug_targets_reference.csv'), index_col=0)
    th = TT['Hepatocyte']; H9['t_snRNA_hepatocyte'] = th.reindex(H9.gene_name.values).values; R9['t_snRNA_hepatocyte'] = th.reindex(R9.gene_name.values).values
    h = H9.dropna(subset=['t_snRNA_hepatocyte']); agree = (np.sign(h.t_snRNA_hepatocyte) == np.sign(h.t_bulk_comp_adjusted))
    print('\nhepatocyte class: sign agreement %d/%d (binomial P = %.3g); rho = %.2f' % (agree.sum(), len(h), stats.binomtest(int(agree.sum()), len(h), 0.5, alternative='greater').pvalue, stats.spearmanr(h.t_snRNA_hepatocyte, h.t_bulk_comp_adjusted)[0]))
    H9.round(3).to_csv(tab('Table_S18d_snRNA_hepatocyte_class.csv')); R9.round(3).to_csv(tab('Table_S18e_snRNA_hepatocyte_drug_receptors.csv'))
    print(R9[['gene_name', 't_bulk_stage', 't_bulk_comp_adjusted', 't_snRNA_hepatocyte']].round(2).to_string(index=False))
    # ---- C5 neighbour coupling within populations
    X0, keep = pd.read_pickle(inter('expr.pkl')); order = keep.sort_values(['chr', 'grid_index']).assign(gene_name=lambda d: d.gene_name.astype(str))
    rng = np.random.default_rng(0); rows = []
    for key, ts in TT.items():
        o = order[order.gene_name.isin(ts.index)].copy(); o['t'] = ts.reindex(o.gene_name).values
        byc = [o[o.chr == c].t.values for c in CHR]; byc = [v for v in byc if len(v) > 20]
        a1 = acf_by_chr(byc, (1,))[0]; nul = np.array([acf_by_chr([rng.permutation(v) for v in byc], (1,))[0] for _ in range(200)])
        rows.append({'population': key, 'genes': len(o), 'acf_lag1': a1, 'null_p97.5': np.quantile(nul, .975), 'z': (a1 - nul.mean()) / nul.std()})
    AC = pd.DataFrame(rows); AC.round(4).to_csv(tab('Table_S18f_snRNA_neighbour_coupling.csv'), index=False); print(); print(AC.round(3).to_string(index=False))
    # ---- C6 the F4 switch within populations (programme scores and the 291 switch genes)
    SW = pd.read_csv(tab('Table_S17a_F4_switch_per_gene.csv'), index_col=0); swg = set(SW[SW.switch_gene].gene_name.astype(str)); rows = []
    for t in types:
        for inc in (False, True):
            Td, Cd, Pd = donor_matrix(D, recs, genes, t, include_end_stage=inc)
            if len(Td) < 12: continue
            cpm = np.log2(Cd / Cd.sum(0) * 1e6 + 1); ok = (Cd >= 5).mean(1) >= 0.5; z = cpm[ok]; z = z.sub(z.mean(1), axis=0).div(z.std(1) + 1e-9, axis=0)
            Pd = Pd.assign(switch_genes=z[z.index.isin(swg)].mean(0).reindex(Pd.index).values)
            o = Td.orden.values.astype(float); stp = (o >= 5).astype(float)
            for k in Pd.columns:
                y = Pd[k].values.astype(float); X1 = np.column_stack([np.ones_like(o), o]); X2 = np.column_stack([X1, stp])
                b2, res2, *_ = np.linalg.lstsq(X2, y, rcond=None); r2 = y - X2 @ b2; df2 = len(y) - 3
                se = np.sqrt((r2 ** 2).sum() / df2 * np.linalg.pinv(X2.T @ X2)[2, 2])
                rows.append({'population': t + (' (+end stage)' if inc else ''), 'score': k, 'donors': len(y), 'F4_donors': int(stp.sum()), 'step_F4': b2[2], 't_step_F4': b2[2] / se,
                             **{f'mean_{s}': float(np.mean(y[Td.stage.values == s])) if (Td.stage.values == s).any() else np.nan for s in ['Normal', 'F0', 'F1', 'F2', 'F3', 'F4']}})
    SWt = pd.DataFrame(rows); SWt.round(4).to_csv(tab('Table_S18g_snRNA_F4_switch_by_population.csv'), index=False)
    print(); print(SWt[SWt.score.isin(['YAP_TAZ_targets', 'switch_genes'])][['population', 'score', 'donors', 'F4_donors', 't_step_F4', 'mean_F3', 'mean_F4']].round(2).to_string(index=False))

def step_d():
    """Supplementary Figure: validation in snRNA-seq."""
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, scienceplots
    plt.style.use(['science', 'nature', 'no-latex'])
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial'], 'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                         'legend.fontsize': 6.5, 'axes.linewidth': 0.5, 'xtick.top': False, 'ytick.right': False, 'xtick.minor.visible': False, 'ytick.minor.visible': False,
                         'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False, 'pdf.fonttype': 42})
    OI = {'blue': '#0072B2', 'orange': '#E69F00', 'green': '#009E73', 'red': '#D55E00', 'purple': '#CC79A7', 'sky': '#56B4E9', 'grey': '#7F7F7F'}
    def lab(ax, t, dx=-0.22, dy=1.14): ax.text(dx, dy, t, transform=ax.transAxes, fontsize=8, fontweight='bold', va='top')
    CO = pd.read_csv(tab('Table_S18c_snRNA_concordance_with_bulk.csv')); AC = pd.read_csv(tab('Table_S18f_snRNA_neighbour_coupling.csv'))
    H9 = pd.read_csv(tab('Table_S18d_snRNA_hepatocyte_class.csv'), index_col=0); R9 = pd.read_csv(tab('Table_S18e_snRNA_hepatocyte_drug_receptors.csv'), index_col=0)
    SWt = pd.read_csv(tab('Table_S18g_snRNA_F4_switch_by_population.csv'))
    W = 183 / 25.4; fig = plt.figure(figsize=(W, W * 0.62)); gs = fig.add_gridspec(2, 3, hspace=0.72, wspace=0.55)
    pops = ['Hepatocyte', 'Cholangiocyte', 'HSC', 'Endothelium', 'Macrophage', 'T/NK', 'B/Plasma']
    # A concordance with the bulk composition-adjusted effect
    ax = fig.add_subplot(gs[0, 0]); c1 = CO.set_index('population').reindex(pops); c2 = CO.set_index('population').reindex([p + ' (+end stage)' for p in pops])
    y = np.arange(len(pops)); ax.barh(y + 0.19, c1.rho_hepspec_vs_bulk_adjusted, 0.38, color=[OI['green'] if p == 'Hepatocyte' else '0.6' for p in pops], lw=0, label='Biopsies')
    ax.barh(y - 0.19, c2.rho_hepspec_vs_bulk_adjusted.values, 0.38, color=[OI['green'] if p == 'Hepatocyte' else '0.8' for p in pops], lw=0, alpha=0.55, label='+ end-stage explants')
    ax.set_yticks(y); ax.set_yticklabels(pops); ax.invert_yaxis(); ax.axvline(0, color='k', lw=0.5)
    ax.set_xlabel('ρ with bulk intrinsic effect\n(hepatocyte-specific genes)'); ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, columnspacing=1.0); lab(ax, 'A', -0.42, 1.2)
    # B hepatocyte class: bulk intrinsic vs snRNA hepatocyte nuclei
    ax = fig.add_subplot(gs[0, 1]); h = H9.dropna(subset=['t_snRNA_hepatocyte'])
    ax.axhline(0, color='0.75', lw=0.4); ax.axvline(0, color='0.75', lw=0.4)
    ax.scatter(h.t_bulk_comp_adjusted, h.t_snRNA_hepatocyte, s=11, color=OI['green'], lw=0.3, edgecolor='white', label='Hepatocyte-specific genes')
    r = R9.dropna(subset=['t_snRNA_hepatocyte']); ax.scatter(r.t_bulk_comp_adjusted, r.t_snRNA_hepatocyte, s=14, marker='D', color=OI['orange'], lw=0.3, edgecolor='k', label='Drug receptors')
    offs = {'THRB': (-22, -9), 'NR1H4': (4, -9), 'KLB': (-14, 5), 'SCD': (4, 4), 'PNPLA3': (4, 4), 'HSD17B13': (4, -8)}
    for _, q in r.iterrows():
        if q.gene_name in offs: ax.annotate(q.gene_name, (q.t_bulk_comp_adjusted, q.t_snRNA_hepatocyte), xytext=offs[q.gene_name], textcoords='offset points', fontsize=6, style='italic')
    for _, q in h.iterrows():
        if q.gene_name in ('FGF21', 'MAT1A', 'HAAO', 'CHI3L1'): ax.annotate(q.gene_name, (q.t_bulk_comp_adjusted, q.t_snRNA_hepatocyte), xytext=(3, -7), textcoords='offset points', fontsize=6, style='italic')
    agree = (np.sign(h.t_snRNA_hepatocyte) == np.sign(h.t_bulk_comp_adjusted)).sum()
    ax.text(0.03, 0.97, f'sign agreement {agree}/{len(h)}', transform=ax.transAxes, va='top', fontsize=6.5)
    ax.set_xlabel('Bulk, composition-adjusted (t)'); ax.set_ylabel('snRNA-seq hepatocyte nuclei (t)'); ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, handletextpad=0.2, columnspacing=1.0); lab(ax, 'B', -0.22, 1.2)
    # C neighbour coupling within populations
    ax = fig.add_subplot(gs[0, 2]); a1 = AC.set_index('population').reindex(pops)
    ax.barh(y, a1.acf_lag1, color=[OI['green'] if p == 'Hepatocyte' else OI['red'] for p in pops], lw=0, height=0.62)
    ax.errorbar([0] * len(pops), y, xerr=a1['null_p97.5'], fmt='none', color='k', capsize=2, lw=0.6)
    ax.set_yticks(y); ax.set_yticklabels(pops); ax.invert_yaxis(); ax.set_xlabel('Neighbour autocorrelation of\nthe stage effect (biopsies)'); lab(ax, 'C', -0.42, 1.2)
    # D–F the F4 switch within populations
    for k, (pop, L) in enumerate([('Hepatocyte', 'D'), ('Endothelium', 'E'), ('Cholangiocyte', 'F')]):
        ax = fig.add_subplot(gs[1, k]); rr = SWt[(SWt.population == pop + ' (+end stage)') & (SWt.score == 'switch_genes')].iloc[0]
        st = ['Normal', 'F0', 'F1', 'F2', 'F3', 'F4']; vals = [rr[f'mean_{s}'] for s in st]
        ax.plot(range(6), vals, 'o-', color=OI['red'], ms=3, lw=1.0)
        rb = SWt[(SWt.population == pop) & (SWt.score == 'switch_genes')]
        ax.set_xticks(range(6)); ax.set_xticklabels(st); ax.axhline(0, color='0.8', lw=0.4); ax.set_xlabel('Fibrosis stage'); ax.set_ylabel('F4 switch-gene score\n(mean z, pseudobulk)')
        txt = f'step at F4: t = {rr.t_step_F4:.1f} (with explants)'
        if len(rb): txt += f'\n{rb.iloc[0].t_step_F4:.1f} (biopsies only, {int(rb.iloc[0].F4_donors)} F4)'
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va='top', fontsize=6.2); ax.set_title(pop, fontsize=7, pad=2); lab(ax, L, -0.3, 1.18)
    fig.savefig(os.path.join(FIG, 'jhep', 'SuppFig_snRNAseq_validation.pdf'), bbox_inches='tight'); fig.savefig(os.path.join(FIG, 'jhep', 'SuppFig_snRNAseq_validation.png'), dpi=400, bbox_inches='tight')
    fig.savefig(os.path.join(FIG, 'jhep', 'SuppFig_snRNAseq_validation.tiff'), dpi=300, bbox_inches='tight', pil_kwargs={'compression': 'tiff_lzw'}); print('figure ok')

def step_e():
    """Neighbour co-expression across single nuclei, per sample × population (≥150 nuclei), pooled per donor."""
    out = tab('Table_S18h_snRNA_neighbour_coexpression.csv')
    D = meta(); X0, keep = pd.read_pickle(inter('expr.pkl'))
    order = keep.sort_values(['chr', 'grid_index']).assign(gene_name=lambda d: d.gene_name.astype(str))
    pairs = []
    for c, g in order.groupby('chr'):
        n = g.gene_name.values; pairs += list(zip(n[:-1], n[1:]))
    allS = np.hstack([np.load(f)['S'] for f in sorted(glob.glob(os.path.join(CACHE, '*_scores.npz')))]); q99 = np.percentile(allS, 99, axis=1, keepdims=True)
    types = np.array(list(MARK)); rng = np.random.default_rng(0)
    part = os.path.join(CACHE, 'coexpr_partial.csv'); done = set(pd.read_csv(part).gsm) if os.path.exists(part) else set()
    for f in sorted(glob.glob(os.path.join(DD, 'GSM*raw_counts_csv.gz'))):
        gsm = os.path.basename(f).split('_')[0]
        if gsm in done: continue
        X, genes, L, gi = _read(f); Z = np.load(os.path.join(CACHE, gsm + '_scores.npz'))['S'] / q99
        srt = np.sort(Z, 0); ct = types[Z.argmax(0)]; ct[(srt[-1] < 0.3) | ((srt[-1] - srt[-2]) < 0.1)] = 'Unassigned'
        res = []
        for t in MARK:
            m = ct == t
            if m.sum() < 150: continue
            A = L[:, m]; ex = (A > 0).mean(1) >= 0.10; sub = genes[ex]; pos = {g: j for j, g in enumerate(sub)}
            if len(sub) < 300: continue
            A = A[ex]; A = (A - A.mean(1, keepdims=True)) / (A.std(1, keepdims=True) + 1e-9)
            pp = [(a, b) for a, b in pairs if a in pos and b in pos]
            if len(pp) < 50: continue
            r = np.einsum('ij,ij->i', A[[pos[a] for a, _ in pp]], A[[pos[b] for _, b in pp]]) / A.shape[1]
            i1 = rng.integers(0, len(sub), len(pp)); i2 = rng.integers(0, len(sub), len(pp)); ok = i1 != i2
            rn = np.einsum('ij,ij->i', A[i1[ok]], A[i2[ok]]) / A.shape[1]
            res.append({'gsm': gsm, 'donor': D.loc[gsm, 'patient id'], 'type': t, 'nuclei': int(m.sum()), 'pairs': len(pp), 'r_neighbours': r.mean(), 'r_random': rn.mean()})
        pd.DataFrame(res if res else [{'gsm': gsm}]).to_csv(part, mode='a', header=not os.path.exists(part), index=False); print(gsm, len(res), flush=True)
    P = pd.read_csv(part).dropna(subset=['type'])
    Dn = P.groupby(['donor', 'type'])[['r_neighbours', 'r_random']].mean().reset_index()
    S = Dn.groupby('type').agg(donors=('donor', 'nunique'), r_neighbours=('r_neighbours', 'mean'), r_random=('r_random', 'mean')).reset_index()
    S['ratio'] = S.r_neighbours / S.r_random
    S['wilcoxon_P'] = [stats.wilcoxon(g.r_neighbours, g.r_random).pvalue if len(g) >= 6 else np.nan for _, g in Dn.groupby('type')]
    S.round(5).to_csv(out, index=False); print(S.round(4).to_string(index=False))

if __name__ == '__main__':
    stp = sys.argv[1] if len(sys.argv) > 1 else 'ABCDE'
    if 'A' in stp: step_a()
    if 'B' in stp: step_b()
    if 'C' in stp: step_c()
    if 'D' in stp: step_d()
    if 'E' in stp: step_e()
