"""15 — MASLD GWAS loci and cis-coupled neighbourhoods. Curated list of replicated MASLD/NAFLD loci (candidate gene per
locus); replace/extend with the GWAS Catalog export (data/external/gwas_catalog_MASLD.tsv) when available.
Tests: (i) stage effect of each GWAS gene and its ±2 neighbours; (ii) enrichment of GWAS genes in coupled neighbour pairs;
(iii) spatial autocorrelation of the stage effect in ±5-gene windows around loci vs random windows; (iv) shared eQTLs.
Outputs: Table_S16 (curated loci), S16b (curated-locus enrichment), S16c (GWAS Catalog enrichment by trait), S16d (per-locus table)
"""
import numpy as np, pandas as pd
from scipy import stats
from common import *
X, keep = pd.read_pickle(inter('expr.pkl')); T = pd.read_csv(tab('Table_S6c_stage_effect_per_gene_with_without_composition.csv')).set_index('gene_id')
P = pd.read_csv(tab('Table_S5f_adjacent_pairs_distance_orientation_TAD.csv')); V = pd.read_csv(tab('Table_S5c_contiguous_DE_neighbourhoods.csv')); Cn = pd.read_csv(tab('Table_S9a_candidate_targets_lineage_intrinsic.csv'))
E = pd.read_csv(tab('Table_S13d_adjacent_pairs_shared_GTEx_signif_variants.csv'))
GWAS = {'PNPLA3': 'rs738409', 'TM6SF2': 'rs58542926', 'GCKR': 'rs1260326', 'MBOAT7': 'rs641738', 'TMC4': 'rs641738', 'HSD17B13': 'rs72613567', 'MTARC1': 'rs2642438', 'GPAM': 'rs2792751', 'APOE': 'rs429358', 'TRIB1': 'rs2954021', 'LYPLAL1': 'rs12137855', 'PPP1R3B': 'rs4240624', 'COBLL1': 'rs12052337', 'MERTK': 'rs4374383', 'SUGP1': 'rs10401969', 'INHBE': 'rs35781401', 'CIDEB': 'rs6063', 'APOH': 'rs1801689', 'FTO': 'rs1421085', 'ADH1B': 'rs1229984', 'SERPINA1': 'rs28929474', 'TOR1B': 'rs7029757', 'LEPR': 'rs1137101', 'MTTP': 'rs3816873', 'ERLIN1': 'rs2862954', 'GATAD2A': 'rs4808199', 'HFE': 'rs1800562', 'TMEM106B': 'rs1990622', 'HNF1A': 'rs1169288'}
names = keep.set_index('gene_id').gene_name.astype(str); g2id = {n: i for i, n in names.items()}
order = keep.sort_values(['chr', 'grid_index']).copy(); order['t'] = T.t_stage.reindex(order.gene_id).values; order['tc'] = T.t_stage_comp_adjusted.reindex(order.gene_id).values; order['pos'] = np.arange(len(order))
present = [g for g in GWAS if g in g2id]; oi = order.set_index('gene_name'); rows = []
for g in present:
    p = int(oi.loc[g, 'pos']); c = oi.loc[g, 'chr']; w = order[(order.pos >= p - 2) & (order.pos <= p + 2) & (order.chr == c)]; nb = w[w.gene_name != g]
    rows.append({'gene': g, 'snp': GWAS[g], 'chr': c, 't_stage': oi.loc[g, 't'], 't_comp_adj': oi.loc[g, 'tc'], 'neighbours': ','.join(nb.gene_name.astype(str)), 'neighbour_t': ','.join(f'{x:.1f}' for x in nb.t),
                 'n_neighbours_concordant': int(((np.sign(nb.t) == np.sign(oi.loc[g, 't'])) & (nb.t.abs() > 2)).sum()), 'in_DE_neighbourhood': any(g in str(s).split(',') for s in V.genes), 'lineage_intrinsic_target': g in set(Cn.gene_name)})
G = pd.DataFrame(rows).sort_values('t_stage'); G.round(3).to_csv(tab('Table_S16_GWAS_loci_neighbourhoods.csv'), index=False)
P['coupled'] = (np.sign(P.t1) == np.sign(P.t2)) & (P.t1.abs() > 3) & (P.t2.abs() > 3); cg = set(P[P.coupled].n1) | set(P[P.coupled].n2); allg = set(P.n1) | set(P.n2)
a = sum(g in cg for g in present); n = sum(g in allg for g in present); base = len(cg) / len(allg)
print(f'GWAS genes in a coupled pair: {a}/{n} ({100*a/n:.0f}%) vs {100*base:.0f}% background; binomial P={stats.binomtest(a, n, base, alternative="greater").pvalue:.3f}')
def acf_win(centers, L=1):
    num = den = 0
    for p in centers:
        w = order[(order.pos >= p - 5) & (order.pos <= p + 5)]; v = w.t.fillna(0).values - w.t.fillna(0).mean(); num += np.sum(v[:-L] * v[L:]); den += np.sum(v * v) * (len(v) - L) / len(v)
    return num / den
cen = [int(oi.loc[g, 'pos']) for g in present]; obs = acf_win(cen); rng = np.random.default_rng(0); nl = np.array([acf_win(rng.integers(5, len(order) - 5, len(cen))) for _ in range(500)])
print(f'lag-1 ACF around GWAS loci {obs:.3f} vs random windows {nl.mean():.3f}±{nl.std():.3f} (P={np.mean(nl >= obs):.3f})')
pd.DataFrame({'acf_gwas_windows': [obs], 'null_mean': [nl.mean()], 'null_sd': [nl.std()], 'p': [np.mean(nl >= obs)], 'frac_gwas_in_coupled_pair': [a / n], 'background': [base]}).to_csv(tab('Table_S16b_GWAS_enrichment_tests.csv'), index=False)

# ---- GWAS Catalog version (requires data/external/gwas-catalog-download-associations-v1.0-full.tsv; genome-wide significant P<5e-8)
cat = os.path.join(EXT, 'gwas-catalog-download-associations-v1.0-full.tsv')
if not os.path.exists(cat):
    raise FileNotFoundError(f'{cat} is required for Supplementary Table S16c-d and Supplementary Fig. S7; '
                            'download "All associations v1.0" from https://www.ebi.ac.uk/gwas/docs/file-downloads')
cols = ['DISEASE/TRAIT', 'CHR_ID', 'CHR_POS', 'MAPPED_GENE', 'SNP_GENE_IDS', 'UPSTREAM_GENE_ID', 'DOWNSTREAM_GENE_ID', 'SNPS', 'P-VALUE', 'PUBMEDID', 'INTERGENIC']
Gc = pd.read_csv(cat, sep='\t', usecols=cols, low_memory=False, encoding='utf-8'); Gc['P'] = pd.to_numeric(Gc['P-VALUE'], errors='coerce'); Gc = Gc[Gc.P < 5e-8]
t = Gc['DISEASE/TRAIT'].str.lower()
traits = {'MASLD/NAFLD': t.str.contains(r'non-?alcoholic fatty liver|nafld|steatotic liver|metabolic dysfunction.associated steatotic|fatty liver|hepatic steatosis|liver fat', regex=True),
          'liver fibrosis/cirrhosis': t.str.contains(r'liver fibrosis|cirrhosis|fibrosis.*liver|liver stiffness|nash|steatohepatitis', regex=True) & ~t.str.contains('primary biliary|sclerosing', regex=True),
          'ALT (liver enzyme)': t.str.contains(r'alanine aminotransferase', regex=True),
          'height (control)': t.str.fullmatch(r'height|body height|standing height'),
          'educational attainment (control)': t.str.contains(r'educational attainment', regex=True)}
Tg = pd.read_csv(tab('Table_S6c_stage_effect_per_gene_with_without_composition.csv')).set_index('gene_id')
Pp = pd.read_csv(tab('Table_S5f_adjacent_pairs_distance_orientation_TAD.csv')); Vn = pd.read_csv(tab('Table_S5c_contiguous_DE_neighbourhoods.csv'))
od = keep.sort_values(['chr', 'grid_index']).reset_index(drop=True); od['t'] = Tg.t_stage.reindex(od.gene_id).fillna(0).values; od['pos'] = np.arange(len(od))
gid2pos = dict(zip(od.gene_id, od.pos)); name2id = dict(zip(od.gene_name.astype(str), od.gene_id)); tv = od.t.values
def locus_genes(df):
    ids = set()
    for g, u, d in zip(df.SNP_GENE_IDS.fillna('').astype(str), df.UPSTREAM_GENE_ID.fillna('').astype(str), df.DOWNSTREAM_GENE_ID.fillna('').astype(str)):
        if g.strip(): ids |= set(x.strip() for x in g.split(','))
        else:
            for v in (u, d):
                if v.strip(): ids.add(v.strip())
    return {i for i in ids if i in gid2pos}
Pp['coupled'] = (np.sign(Pp.t1) == np.sign(Pp.t2)) & (Pp.t1.abs() > 3) & (Pp.t2.abs() > 3)
cg = set(Pp[Pp.coupled].g1) | set(Pp[Pp.coupled].g2); allg = set(Pp.g1) | set(Pp.g2); base = len(cg) / len(allg)
run_genes = set(name2id[g] for s_ in Vn.genes for g in str(s_).split(',') if g in name2id); base_run = len(run_genes) / len(od)
def acf_win(centers, half=5):
    num = den = 0.0
    for p in centers:
        lo, hi = max(0, p - half), min(len(tv), p + half + 1); v = tv[lo:hi]; v = v - v.mean(); num += np.sum(v[:-1] * v[1:]); den += np.sum(v * v) * (len(v) - 1) / len(v)
    return num / den
rng = np.random.default_rng(0); rows = []; detail = {}
for k, m in traits.items():
    df = Gc[m]; genes = locus_genes(df); detail[k] = genes; gl = [g for g in genes if g in allg]; a = sum(g in cg for g in gl); r = sum(g in run_genes for g in genes)
    cen = np.array(sorted(gid2pos[g] for g in genes))   # all loci, sorted: deterministic (no subsampling)
    obs = acf_win(cen); nl = np.array([acf_win(rng.integers(5, len(od) - 5, len(cen))) for _ in range(200)])
    nbt = np.concatenate([np.abs(tv[[p - 1, p + 1]]) for p in cen if 0 < p < len(tv) - 1])
    rows.append({'trait_set': k, 'n_assoc': len(df), 'locus_genes_on_grid': len(genes), 'frac_in_coupled_pair': a / max(1, len(gl)), 'bg_coupled': base,
                 'P_coupled': stats.binomtest(a, len(gl), base, alternative='greater').pvalue,
                 'frac_in_DE_neighbourhood': r / len(genes), 'bg_nbhd': base_run, 'P_nbhd': stats.binomtest(r, len(genes), base_run, alternative='greater').pvalue,
                 'acf_windows': obs, 'acf_random': nl.mean(), 'acf_random_sd': nl.std(), 'P_acf': np.mean(nl >= obs), 'mean_|t|_neighbours': nbt.mean(), 'mean_|t|_genome': np.abs(tv).mean()})
R16 = pd.DataFrame(rows); R16.round(4).to_csv(tab('Table_S16c_GWAS_catalog_enrichment_by_trait.csv'), index=False); print(R16.round(3).to_string(index=False))
mg = detail['MASLD/NAFLD'] | detail['liver fibrosis/cirrhosis']; oi = od.set_index('gene_id'); rows = []
for g in mg:
    p = int(oi.loc[g, 'pos']); c = oi.loc[g, 'chr']; w = od[(od.pos >= p - 2) & (od.pos <= p + 2) & (od.chr == c)]; nb = w[w.gene_id != g]
    rows.append({'gene_id': g, 'gene': oi.loc[g, 'gene_name'], 'chr': c, 't_stage': oi.loc[g, 't'], 'neighbours': ','.join(nb.gene_name.astype(str)), 'neighbour_t': ','.join(f'{x:.1f}' for x in nb.t),
                 'in_coupled_pair': g in cg, 'in_DE_neighbourhood': g in run_genes, 'trait': 'MASLD' if g in detail['MASLD/NAFLD'] else 'fibrosis/cirrhosis'})
pd.DataFrame(rows).sort_values('t_stage').round(3).to_csv(tab('Table_S16d_GWAS_catalog_MASLD_fibrosis_loci_neighbourhoods.csv'), index=False)
