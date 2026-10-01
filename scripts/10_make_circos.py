"""10 — Figures 6 (genome circos: DE, composition-adjusted DE, TF targets, neighbourhoods, pathway links) and 7
(chord: intrinsic TFs x Hallmark pathways). Run after scripts 04, 05, 07, 08."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from pycirclize import Circos
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from common import *

OI = {'blue': '#0072B2', 'orange': '#E69F00', 'green': '#009E73', 'red': '#D55E00', 'purple': '#CC79A7', 'sky': '#56B4E9', 'grey': '#7F7F7F'}
X, keep = pd.read_pickle(inter('expr.pkl')); de, de2 = pd.read_pickle(inter('de.pkl')); grid = load_grid(); N = chrom_lengths(grid)
T = pd.read_csv(tab('Table_S6c_stage_effect_per_gene_with_without_composition.csv')).set_index('gene_id'); V = pd.read_csv(tab('Table_S5c_contiguous_DE_neighbourhoods.csv')); V['chr'] = V.chr.astype(str)
D = pd.read_csv(tab('Table_S10c_dorothea_ABC_regulons_used.csv')); G = pd.read_csv(tab('Table_S10a_GSEA_hallmark_ordinal.csv')); Cn = pd.read_csv(tab('Table_S9a_candidate_targets_lineage_intrinsic.csv'))
TF = pd.read_csv(tab('Table_S10d_DoRothEA_TF_activity_vs_stage.csv')); Mx = pd.read_csv(tab('Table_S10h_TF_x_Hallmark_overlap.csv'), index_col=0)
K = keep.set_index('gene_id').join(de2[['stat']]).join(T[['t_stage_comp_adjusted']]); K['stat'] = K.stat.fillna(0); K['gene_name'] = K.gene_name.astype(str)
hubs = ['EPITHELIAL_MESENCHYMAL_TRANSITION', 'TNFA_SIGNALING_VIA_NFKB', 'INFLAMMATORY_RESPONSE', 'ANGIOGENESIS', 'INTERFERON_GAMMA_RESPONSE', 'P53_PATHWAY', 'BILE_ACID_METABOLISM', 'FATTY_ACID_METABOLISM', 'OXIDATIVE_PHOSPHORYLATION']
short = {'EPITHELIAL_MESENCHYMAL_TRANSITION': 'EMT', 'TNFA_SIGNALING_VIA_NFKB': 'TNFα/NF-κB', 'INFLAMMATORY_RESPONSE': 'Inflammation', 'ANGIOGENESIS': 'Angiogenesis', 'INTERFERON_GAMMA_RESPONSE': 'IFN-γ', 'P53_PATHWAY': 'p53', 'BILE_ACID_METABOLISM': 'Bile acid', 'FATTY_ACID_METABOLISM': 'Fatty acid', 'OXIDATIVE_PHOSPHORYLATION': 'OXPHOS'}
lead = {r.Term: set(str(r.Lead_genes).split(';')) for r in G.itertuples()}; nes = dict(zip(G.Term, G.NES)); hcol = dict(zip(hubs, plt.cm.tab10(np.linspace(0, 1, 10))[:len(hubs)]))
links = []
for _, v in V.iterrows():
    genes = set(str(v.genes).split(',')); mid = (v.start_grid + v.end_grid) / 2
    for h in hubs:
        if genes & lead.get(h, set()): links.append({'chr': v.chr, 'pos': mid, 'hub': h})
L = pd.DataFrame(links); L.to_csv(tab('Table_S5j_neighbourhood_pathway_links.csv'), index=False)
tot = sum(N.values()); hubsize = int(0.045 * tot); nfkb = set(D[D.tf.isin(['NFKB1', 'RELA'])].target); tp53 = set(D[D.tf == 'TP53'].target); cand = set(Cn.gene_name)
circos = Circos({**{c: N[c] for c in CHR}, **{short[h]: hubsize for h in hubs}}, space=1.2, start=0, end=350)
for sec in circos.sectors:
    c = sec.name
    if c in CHR:
        sec.text(c, r=109, size=6); k = K[K.chr == c].sort_values('grid_index'); x = k.grid_index.values.astype(float); y = k.stat.values
        tr = sec.add_track((86, 100), r_pad_ratio=0.05); tr.axis(fc='#f2f2f2', ec='none'); tr.bar(x[y > 0], np.clip(y[y > 0], 0, 12), width=1, color=OI['red'], lw=0, vmin=0, vmax=12, bottom=0); tr.bar(x[y < 0], np.clip(-y[y < 0], 0, 12), width=1, color=OI['blue'], lw=0, vmin=0, vmax=12, bottom=0)
        tr2 = sec.add_track((76, 84), r_pad_ratio=0.05); tr2.axis(fc='#fafafa', ec='none'); ya = k.t_stage_comp_adjusted.fillna(0).values
        tr2.bar(x[ya > 0], np.clip(ya[ya > 0], 0, 6), width=1, color='#b05a2c', lw=0, vmin=0, vmax=6, bottom=0); tr2.bar(x[ya < 0], np.clip(-ya[ya < 0], 0, 6), width=1, color='#2c5f8c', lw=0, vmin=0, vmax=6, bottom=0)
        tr3 = sec.add_track((69, 74)); tr3.axis(fc='white', ec='0.85', lw=0.3)
        gi = k[k.gene_name.isin(nfkb)].grid_index.values; tr3.scatter(gi, [0.75] * len(gi), s=0.8, color=OI['orange'], vmin=0, vmax=1, lw=0)
        gi = k[k.gene_name.isin(tp53)].grid_index.values; tr3.scatter(gi, [0.25] * len(gi), s=0.8, color=OI['purple'], vmin=0, vmax=1, lw=0)
        tr4 = sec.add_track((62, 68)); tr4.axis(fc='white', ec='0.85', lw=0.3)
        for _, v in V[V.chr == c].iterrows(): tr4.rect(v.start_grid, v.end_grid + 1, fc=OI['red'] if v.direction == 'up' else OI['blue'], ec='none')
        for g_ in k[k.gene_name.isin(cand)].grid_index.values: tr4.line([g_, g_], [0, 1], color=OI['green'], lw=0.5, vmin=0, vmax=1)
    else:
        h = [hh for hh in hubs if short[hh] == c][0]; sec.text(f"{c}\nNES {nes.get(h, 0):+.1f}", r=112, size=5.5, color=hcol[h], fontweight='bold'); tr = sec.add_track((62, 100)); tr.axis(fc=hcol[h], ec='none', alpha=0.3)
for _, r in L.iterrows():
    hs = short[r.hub]; tsec = circos.get_sector(hs); e = tsec.size * ((CHR.index(r.chr) + 0.5) / len(CHR))
    circos.link((r.chr, r.pos - 1.5, r.pos + 1.5), (hs, e - hubsize * 0.015, e + hubsize * 0.015), color=hcol[r.hub], alpha=0.45, lw=0.25, r1=61, r2=61)
fig = circos.plotfig(figsize=(7.5, 8.2), dpi=300)
handles = [Patch(fc=OI['red'], label='Stage effect per gene, up (outer: DESeq2; middle: composition-adjusted)'), Patch(fc=OI['blue'], label='Stage effect per gene, down'), Line2D([], [], marker='o', ls='', color=OI['orange'], ms=4, label='NF-κB (NFKB1/RELA) target'), Line2D([], [], marker='o', ls='', color=OI['purple'], ms=4, label='TP53 target'), Patch(fc=OI['red'], alpha=.6, label='Contiguous DE neighbourhood (≥3 genes)'), Line2D([], [], color=OI['green'], lw=1.2, label='Lineage-intrinsic candidate target'), Line2D([], [], color='0.5', lw=1, label='Link: neighbourhood → Hallmark pathway (leading edge)')]
fig.legend(handles=handles, loc='lower center', ncol=2, fontsize=6, frameon=False, bbox_to_anchor=(0.5, -0.06))
fig.savefig(f'{FIG}/Fig6_circos_genome_pathways_TF.png', dpi=400, bbox_inches='tight'); fig.savefig(f'{FIG}/Fig6_circos_genome_pathways_TF.pdf', bbox_inches='tight'); plt.close(fig)

# ---- Fig 7
sig = TF[TF.q_adj_comp < 0.05]; up = set(sig[sig.t_stage_comp_adjusted > 0].tf)
shortp = {'EPITHELIAL_MESENCHYMAL_TRANSITION': 'EMT', 'TNFA_SIGNALING_VIA_NFKB': 'TNFα/NF-κB', 'INFLAMMATORY_RESPONSE': 'Inflammatory response', 'ANGIOGENESIS': 'Angiogenesis', 'INTERFERON_GAMMA_RESPONSE': 'IFN-γ response', 'P53_PATHWAY': 'p53 pathway', 'BILE_ACID_METABOLISM': 'Bile acid metabolism', 'FATTY_ACID_METABOLISM': 'Fatty acid metabolism', 'OXIDATIVE_PHOSPHORYLATION': 'OXPHOS', 'ALLOGRAFT_REJECTION': 'Allograft rejection', 'APOPTOSIS': 'Apoptosis', 'IL6_JAK_STAT3_SIGNALING': 'IL6–JAK–STAT3', 'HYPOXIA': 'Hypoxia', 'APICAL_JUNCTION': 'Apical junction', 'KRAS_SIGNALING_UP': 'KRAS up', 'IL2_STAT5_SIGNALING': 'IL2–STAT5', 'G2M_CHECKPOINT': 'G2M checkpoint', 'TGF_BETA_SIGNALING': 'TGF-β', 'E2F_TARGETS': 'E2F targets', 'XENOBIOTIC_METABOLISM': 'Xenobiotic metabolism', 'PEROXISOME': 'Peroxisome', 'COMPLEMENT': 'Complement'}
Mx.index = [shortp.get(p, p) for p in Mx.index]
pc = dict(zip(Mx.index, plt.cm.tab20(np.linspace(0, 1, 20))[:len(Mx.index)])); cmap = {**pc, **{t: (OI['red'] if t in up else OI['blue']) for t in Mx.columns}}
circos2 = Circos.chord_diagram(Mx, space=2.5, cmap=cmap, label_kws=dict(size=6.5, r=107, orientation='vertical'), link_kws=dict(ec='none', alpha=0.55, direction=1))
fig2 = circos2.plotfig(figsize=(7.5, 7.5), dpi=300)
handles = [Patch(fc=OI['red'], label='TF: activity increases with stage\n(DoRothEA A–C; composition-adjusted q<0.05)'), Patch(fc=OI['blue'], label='TF: activity decreases with stage'), Patch(fc='0.6', label='Ribbon width = regulon ∩ pathway leading edge (genes); colour = pathway')]
fig2.legend(handles=handles, loc='upper left', fontsize=6.5, frameon=False, bbox_to_anchor=(0.01, 0.99))
fig2.savefig(f'{FIG}/Fig7_circos_TF_Hallmark.png', dpi=400, bbox_inches='tight'); fig2.savefig(f'{FIG}/Fig7_circos_TF_Hallmark.pdf', bbox_inches='tight'); print('ok')
