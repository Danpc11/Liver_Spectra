"""Shared paths, constants and helpers for the positional-spectra pipeline."""
import os, glob, re
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, 'data', 'raw')
EXT = os.path.join(ROOT, 'data', 'external')
TAB = os.path.join(ROOT, 'results', 'tables')
FIG = os.path.join(ROOT, 'results', 'figures')
INTER = os.path.join(ROOT, 'results', 'intermediate')
for d in (TAB, FIG, INTER, os.path.join(FIG, 'genome_research'), os.path.join(FIG, 'jhep')):
    os.makedirs(d, exist_ok=True)

COHORTS = ['GSE130970', 'GSE135251', 'GSE162694']
FULL = ['GSE130970', 'GSE162694']            # cohorts containing all six conditions
STAGES = ['Normal', 'F0', 'F1', 'F2', 'F3', 'F4']
ORDER_MAP = {'Normal': 0, 'Control': 0, 'F0': 1, 'F1': 2, 'F2': 3, 'F3': 4, 'F4': 5}
CHR = [str(i) for i in range(1, 23)] + ['X']
BANDS = [1, 2.5, 3, 5, 10, 20, 50, 100, 300, 5000]
EULER = 0.5772156649  # -E[log Exp(1)]


def tab(name):
    return os.path.join(TAB, name)


def inter(name):
    return os.path.join(INTER, name)


def load_grid():
    """Positional grid: union of the per-cohort rejilla_genes.tsv files (same gene -> same grid_index)."""
    files = glob.glob(os.path.join(RAW, 'grid', '*_gene_grid.tsv'))
    g = pd.concat([pd.read_csv(f, sep='\t') for f in files]).drop_duplicates('gene_id')
    g['chr'] = g['chr'].astype(str)
    return g


def chrom_lengths(grid):
    return {c: int(grid[grid.chr == c].grid_index.max()) for c in CHR}


def whiten_matrix(P, n):
    """Whiten a (freq x samples) power matrix against a per-chromosome log-log linear fit; unit mean."""
    k = np.arange(1, P.shape[0] + 1)
    lx = np.log10(k / n)
    L = np.log10(P + 1e-12)
    b = np.polyfit(lx, L, 1)
    fit = np.outer(lx, b[0]) + b[1]
    w = P / 10 ** fit
    return w / w.mean(axis=0)


def spectra(mat, keep, N, tapers=None, rng=None):
    """Per-sample whitened positional spectra. mat: genes x samples (index gene_id).
    tapers: optional list of DPSS tapers per chromosome length (dict n -> array K x n)."""
    cols, names = [], []
    for c in CHR:
        g = keep[keep.chr == c].sort_values('grid_index')
        vals = mat.loc[g.gene_id].values
        if rng is not None:
            vals = vals[rng.permutation(len(g))]
        n = N[c]
        S = np.zeros((n, vals.shape[1]))
        S[g.grid_index.values - 1] = vals - vals.mean(axis=0)
        if tapers is None:
            P = np.abs(np.fft.rfft(S, axis=0)[1:n // 2 + 1]) ** 2
        else:
            T = tapers[n]
            P = np.zeros((n // 2, vals.shape[1]))
            for t in T:
                P += np.abs(np.fft.rfft(S * t[:, None], axis=0)[1:n // 2 + 1]) ** 2
            P /= len(T)
        cols.append(whiten_matrix(P, n))
        names += [(c, k) for k in range(1, P.shape[0] + 1)]
    return pd.DataFrame(np.vstack(cols), index=pd.MultiIndex.from_tuples(names, names=['chr', 'k']), columns=mat.columns)


def design(M, s, extra=None):
    """Design matrix: cohort dummies + sex + ordinal stage (+ optional extra columns)."""
    D = pd.get_dummies(M.loc[s, ['cohorte']], drop_first=True).astype(float)
    D['sexoM'] = (M.loc[s, 'sexo'] == 'M').astype(float)
    D['orden'] = M.loc[s, 'orden'].astype(float)
    D.insert(0, 'const', 1.0)
    if extra is not None:
        D = D.join(extra.loc[s])
    return D


def ols_t(Y, D, col='orden'):
    """t statistics of coefficient `col` for every column of Y (samples x features)."""
    Dv = D.values
    XtXi = np.linalg.pinv(Dv.T @ Dv)
    B = XtXi @ Dv.T @ Y
    R = Y - Dv @ B
    df = Y.shape[0] - np.linalg.matrix_rank(Dv)
    j = list(D.columns).index(col)
    se = np.sqrt((R ** 2).sum(0) / df * XtXi[j, j])
    return B[j] / se, df


def acf_by_chr(values_by_chr, lags=(1, 2, 3, 5, 10)):
    out = []
    for L in lags:
        num = den = 0.0
        for v in values_by_chr:
            v = v - v.mean()
            if len(v) > L + 5:
                num += np.sum(v[:-L] * v[L:])
                den += np.sum(v * v) * (len(v) - L) / len(v)
        out.append(num / den)
    return np.array(out)


MARKERS = {
    'Hepatocito': r'^(ALB|APO[ABCEH]|CYP[1-4][A-Z][0-9]|UGT[12]|SERPINA[1-3]|FG[ABG]|TTR|TF|HP|AHSG|HPX|APOH|F2|F9|PCK1|G6PC1?|CPS1|ARG1|OTC|ASS1|HAO1|GLYAT|ACSM[235]|BAAT|SLC22A|SLCO1B|ADH[1-6]|ALDOB|RBP4|AMBP|VTN|ITIH|KNG1|HRG|AGXT|TAT|HPD|FAH|GLS2|CYP8B1|SLC10A1|ABCB11)',
    'HSC': r'^(COL1A[12]|COL3A1|COL5A[12]|COL6A[1-3]|ACTA2|TAGLN|PDGFR[AB]|LUM|DCN|BGN|VCAN|POSTN|SPARC|FN1|LOX|LOXL[1-4]|TIMP[1-3]|MMP2|CTHRC1|FAP|LRRC15|RGS5|RELN|LRAT|PTH1R|HHIP|COLEC1[01]|VIPR1|IGFBP[3567]|TGFB1|THBS[12]|ELN|FBLN[125]|MFAP[45]|MGP|CCN2|SFRP[14])',
    'Colangiocito': r'^(KRT7|KRT19|KRT8|KRT18|KRT23|EPCAM|SOX9|SPP1|CFTR|HNF1B|ONECUT1|PROM1|CD24|TACSTD2|CLDN4|CLDN10|MUC1|MUC5B|ANXA4|DCDC2|FGFR2|JAG1|SLC4A2|AQP1|TFF[123]|LGALS4|BICC1|PKHD1|SCTR|GGT1|DEFB1|PIGR|CXCL6)$',
    'Macrofago': r'^(CD68|CD163|CD14|MARCO|VSIG4|CLEC4F|TIMD4|C1Q[ABC]|LYZ|AIF1|TYROBP|FCER1G|FCGR[1-3][AB]?|CSF1R|MS4A[467]A?|TREM2|GPNMB|MPEG1|CD5L|HMOX1|SLC40A1|CTSB|CTSD|CTSS|LAPTM5|S100A[89]|LST1|FCN1|ITGAM|ITGAX|ADGRE1|CD86|MRC1|MSR1|SIGLEC1|LILR[AB][1-6])$',
    'Linfocito': r'^(CD3[DEG]|CD2|CD247|CD8[AB]|CD4|TRAC|TRBC[12]|IL7R|CCR7|LTB|CD79[AB]|MS4A1|CD19|CD22|IGH[AGMD]|IGKC|IGLC[1-7]|JCHAIN|MZB1|XBP1|GZM[ABHK]|NKG7|GNLY|KLR[BDFK][1-3]?|PRF1|CST7|CTSW|PTPRC|CORO1A|RAC2|LCK|ZAP70|ITK|CD69|CXCR4|SELL|CCL5|IL32|TRDC|TRGC[12])$',
    'Endotelio': r'^(PECAM1|VWF|CDH5|KDR|FLT1|TEK|TIE1|CLEC4G|CLEC4M|CLEC1B|STAB[12]|LYVE1|OIT3|FCN[23]|DNASE1L3|CRHBP|GPR182|LIFR|ENG|PLVAP|ESAM|EGFL7|RAMP2|ROBO4|EMCN|CD34|SOX18|ERG|FLT4|ACKR1|SELE|SELP|ICAM2|MMRN[12]|CLDN5|GJA[45]|TM4SF1|ADGRL4|SPARCL1|ID1|EHD3|CAV1)$',
}

PROGRAMMES = {
    'ECM/fibrosis': r'^(COL[0-9]+A|MMP[0-9]|TIMP[1-4]|LOX(L[1-4])?$|SPARC|FN1$|ACTA2|PDGFR[AB]|TGFB[1-3]$|TGFBI|THBS[1-4]|LUM$|DCN$|BGN$|VCAN|POSTN|CCN2|ELN$|FBN[12]|LAM[ABC][0-9]|ITGA|ITGB|VIM$|TNC$|PLOD|P4HA|FAP$|LRRC15|CTHRC1)',
    'Hepatocito': r'^(CYP[0-9]|APO[A-Z]|ALB$|UGT[12]|SLC22A|SLCO1|ADH[1-7]|ALDH|HSD|AKR1|SERPINA|FG[ABG]$|TTR$|TF$|HP$|APOH|AHSG|F2$|F9$|PCK1|G6PC|CPS1|ARG1|OTC$|ASS1|HAO|GLYAT|ACSM|BAAT)',
    'Inmune': r'^(CD[0-9]+[A-Z]?$|HLA-|CXCL|CCL[0-9]|CCR[0-9]|CXCR|IL[0-9]+[A-Z]?$|IL[0-9]+R|TNF|PTPRC|LCP|FCGR|FCER|LYZ$|C1Q|CSF[1-3]R|TLR[0-9]|IRF|STAT[1-6]|IFI|ISG|GBP[1-7]|MS4A|TREM|TYROBP|AIF1|S100A[89]|CTSS|LAPTM5|CORO1A|RAC2|B2M$)',
    'Colangiocito/progenitor': r'^(KRT7$|KRT19$|KRT23|EPCAM|SOX9|SPP1$|CFTR|HNF1B|ONECUT|PROM1|CD24$|TACSTD2|CLDN4|MUC1|ANXA4|DCDC2|FGFR2|JAG1|NOTCH2)',
    'Ciclo celular': r'^(MKI67|TOP2A|CCN[ABDE][0-9]|CDK[0-9]$|CDC[0-9]|BUB1|AURK|PLK[1-4]|KIF[0-9]|CENP|MCM[0-9]|PCNA|TYMS|RRM[12]|E2F|FOXM1|BIRC5|UBE2C|NUSAP1|TPX2)',
    'Ribosoma': r'^(RP[LS][0-9]|EEF|EIF[1-6]|RPLP|RPSA)',
}
