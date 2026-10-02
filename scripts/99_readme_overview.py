"""99 — Overview figure for the README (assets/overview.png), drawn from the released tables."""
import os, numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from common import *
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial'], 'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.linewidth': 0.7, 'legend.frameon': False})
OI = {'blue': '#0072B2', 'orange': '#E69F00', 'red': '#D55E00', 'grey': '#7F7F7F', 'green': '#009E73'}
ST = ['Normal', 'F0', 'F1', 'F2', 'F3', 'F4']
D = pd.read_csv(tab('Table_S7c_F4_likeness_per_sample.csv'), index_col=0); D = D[D.estadio.isin(ST)]
fig, axs = plt.subplots(1, 2, figsize=(10, 3.4), gridspec_kw={'width_ratios': [1.35, 1]}); plt.subplots_adjust(wspace=0.32)
ax = axs[0]
g = D.groupby('estadio')['F4_likeness'].agg(['mean', 'sem']).reindex(ST); lo, hi = g['mean'].min(), g['mean'].max()
ax.errorbar(range(6), (g['mean'] - lo) / (hi - lo), yerr=g['sem'] / (hi - lo), fmt='o-', color=OI['blue'], lw=1.8, ms=5, capsize=2, label='Cell replacement (composition)')
# YAP/TAZ programme beyond composition: residual of a fit on composition scores, estimated in Normal-F3 only
cc = ['Hepatocito', 'HSC', 'Colangiocito']; fit = D[D.estadio != 'F4']
X = np.column_stack([np.ones(len(fit))] + [fit[c].values for c in cc]); b = np.linalg.lstsq(X, fit.YAP_TAZ_targets.values, rcond=None)[0]
D['yap_beyond'] = D.YAP_TAZ_targets - np.column_stack([np.ones(len(D))] + [D[c].values for c in cc]) @ b
g = D.groupby('estadio')['yap_beyond'].agg(['mean', 'sem']).reindex(ST); sc_ = g['mean'].max() - g['mean'].min()
ax.errorbar(range(6), (g['mean'] - g['mean'].iloc[:5].mean()) / sc_, yerr=g['sem'] / sc_, fmt='o-', color=OI['red'], lw=1.8, ms=5, capsize=2, label='YAP/TAZ programme beyond composition')
ax.axvspan(4.5, 5.5, color=OI['red'], alpha=0.08, lw=0); ax.text(5, -0.12, 'switch', ha='center', color=OI['red'], fontsize=9)
ax.axhline(0, color='0.85', lw=0.6)
ax.set_xticks(range(6)); ax.set_xticklabels(ST); ax.set_ylabel('Scaled score'); ax.set_xlabel('Fibrosis stage (437 biopsies)')
ax.set_title('Two layers: linear cell replacement, then a switch at F4', fontsize=10, loc='left'); ax.legend(loc='upper left', fontsize=8.5)
ax = axs[1]; thr = D[D.estadio == 'F4'].F4_likeness.quantile(0.25); f3 = D[D.estadio == 'F3']
grp = {'F3\nother': f3[f3.F4_likeness <= thr], 'F3 with F4-like\ncomposition': f3[f3.F4_likeness > thr], 'F4': D[D.estadio == 'F4']}
rng = np.random.default_rng(1)
for i, (k, d) in enumerate(grp.items()):
    v = d.YAP_TAZ_targets.values; ax.scatter(i + rng.uniform(-0.17, 0.17, len(v)), v, s=9, color=[OI['grey'], OI['blue'], OI['red']][i], alpha=0.7, lw=0)
    ax.errorbar(i, v.mean(), yerr=v.std() / np.sqrt(len(v)), fmt='_', color='k', ms=18, mew=1.6, elinewidth=1.2)
ax.set_xticks(range(3)); ax.set_xticklabels(list(grp)); ax.set_ylabel('YAP/TAZ target score'); ax.axhline(0, color='0.8', lw=0.6)
ax.set_title('Composition does not define cirrhosis; the switch does', fontsize=10, loc='left')
os.makedirs(os.path.join(ROOT, 'assets'), exist_ok=True); fig.savefig(os.path.join(ROOT, 'assets', 'overview.png'), dpi=200, bbox_inches='tight'); print('assets/overview.png')
