"""09 — Figures 1–5 (Nature style). Run after scripts 01–08."""
import pandas as pd, numpy as np, glob, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, scienceplots, seaborn as sns
from matplotlib.lines import Line2D
from adjustText import adjust_text
plt.style.use(['science','nature','no-latex'])
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Liberation Sans','Arial','Helvetica'],'font.size':7,'axes.titlesize':7.5,'axes.labelsize':7,
 'xtick.labelsize':6.5,'ytick.labelsize':6.5,'legend.fontsize':6,'axes.linewidth':0.6,'xtick.major.width':0.6,'ytick.major.width':0.6,
 'xtick.direction':'out','ytick.direction':'out','xtick.top':False,'ytick.right':False,'xtick.minor.visible':False,'ytick.minor.visible':False,
 'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,'pdf.fonttype':42,'ps.fonttype':42,'figure.dpi':150})
MM=1/25.4; W2=183*MM
from common import *
FO=FIG; OUT=TAB
ST=['Normal','F0','F1','F2','F3','F4']; stage_pal=dict(zip(ST,sns.color_palette('viridis',6)))
OI={'blue':'#0072B2','orange':'#E69F00','green':'#009E73','red':'#D55E00','purple':'#CC79A7','sky':'#56B4E9','yellow':'#F0E442','grey':'#7F7F7F'}
Wz,Wl=pd.read_pickle(inter('W_adj.pkl')); M=pd.read_pickle(inter('meta.pkl'))
grid=load_grid(); N=chrom_lengths(grid)
def save(fig,name): fig.savefig(f'{FO}/{name}.pdf',bbox_inches='tight'); fig.savefig(f'{FO}/{name}.png',dpi=400,bbox_inches='tight'); plt.close(fig)
def lab(ax,s,dx=-0.22,dy=1.16): ax.text(dx,dy,s,transform=ax.transAxes,fontsize=9,fontweight='bold',va='top')
bandlab=lambda cs:[str(c).replace('(','').replace(']','').replace(', ','–').replace('.0','') for c in cs]
def box(ax,data,labels,pal):
    bp=ax.boxplot(data,tick_labels=labels,patch_artist=True,showfliers=False,widths=0.62,medianprops=dict(color='white',lw=1.2),whiskerprops=dict(lw=0.6),capprops=dict(lw=0.6),boxprops=dict(lw=0.5))
    for p,s in zip(bp['boxes'],labels): p.set_facecolor(pal[s]); p.set_edgecolor('none')

# ---------------- Fig 1 (genome-wide a; b,c,d; e,f)
Lw=np.log(Wl); m=Lw.mean(axis=1); frac=(Wl>3).mean(axis=1); per=pd.Series([N[c]/k for c,k in Wl.index],index=Wl.index); u=frac>=0.9
x=np.arange(len(m)); chr_of=np.array([c for c,_ in Wl.index]); xt=[np.where(chr_of==c)[0].mean() for c in CHR]
fig=plt.figure(figsize=(W2,W2*0.9)); gs=fig.add_gridspec(3,6,height_ratios=[1.25,1,1],hspace=0.75,wspace=0.9)
ax=fig.add_subplot(gs[0,:])
for i,c in enumerate(CHR):
    sel=chr_of==c; ax.plot(x[sel],np.exp(m[sel]+0.5772),color=OI['grey'] if i%2 else '#4d4d4d',lw=0.35)
ax.axhline(1,color='k',lw=0.5,ls=':'); ax.scatter(x[u],np.exp(m[u]+0.5772),s=4,color=OI['red'],zorder=5,lw=0,label=f'Universal peaks (≥90% of 437 samples, n={u.sum()})')
for c,k,nm,dy in [('7',131,'chr7 k=131, T≈7 genes',12),('1',117,'chr1 k=117, T≈18',12),('19',565,'chr19 k=565, T≈2.6',12),('12',201,'chr12 k=201, T≈5',12),('2',168,'chr2 k=168, T≈7.5',26)]:
    i=Wl.index.get_loc((c,k)); ax.annotate(nm,(x[i],np.exp(m.iloc[i]+0.5772)),fontsize=5.5,xytext=(0,dy),textcoords='offset points',ha='center',arrowprops=dict(arrowstyle='-',lw=0.4,color='0.4'))
ax.set_xticks(xt); ax.set_xticklabels([c if (c=='X' or int(c)%2==1) else '' for c in CHR])
ax.set_yscale('log'); ax.set_ylim(0.3,60); ax.set_ylabel('Power / 1/f background'); ax.set_xlabel('Chromosome (positional frequency k increases left→right within each; odd chromosomes labelled)'); ax.legend(loc='lower left',bbox_to_anchor=(0.0,1.0)); lab(ax,'a',-0.05,1.22)
# b: collapsed on period
ax=fig.add_subplot(gs[1,0:2]); ax.hexbin(np.log10(per.values),np.log10(np.exp(m.values+0.5772)),gridsize=40,cmap='Greys',mincnt=1,linewidths=0)
ax.scatter(np.log10(per[u]),np.log10(np.exp(m[u]+0.5772)),s=4,color=OI['red'],lw=0,zorder=5)
ax.set_xticks([0.3,1,2,3]); ax.set_xticklabels(['2','10','100','1000']); ax.set_yticks([0,0.5,1]); ax.set_yticklabels(['1','3','10']); ax.set_xlabel('Period (genes)'); ax.set_ylabel('Mean power / background'); ax.set_title('All 10,007 frequencies',fontsize=6.5,pad=2); lab(ax,'b')
# c: reproducibility
ax=fig.add_subplot(gs[1,2:4]); ax.hist(frac,bins=50,color=OI['grey'],lw=0); ax.set_yscale('log'); ax.set_xlabel('Fraction of samples with power >3×'); ax.set_ylabel('Number of frequencies'); ax.axvline(0.9,color=OI['red'],ls='--',lw=0.8); ax.text(0.88,ax.get_ylim()[1]*0.6,f'n = {u.sum()}',ha='right',fontsize=6,color=OI['red']); lab(ax,'c')
# d: band fraction
ax=fig.add_subplot(gs[1,4:6]); edges=[1,2.5,3,5,10,20,50,100,300,5000]; cnt=pd.cut(per[u],edges).value_counts().sort_index(); tot=pd.cut(per,edges).value_counts().sort_index()
ax.bar(range(len(cnt)),100*cnt.values/tot.values,color=OI['red'],lw=0); ax.set_xticks(range(len(cnt))); ax.set_xticklabels(bandlab(cnt.index),rotation=45,ha='right'); ax.set_xlabel('Period (genes)'); ax.set_ylabel('% of band that is a\nuniversal peak'); lab(ax,'d')
# e: per chromosome
ax=fig.add_subplot(gs[2,0:3]); pc=pd.Series([c for c,_ in Wl.index[u]]).value_counts().reindex(CHR).fillna(0); tot_c=pd.Series([c for c,_ in Wl.index]).value_counts().reindex(CHR); ax.bar(range(len(CHR)),100*pc.values/tot_c.values,color=OI['grey'],lw=0); ax.set_xticks(range(len(CHR))); ax.set_xticklabels(CHR,fontsize=6); ax.set_xlabel('Chromosome'); ax.set_ylabel('Universal peaks\nper 100 frequencies'); lab(ax,'e',-0.14,1.18)
# f: peak by stage
ax=fig.add_subplot(gs[2,3:6]); d=Wl.loc[('7',131)]; box(ax,[d[M.index[M.estadio==s]].values for s in ST],ST,stage_pal); ax.axhline(4.6,color='0.4',ls=':',lw=0.6); ax.text(5.45,4.75,'P < 0.01',fontsize=5.5,color='0.4',ha='right'); ax.set_ylabel('Power / background\n(chr7 k=131)'); ax.set_xlabel('Fibrosis stage'); lab(ax,'f',-0.14,1.18)
save(fig,'Fig1_invariant_architecture')

# ---------------- Fig 2
P=pd.read_csv(f'{OUT}/Table_S2e_band_profile_by_condition_log.csv',index_col=0).reindex(ST); E=pd.read_csv(f'{OUT}/Table_S2d_stage_effect_per_band_log.csv').rename(columns={'band':'banda','t_stage':'t'})
SP=pd.read_csv(f'{OUT}/Table_S3a_specparam_per_sample_chromosome.csv').rename(columns={'sample':'muestra','exponent':'exponente','n_peaks':'n_picos'}); EC=pd.read_csv(f'{OUT}/Table_S3d_exponent_per_chromosome_vs_stage.csv'); EC['chr']=EC.chr.astype(str)
fig=plt.figure(figsize=(W2,W2*0.9)); gs=fig.add_gridspec(3,6,hspace=0.85,wspace=1.6,height_ratios=[1,1,0.8])
axs=[fig.add_subplot(gs[0,0:2]),fig.add_subplot(gs[0,2:4]),fig.add_subplot(gs[0,4:6]),fig.add_subplot(gs[1,0:3]),fig.add_subplot(gs[1,3:6]),fig.add_subplot(gs[2,:])]
ax=axs[0]
for s in ST: ax.plot(range(P.shape[1]),np.exp(P.loc[s].values.astype(float)+0.5772),'o-',color=stage_pal[s],ms=2.5,lw=0.9,label=s)
ax.axhline(1,color='0.5',ls=':',lw=0.6); ax.set_xticks(range(P.shape[1])); ax.set_xticklabels(bandlab(P.columns),rotation=45,ha='right'); ax.set_xlabel('Period (genes)'); ax.set_ylabel('Power / null'); ax.legend(ncol=3,handlelength=1.0,columnspacing=0.8,loc='upper left',bbox_to_anchor=(0,1.02)); ax.set_ylim(0.85,2.1); lab(ax,'a')
ax=axs[1]; ax.bar(range(len(E)),E.t,color=[OI['red'] if t>0 else OI['blue'] for t in E.t],lw=0); ax.axhline(0,color='k',lw=0.5)
for y in (2,-2): ax.axhline(y,color='0.5',ls='--',lw=0.5)
ax.set_xticks(range(len(E))); ax.set_xticklabels(bandlab(E.banda),rotation=45,ha='right'); ax.set_ylabel('t of stage\n(adj. cohort, sex)'); ax.set_xlabel('Period (genes)'); lab(ax,'b')
sm=SP.groupby('muestra').agg(offset=('offset','mean'),exponente=('exponente','mean'),n_picos=('n_picos','mean')).join(M.estadio)
for ax,var,yl,L,note in [(axs[2],'offset','Aperiodic offset','c',r'$P = 5\times10^{-15}$'),(axs[3],'n_picos','Periodic peaks per chromosome','d',r'$P = 2\times10^{-4}$'),(axs[4],'exponente','Aperiodic exponent','e',r'$P = 0.10$')]:
    box(ax,[sm[sm.estadio==s][var].values for s in ST],ST,stage_pal); ax.set_ylabel(yl); ax.set_xlabel('Fibrosis stage'); ax.set_title(note,loc='right',fontsize=6,pad=2); ax.tick_params(axis='x',labelsize=5.8); lab(ax,L)
ax=axs[5]; EC=EC.set_index('chr').reindex(CHR).reset_index(); c_=[OI['red'] if (q<0.05 and t>0) else OI['blue'] if (q<0.05 and t<0) else '0.75' for t,q in zip(EC.t,EC.q)]
ax.bar(range(len(EC)),EC.t,color=c_,lw=0); ax.set_xticks(range(len(EC))); ax.set_xticklabels(CHR,fontsize=5.5,rotation=90); ax.tick_params(axis='x',length=1.5); ax.axhline(0,color='k',lw=0.5); ax.set_ylabel('t of stage on exponent'); ax.set_xlabel('Chromosome'); ax.tick_params(axis='x',labelsize=6.5,rotation=0); lab(ax,'f',-0.07)
save(fig,'Fig2_stage_effects_spectrum')

# ---------------- Fig 3
Pp=pd.read_pickle(inter('pairs.pkl')).rename(columns={'intergenic_bp':'dist_bp','same_tad':'mismo_tad'}); Pp['lag']=1
V=pd.read_csv(f'{OUT}/Table_S5c_contiguous_DE_neighbourhoods.csv'); Kt=pd.read_pickle(inter('K_tad.pkl'))
V['mismo_TAD']=[ (lambda g:(g.tad.nunique()==1) if len(g)>=3 else np.nan)(Kt[(Kt.chr==str(r.chr))&(Kt.grid_index>=r.start_grid)&(Kt.grid_index<=r.end_grid)&(Kt.tad>=0)]) for _,r in V.iterrows()]
R=pd.read_csv(f'{OUT}/Table_S5d_acf_with_without_composition.csv').rename(columns={'acf_unadjusted':'sin_ajustar','acf_composition_adjusted':'ajustada_composicion'})
Pq=pd.read_pickle(inter('pairs.pkl')).rename(columns={'paralog_family':'familia','intergenic_bp':'inter_bp'}); Pq['dbin']=pd.cut(Pq.inter_bp,[-1e9,0,1e3,1e4,5e4,1e5,5e5,1e9],labels=['solapan','<1 kb','1–10 kb','10–50 kb','50–100 kb','100–500 kb','>500 kb'])
fig,axs=plt.subplots(2,2,figsize=(W2*0.8,W2*0.72)); axs=axs.ravel(); plt.subplots_adjust(hspace=0.65,wspace=0.45)
ax=axs[0]; ax.axhspan(-0.017,0.017,color='0.9',lw=0,label='95% permutation null'); ax.plot(R.lag,R.sin_ajustar,'o-',color='k',ms=3,lw=0.9,label='Stage effect per gene'); ax.plot(R.lag,R.ajustada_composicion,'s--',color=OI['orange'],ms=3,lw=0.9,label='Adjusted for cell composition'); ax.set_xscale('log'); ax.set_xlabel('Distance (genes)'); ax.set_ylabel('Spatial autocorrelation'); ax.legend(); lab(ax,'a')
ax=axs[1]; q=Pq[~Pq.familia]; od=['solapan','<1 kb','1–10 kb','10–50 kb','50–100 kb','100–500 kb','>500 kb']; t=q.groupby('dbin').apply(lambda h:pd.Series({'r':np.corrcoef(h.t1,h.t2)[0,1],'rc':np.corrcoef(h.t1c,h.t2c)[0,1],'co':h.coexpr.mean()})).reindex(od)
xx=list(range(len(od))); ax.plot(xx,t.r,'o-',color='k',ms=3,lw=0.9,label='Concordance of stage effect'); ax.plot(xx,t.rc,'o--',color=OI['orange'],ms=3,lw=0.9,label='…adjusted for composition'); ax.plot(xx,t.co,'s-',color=OI['blue'],ms=3,lw=0.9,label='Co-expression across samples')
ax.set_xticks(xx); ax.set_xticklabels(['overlap','<1','1–10','10–50','50–100','100–500','>500'],rotation=45,ha='right'); ax.set_xlabel('Intergenic distance (kb)\nnon-paralog neighbours'); ax.set_ylabel('Correlation'); ax.axhline(0,color='0.6',lw=0.5); ax.legend(loc='upper right',bbox_to_anchor=(1.0,1.02)); ax.set_ylim(-0.08,0.5); lab(ax,'b')
ax=axs[2]; q2=Pp[Pp.lag==1].dropna(subset=['stat1','stat2']).copy(); bins=[0,2e4,5e4,1e5,2e5,5e5,1e6,5e6]; q2['db']=pd.cut(q2.dist_bp,bins)
for same,col,nm in [(True,OI['red'],'Same TAD'),(False,OI['blue'],'Different TADs')]:
    g=q2[q2.mismo_tad==same].groupby('db',observed=True); r=g.apply(lambda h:np.corrcoef(h.stat1,h.stat2)[0,1] if len(h)>=60 else np.nan); xs=np.array([b.right/1e3 for b in r.index]); ax.plot(xs,r.values,'o-',color=col,ms=3,lw=0.9,label=nm)
ax.set_xscale('log'); ax.set_xlabel('Distance between neighbours (kb)'); ax.set_ylabel('Concordance of stage effect'); ax.axhline(0,color='0.6',lw=0.5); ax.set_ylim(-0.2,0.45); ax.legend(title='Liver TADs (Leung et al. 2015)',loc='upper right',bbox_to_anchor=(1.0,0.9)); ax.text(0.98,0.98,'Δ at equal distance: P ≈ 0.5',transform=ax.transAxes,ha='right',va='top',fontsize=6); lab(ax,'c')
ax=axs[3]; obs=[V[V.n_genes>=n].mismo_TAD.astype(float).mean()*100 for n in (3,4,5)]; ax.axhspan(77.4-7.8,77.4+7.8,color='0.9',lw=0,zorder=0); ax.axhline(77.4,color='k',ls='--',lw=0.7,label='Chance (77 ± 4%)',zorder=1); ax.bar([0,1,2],obs,color=OI['grey'],width=.6,lw=0,zorder=2)
for i,o in enumerate(obs): ax.text(i,o-3,f'{o:.0f}%',ha='center',va='top',fontsize=6,color='white',fontweight='bold')
ax.set_xticks([0,1,2]); ax.set_xticklabels(['≥3','≥4','≥5']); ax.set_xlabel('Genes per DE neighbourhood'); ax.set_ylabel('% fully within one TAD'); ax.set_ylim(50,100); ax.legend(loc='upper right'); lab(ax,'d')
save(fig,'Fig3_spatial_coupling_TADs')

# ---------------- Fig 4
C=pd.read_csv(f'{OUT}/Table_S6b_composition_scores_per_sample.csv',index_col=0).rename(columns={'n_peaks':'n_picos'}); B=pd.read_csv(f'{OUT}/Table_S7a_threshold_vs_linear.csv'); B.columns=[c.replace('dAIC_step_at_','ΔAIC_salto_en_') for c in B.columns]; B['variable']=B.variable.replace({'n_peaks':'n_picos'})
fig=plt.figure(figsize=(W2*0.85,W2*0.72)); gs=fig.add_gridspec(2,2,hspace=0.55,wspace=0.45,height_ratios=[1,0.95])
axs=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,:])]
ax=axs[0]; ct_cols={'Hepatocito':('Hepatocyte',OI['green']),'HSC':('Stellate / myofibroblast',OI['red']),'Colangiocito':('Cholangiocyte',OI['purple']),'Macrofago':('Macrophage',OI['orange']),'Linfocito':('Lymphocyte',OI['blue']),'Endotelio':('Endothelium',OI['grey'])}
for ct,(nm,col) in ct_cols.items(): g=C.groupby('estadio')[ct].agg(['mean','sem']).reindex(ST); ax.errorbar(range(6),g['mean'],yerr=g['sem'],fmt='o-',color=col,ms=2.5,lw=0.9,capsize=1.5,elinewidth=0.6,label=nm)
ax.set_xticks(range(6)); ax.set_xticklabels(ST); ax.axhline(0,color='0.7',lw=0.5); ax.set_ylabel('Composition score (marker z)'); ax.set_xlabel('Fibrosis stage'); ax.set_ylim(-0.9,1.6); ax.legend(ncol=2,handlelength=1.2,loc='upper left',columnspacing=0.8,fontsize=5.5); lab(ax,'a')
ax=axs[1]
for st in ST: g=C[C.estadio==st]; ax.scatter(g.Hepatocito,g.offset,s=6,color=stage_pal[st],label=st,alpha=.85,lw=0)
ax.set_xlabel('Hepatocyte score'); ax.set_ylabel('Aperiodic offset (spectral variance)'); ax.set_ylim(3.09,3.41); ax.text(0.03,0.98,'Stage effect on offset:\nt = −8.0 → +1.4 after\nadjusting for composition',transform=ax.transAxes,fontsize=6,va='top'); ax.set_xlim(-2.3,2.1); ax.legend(loc='lower right',ncol=1,markerscale=1.5,title='Stage',title_fontsize=6); lab(ax,'b')
ax=axs[2]; H=B.set_index('variable')[[c for c in B.columns if c.startswith('ΔAIC')]]; H.columns=['≥F1','≥F2','≥F3','≥F4']; H=H.loc[['Hepatocito','HSC','Colangiocito','Macrofago','Linfocito','Endotelio','offset','n_picos']]; H.index=['Hepatocyte','Stellate','Cholangiocyte','Macrophage','Lymphocyte','Endothelium','Spectral offset','Peak count']
H=H.T; im=ax.imshow(H.values,cmap='RdBu_r',vmin=-15,vmax=15,aspect='auto'); ax.set_yticks(range(4)); ax.set_yticklabels(H.index); ax.set_xticks(range(H.shape[1])); ax.set_xticklabels(H.columns,rotation=30,ha='right'); ax.set_ylabel('Jump from stage')
for i in range(H.shape[0]):
    for j_ in range(H.shape[1]): ax.text(j_,i,f'{H.values[i,j_]:.0f}',ha='center',va='center',fontsize=5.5,color='white' if abs(H.values[i,j_])>9 else 'k')
cb=plt.colorbar(im,ax=ax,fraction=0.025,pad=0.02); cb.set_label('ΔAIC (jump − linear)',fontsize=6); cb.ax.tick_params(labelsize=5.5); ax.spines[['left','bottom']].set_visible(False); ax.tick_params(length=0); lab(ax,'c',-0.1,1.25)
save(fig,'Fig4_composition_and_thresholds')

# ---------------- Fig 5
R=pd.read_csv(f'{OUT}/Table_S8c_within_type_spatial_autocorrelation.csv').rename(columns={'type':'tipo','acf_lag1':'acf_lag1_DE','null_p97.5':'nulo_p97.5'}); CO=pd.read_csv(f'{OUT}/Table_S8b_within_type_neighbour_coexpression.csv').rename(columns={'type':'tipo','r_neighbours':'r_vecinos','r_random':'r_aleatorio'}); Cn=pd.read_csv(f'{OUT}/Table_S9a_candidate_targets_lineage_intrinsic.csv').rename(columns={'t_bulk_comp_adjusted':'t_bulk_comp_adjusted'})
Cn['development_stage']=Cn.get('development_stage','not curated')
en={'Fagocito_mononuclear':'Mono. phagocyte','Endotelio':'Endothelium','Mesenquima_HSC':'Mesenchyme/HSC','Colangiocito':'Cholangiocyte','Plasma':'Plasma','T_NK':'T/NK','B':'B','pDC':'pDC','Hepatocito':'Hepatocyte'}
fig=plt.figure(figsize=(W2*0.85,W2*0.95)); gs=fig.add_gridspec(2,2,hspace=0.55,wspace=0.7,height_ratios=[0.85,1])
axs=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,:])]
ax=axs[0]; R=R.sort_values('acf_lag1_DE'); ax.barh(range(len(R)),R.acf_lag1_DE,color=OI['red'],lw=0,height=0.65); ax.errorbar([0]*len(R),range(len(R)),xerr=R['nulo_p97.5'],fmt='none',color='k',capsize=2,lw=0.6,label='95% null'); ax.set_yticks(range(len(R))); ax.set_yticklabels([en[t] for t in R.tipo]); ax.set_xlabel('Spatial autocorrelation of\ncirrhosis-vs-healthy change'); ax.legend(loc='lower right'); lab(ax,'a',-0.55,1.12)
ax=axs[1]; co=CO.groupby('tipo')[['r_vecinos','r_aleatorio']].mean().sort_values('r_vecinos'); y=np.arange(len(co)); ax.barh(y+0.18,co.r_vecinos,0.36,color=OI['blue'],lw=0,label='Neighbouring genes'); ax.barh(y-0.18,co.r_aleatorio,0.36,color='0.65',lw=0,label='Random pairs'); ax.set_yticks(y); ax.set_yticklabels([en[t] for t in co.index]); ax.set_xlabel('Co-expression across single cells'); ax.legend(loc='lower right'); lab(ax,'b',-0.55,1.12)
ax=axs[2]; top=Cn.head(28); stg=top.development_stage.fillna('not curated')
colmap={'clinical, liver':OI['red'],'clinical, inflammasome':OI['orange'],'clinical, other fibrosis':OI['orange'],'clinical, other':OI['orange'],'approved, other':OI['orange'],'preclinical':OI['sky'],'tool compounds':OI['sky'],'biomarker':OI['green'],'none':'0.7','not curated':'0.85'}
ax.scatter(top.t_bulk_comp_adjusted,top.t_sc_max,s=22,c=[colmap.get(s,'0.8') for s in stg],lw=0.3,edgecolor='k',zorder=3)
texts=[ax.text(r.t_bulk_comp_adjusted,r.t_sc_max,r.gene_name,fontsize=5) for _,r in top.iterrows()]; ax.set_xlim(0.8,8.5); ax.set_ylim(1.7,5.6); adjust_text(texts,ax=ax,arrowprops=dict(arrowstyle='-',color='0.5',lw=0.3),expand=(1.4,1.6),force_text=(0.5,0.8))
hd=[Line2D([],[],marker='o',ls='',color=c,markeredgecolor='k',markeredgewidth=0.3,ms=4,label=l) for l,c in [('Clinical, liver',OI['red']),('Clinical, other indication',OI['orange']),('Preclinical / tool compound',OI['sky']),('Biomarker only',OI['green']),('No modulator','0.7')]]
ax.legend(handles=hd,loc='upper right',ncol=2); ax.set_xlabel('Bulk t (stage, composition-adjusted)'); ax.set_ylabel('Max t within a cell type'); lab(ax,'c',-0.1,1.12)
save(fig,'Fig5_single_cell_and_targets'); print('ok')
