from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy.stats import rankdata
from openpyxl import load_workbook
R=Path(__file__).resolve().parent;O=R/'revision_figures';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':240})
def save(n,fig):fig.savefig(O/f'Figure_{n}.png',bbox_inches='tight',facecolor='white');plt.close(fig)
w=load_workbook(R/'data/raw/historical_benchmark_report.xlsx',data_only=True)
tab=[]
for s in list(w)[:2]:
 vals=list(s.values);hdr=vals[25];tab.append([{hdr[i]:r[i] for i in range(11)} for r in vals[26:] if r[1] is not None])
(R/'revision_analysis/historical_benchmark_metrics.json').write_text(json.dumps(tab,indent=2))
colors=['#276A9A','#A65628','#6A5ACD','#64825A','#D57C23','#9B6F9D','#808080'];cs={r['Model']:c for r,c in zip(tab[0],colors)}
fig,ax=plt.subplots(1,2,figsize=(9,3.1));ions=['Ca','none','Ba','Sr','Cu','Mg'];counts=[98,18,13,11,8,4]
ax[0].bar(ions,counts,color=['#276A9A','#888888']+['#BD6B20']*4);ax[0].set(title='(a) Historical subset B',ylabel='Records');ax[0].bar_label(ax[0].containers[0]);ax[0].set_ylim(0,115)
ax[1].bar(['F1','F2','F3','F4','F5'],[45,27,27,27,26],color='#699DC1');ax[1].bar_label(ax[1].containers[0]);ax[1].set(title='(b) Reported grouped outer folds',ylabel='Records',ylim=(0,53));fig.tight_layout();save(1,fig)
fig,ax=plt.subplots(figsize=(9,5.6));ax.set(xlim=(0,10),ylim=(0,6));ax.axis('off')
def box(x,y,w,h,t,fc='#EEF3F6'):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.06',fc=fc,ec='#59778A',lw=1));ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=10)
def arrow(a,b):ax.annotate('',xy=b,xytext=a,arrowprops={'arrowstyle':'->','color':'#52636D','lw':1.2})
box(.2,4.75,4.2,1,'Historical benchmarks\nA: 116 records, seven model families\nB: 152 records, three model families');box(5.25,4.75,4.3,1,'New series and historical data\n34 positive new records, 16 groups\n164 positive historical records')
box(.2,2.9,4.2,1.15,'Grouped historical validation\nDifferent original splits and GP implementations\nRanks are descriptive, not a selection rule');arrow((2.3,4.75),(2.3,4.05))
box(5.25,2.8,4.3,1.3,'Outer composition held out evaluation\nInner grouped tuning, training only scaling\nHistorical composition matches removed\nHeadline: 30 records, 15 groups');arrow((7.4,4.75),(7.4,4.1))
box(5.25,1.05,4.3,1.05,'Separate final GP refit\nAll 198 positive records\nExploratory choice, not the lowest RMSE model','#F9F0DC');arrow((7.4,2.8),(7.4,2.1))
box(.2,.2,4.2,1.7,'Computational candidate screening\nStock and dose search + support flags\nConditional Monte Carlo scenarios\nRecipes require prospective preparation','#EDF4EC');arrow((5.25,1.55),(4.4,1.2));fig.tight_layout();save(2,fig)
fig,axs=plt.subplots(1,2,figsize=(9,3.8),gridspec_kw={'width_ratios':[1.15,1]})
for a,rows,t in zip(axs,tab,['(a) Ca / no ion','(b) Five divalent ions / no ion']):
 names=[r['Model'].replace('Gaussian Process','GP').replace('Polynomial Ridge','Polynomial ridge').replace('Linear regression','Linear') for r in rows];v=[r['RMSE_log10'] for r in rows];a.barh(names,v,color=[cs[r['Model']] for r in rows]);a.invert_yaxis();a.set(xlabel='Out of fold RMSE of log₁₀(G/Pa)',title=t,xlim=(0,max(v)+.24))
 for i,r in enumerate(rows):a.text(r['RMSE_log10']+.008,i,f"{r['RMSE_log10']:.3f}  (R² {r['R2_log10']:.3f})",va='center',fontsize=8)
fig.tight_layout();save(3,fig)
keys=['RMSE_Pa','RMSE_log10','R2_log10','Median_factor','Geomean_factor','Within_factor2'];labels=['RMSE\nPa','RMSE\nlog','R²\nlog','Median\nfactor','GMF','Within\n×2']
fig,axs=plt.subplots(1,2,figsize=(9.2,3.8),gridspec_kw={'width_ratios':[1.25,1]})
for a,rows,t in zip(axs,tab,['(a) Subset A','(b) Subset B']):
 ranks=np.stack([rankdata([(-1 if k in ['R2_log10','Within_factor2'] else 1)*r[k] for r in rows],method='average') for k in keys],axis=1);a.imshow(ranks,cmap='YlGnBu_r',vmin=1,vmax=7,aspect='auto');a.set_xticks(range(6),labels=labels,fontsize=8);a.xaxis.tick_top();a.set_yticks(range(len(rows)),labels=[r['Model'].replace('Gaussian Process','GP').replace('Polynomial Ridge','Poly ridge') for r in rows],fontsize=9);a.set_title(t,pad=42)
 for i in range(len(rows)):
  for j in range(6):a.text(j,i,f'{ranks[i,j]:g}',ha='center',va='center',color='white' if ranks[i,j]<3 else '#222222',fontsize=10)
fig.tight_layout();save(4,fig)
fig,ax=plt.subplots(figsize=(8,2.8))
for i,n in enumerate(['Gaussian Process','RBF-SVR','XGBoost']):
 a=next(r for r in tab[0] if r['Model']==n)['RMSE_log10'];b=next(r for r in tab[1] if r['Model']==n)['RMSE_log10'];ax.plot([a,b],[i,i],color='#B3B9BD',lw=3);ax.scatter(a,i,c=cs[n],s=60);ax.scatter(b,i,facecolor='white',edgecolor=cs[n],s=65);ax.text(a-.01,i+.17,f'{a:.3f}',ha='right',fontsize=9);ax.text(b+.01,i,f'{b:.3f}  (×{b/a:.2f})',va='center',fontsize=9)
ax.set_yticks(range(3),['GP','RBF-SVR','XGBoost']);ax.set_ylim(2.5,-.35);ax.set(xlabel='Out of fold RMSE of log₁₀(G/Pa)',xlim=(.08,.75));fig.tight_layout();save(5,fig)
diag=json.load(open(R/'analysis_v5/diagnostics.json'));rows=diag['rows'];ioncol={'Ca':'#245C8D','Ba':'#C27522','Sr':'#AD9900','Cu':'#32875E','Al':'#8455A9','Fe':'#98634C'}
fig=plt.figure(figsize=(8.4,6));gs=fig.add_gridspec(2,2,height_ratios=[1.5,1]);ax=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,:])]
z=np.geomspace(2000,100000,100);ax[0].fill_between(z,z/2,z*2,color='#E5EEF5');ax[0].plot(z,z,c='#444444',lw=1)
for ion in ['Ca','Ba','Sr','Cu','Al','Fe']:
 rr=[r for r in rows if r['ion']==ion];ax[0].scatter([r['observed'] for r in rr],[r['prediction'] for r in rr],s=20,color=ioncol[ion],label=f'{ion} (n={len(rr)})')
ax[0].set(xscale='log',yscale='log',xlim=(1800,110000),ylim=(1800,110000),xlabel='Observed fitted G (Pa)',ylabel='Held out GP prediction (Pa)',title='(a) New series validation');ax[0].legend(fontsize=6,loc='lower right');ax[0].text(.03,.96,'n = 30; RMSE 0.303\n22/30 within ×2',transform=ax[0].transAxes,va='top',fontsize=8)
families=['ALG+SCLG','ALG+TOCNF','ALG+TOCNF+SCLG','TOCNF']
for i,f in enumerate(families):
 rr=[r for r in rows if r['family']==f]
 for j,r in enumerate(rr):ax[1].scatter(i+(j%3-1)*.06,r['factor'],s=20,c=ioncol[r['ion']])
 ax[1].plot([i-.19,i+.19],[10**np.median([r['abs_log_error'] for r in rr])]*2,c='black',lw=2)
ax[1].axhline(2,c='#AAAAAA',ls='--');ax[1].set(yscale='log',ylabel='Fold error',title='(b) Polymer system',ylim=(.95,5.3));ax[1].set_xticks(range(4),['ALG +\nSCLG\nn=9','ALG +\nTOCNF\nn=9','Ternary\nn=11','TOCNF\nn=1'],fontsize=8);ax[1].set_yticks([1,1.5,2,3,5],labels=['1','1.5','2','3','5'])
for i,r in enumerate(diag['replicates']):ax[2].scatter(r['ratio_KS_old'],i,c=ioncol[r['ion']],s=30);ax[2].text(r['ratio_KS_old']+.035,i,f"×{r['factor']:.2f}",va='center',fontsize=8)
ax[2].axvline(1,c='#AAAAAA');ax[2].set_yticks(range(4),['Ba near repeat','Cu repeat','Fe repeat','Al near repeat'],fontsize=9);ax[2].invert_yaxis();ax[2].set(xlim=(.3,1.35),xlabel='New / earlier G',title='(c) Control pairs');fig.tight_layout();save(6,fig)
mc=json.load(open(R/'revision_analysis/mc_sensitivity.json'));mc15=[r for r in mc if r['target_kPa']==15];G=np.load(R/'revision_analysis/candidate15_MC.npz')['G']/1000
fig,ax=plt.subplots(1,2,figsize=(9,3.5));ax[0].hist(G,bins=32,color='#358D85');ax[0].axvspan(13.5,16.5,color='#E6B46B',alpha=.3,label='15 kPa ±10%');ax[0].axvline(13.6500197,color='#333333',ls='--',label='Nominal prediction');ax[0].legend(fontsize=8);ax[0].set(xlabel='Predicted G (kPa)',ylabel='Monte Carlo draws',title='(a) Conditional preparation spread')
scens=['polymer_only','volumes_only','crosslinker_only','nominal'];vs=[next(r['relative_halfwidth_percent'] for r in mc15 if r['scenario']==s) for s in scens];ax[1].bar(range(4),vs,color='#358D85');ax[1].set_xticks(range(4),['Polymer\nstocks','Dispensed\nvolumes','Crosslinker\nstock','All inputs']);ax[1].bar_label(ax[1].containers[0],fmt='%.2f');ax[1].set(ylim=(0,10),ylabel='Relative 5–95% half width (%)',title='(b) Conditional response half-widths');fig.tight_layout();save(7,fig)
mm=json.load(open(R/'analysis_source/mcmc_validation.json'));rr=[r for r in mm['rows'] if not r['legacy_overlap']];arr=np.load(R/'analysis_source/mcmc/full.npz')['draws'].reshape(-1,3)
fig=plt.figure(figsize=(8.4,6.2));gs=fig.add_gridspec(2,2,height_ratios=[1.1,1]);ax=[fig.add_subplot(gs[0,:]),fig.add_subplot(gs[1,0]),fig.add_subplot(gs[1,1])];z=np.geomspace(1800,150000,100);ax[0].plot(z,z,c='#444444');ax[0].fill_between(z,z/2,z*2,color='#E5EEF5')
for r in rr:
 lo,hi=r['predictive80_Pa'];c='#35875A' if lo<=r['observed']<=hi else '#BE5050';ax[0].plot([r['observed']]*2,[lo,hi],c=c,lw=.8);ax[0].scatter(r['observed'],r['G_Pa'],c=c,s=15)
ax[0].set(xscale='log',yscale='log',xlim=(1800,150000),ylim=(500,170000),xlabel='Observed fitted G (Pa)',ylabel='MCMC prediction (Pa)',title='(a) Nominal 80% intervals');ax[0].text(.03,.97,'Coverage 18/30 (60%)',transform=ax[0].transAxes,va='top',fontsize=8)
ax[1].hist(np.exp(arr[:,0]),bins=35,density=True,color='#699CC1');ax[1].axvline(4,c='#CF7C22');ax[1].set(xlabel='Length scale',ylabel='Posterior density',title='(b) Full data posterior')
ax[2].hist(arr[:,1]/np.log(10),bins=35,density=True,color='#699CC1');ax[2].axvline(np.log10(.3),c='#CF7C22');ax[2].set(xlabel='log₁₀ noise variance',title='(c) Full data posterior');fig.tight_layout();save(8,fig)
