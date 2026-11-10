"""Revision audit. Reuses archived inputs; never refits or fabricates rheology.
Run from the parent of analysis_source, analysis_v4 and analysis_v5.
"""
from pathlib import Path
import json, csv, sys, hashlib, platform
import numpy as np
from scipy.spatial.distance import cdist
from scipy.special import expit
from scipy.stats import norm
from scipy.optimize import brentq
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.ensemble import HistGradientBoostingRegressor
R=Path(__file__).resolve().parent;O=R/'revision_analysis';O.mkdir(exist_ok=True)
sys.path.insert(0,str(R/'analysis_v5'));import train_v5 as v5
sys.path.insert(0,str(R/'analysis_source'));import mcmc_gp as mg
from monte_carlo import simulate,draws
sys.path.insert(0,str(R/'analysis_source/model'));from hydrogel import HydrogelModel
def save(name,x): (O/name).write_text(json.dumps(x,indent=2,default=lambda a:a.item() if hasattr(a,'item') else str(a)))
def csvout(name,rows):
 with (O/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
rs=json.loads((R/'analysis_source/records.json').read_text());ks=v5.ks;old=v5.old;g=v5.g;y=v5.y
ledger=[]
for r in rs:
 ledger.append(dict(record_id=r['id'],source=r['source'],group='|'.join(f'{x:.8f}' for x in r['c']),family=r['family'],ion=r['ion'],c_ALG_wt_pct=r['c'][0],c_TOCNF_wt_pct=r['c'][1],c_SCLG_wt_pct=r['c'][2],ion_final_M=r['dose_M'],ion_stock_M=2*r['dose_M'],G_Pa=r['G_observed_Pa'],gel_formed=int(r['gel_formed']),historical_A=int(r['source']=='legacy' and r['gel_formed'] and r['ion'] in ['Ca','none']),historical_B=int(r['source']=='legacy' and r['gel_formed'] and r['ion'] not in ['Fe','Al']),new_validation=int(r['source']=='KS' and r['gel_formed'] and not r.get('legacy_matches')),repeat_control=int(bool(r.get('legacy_matches'))),legacy_match_ids=';'.join(map(str,r.get('legacy_matches',[]))),preparation_id='not recorded',batch_id='not recorded'))
csvout('record_ledger.csv',ledger)
counts={s:{'records':sum(r['source']==s for r in rs),'positive':sum(r['source']==s and bool(r['gel_formed']) for r in rs),'negative_ids':[r['id'] for r in rs if r['source']==s and not r['gel_formed']],'positive_groups':len({v5.group(r) for r in rs if r['source']==s and r['gel_formed']})} for s in ['legacy','KS']};save('counts.json',counts)
fold_manifest=[];reproduced=[];svr_pred=np.empty(len(ks));hgb_pred=np.empty(len(ks));gp_pred=np.empty(len(ks));mean_pred=np.empty(len(ks))
svrgrid=[(c,ga,e) for c in [1,10,100] for ga in [.1,1.] for e in [.05,.15]]
hgbgrid=[(lr,le,reg) for lr in [.05,.1] for le in [3,7,15] for reg in [1.,10.]]
def features(rows):
 return np.array([[*r['c'],np.log10(1+r['dose_M']/.005),*[float(r['ion']==i) for i in v5.ions],float(r['source']=='KS')] for r in rows])
def altfit(kind,p,train,test):
 sc=StandardScaler().fit(features(train));a=sc.transform(features(train));b=sc.transform(features(test));yy=np.log10([r['G_observed_Pa'] for r in train])
 if kind=='svr': m=SVR(C=p[0],gamma=p[1],epsilon=p[2])
 else:m=HistGradientBoostingRegressor(learning_rate=p[0],max_leaf_nodes=p[1],l2_regularization=p[2],max_iter=150,min_samples_leaf=5,early_stopping=False,random_state=20261101)
 return m.fit(a,yy).predict(b)
def alttune(kind,grid,tr):
 best=None
 for p in grid:
  err=[]
  for aa,bb in GroupKFold(4).split(tr,groups=g[tr]):
   a,b=tr[aa],tr[bb];pred=altfit(kind,p,v5.trainset(a),[ks[i] for i in b]);err.extend((pred-y[b])**2)
  mse=float(np.mean(err))
  if best is None or mse<best[0]:best=(mse,p)
 return best[1]
original=json.loads((R/'analysis_v5/comparison.json').read_text());orig={r['id']:r for r in original['rows']}
for no,(tr,te) in enumerate(LeaveOneGroupOut().split(y,groups=g),1):
 train=v5.trainset(tr);test=[ks[i] for i in te];held={v5.group(r) for r in test}
 assert not held.intersection({v5.group(r) for r in train})
 inner=[]
 for aa,bb in GroupKFold(4).split(tr,groups=g[tr]):
  a,b=tr[aa],tr[bb];traininner=v5.trainset(a);excluded=held|{v5.group(ks[i]) for i in b};assert not excluded.intersection({v5.group(r) for r in traininner})
  inner.append({'training_ids':[r['id'] for r in traininner],'validation_ids':[ks[i]['id'] for i in b]})
 _,gp_param,_=v5.tune('v4_gp',tr);gp_pred[te]=v5.predict(v5.fit('v4_gp',gp_param,tr),test);mean_pred[te]=y[tr].mean()
 sp=alttune('svr',svrgrid,tr);hp=alttune('hgb',hgbgrid,tr)
 svr_pred[te]=altfit('svr',sp,train,test);hgb_pred[te]=altfit('hgb',hp,train,test)
 fold_manifest.append(dict(fold=no,held_ids=[r['id'] for r in test],training_ids=[r['id'] for r in train],inner=inner,gp_param=gp_param,svr_param=sp,hgb_param=hp))
 reproduced.extend(abs(10**gp_pred[i]/orig[ks[i]['id']]['predictions']['v4_gp']-1) for i in te)
 v5.cache.clear();print('Nested audit fold',no,'of 16',flush=True)
save('fold_manifest.json',fold_manifest)
nov=np.array([not bool(r.get('legacy_matches')) for r in ks]);gn=g[nov];yn=y[nov]
preds={k:np.log10([orig[r['id']]['predictions'][k] for r in ks])[nov] for k in original['metrics_novel']};preds.update({'training_mean':mean_pred[nov],'RBF_SVR_added':svr_pred[nov],'hist_gradient_boosting_added':hgb_pred[nov]})
metrics={}
unique=sorted(set(gn))
for k,p in preds.items():
 m=v5.metrics(yn,p);m['group_weighted_RMSE_log10']=float(np.sqrt(np.mean([np.mean((yn[gn==gg]-p[gn==gg])**2) for gg in unique])));m['group_weighted_MAE_log10']=float(np.mean([np.mean(abs(yn[gn==gg]-p[gn==gg])) for gg in unique]));metrics[k]=m
save('new_series_metrics.json',metrics)
rows=[]
for j,r in enumerate(np.array(ks,dtype=object)[nov]):rows.append(dict(id=r['id'],group=gn[j],family=r['family'],ion=r['ion'],observed_Pa=r['G_observed_Pa'],**{k:float(10**v[j]) for k,v in preds.items()}))
csvout('new_series_predictions.csv',rows)
rng=np.random.default_rng(20261101);boots=rng.integers(0,len(unique),(20000,len(unique)));ci={}
for k,p in preds.items():
 a=np.array([np.mean((yn[gn==gg]-p[gn==gg])**2) for gg in unique]);b=np.array([np.mean((yn[gn==gg]-preds['v4_gp'][gn==gg])**2) for gg in unique]);delta=a-b
 ci[k]={'group_mean_MSE_difference_vs_GP':float(delta.mean()),'conditional_cluster_bootstrap95':np.quantile(delta[boots].mean(axis=1),[.025,.975]).tolist()}
save('paired_group_bootstrap.json',ci)
# Retrospective historical-only test. Selection uses historical labels only.
histgroups=np.array(['|'.join(f'{v:.8f}' for v in r['c']) for r in old]);hy=np.log10([r['G_observed_Pa'] for r in old]);grid=[('raw',l,n,.9) for l in [.5,1.5,4] for n in [.03,.3]];scores=[]
for p in grid:
 er=[]
 for tr,te in GroupKFold(4).split(old,groups=histgroups):
  m=v5.gp_fit([old[i] for i in tr],*p);pr=v5.predict(m,[old[i] for i in te]);er.extend((hy[te]-pr)**2)
 scores.append(float(np.mean(er)))
par=grid[int(np.argmin(scores))];m=v5.gp_fit(old,*par);test=[dict(r,source='legacy') for r in ks];frozen=v5.predict(m,test)
save('historical_only.json',{'description':'Retrospective test, not a prospective experiment. No later-series target used in fitting/tuning; same source encoding because between-source correlation is not identifiable from one source.','param':par,'historical_inner_scores':scores,'metrics_novel':v5.metrics(y[nov],frozen[nov]),'predictions':[dict(id=r['id'],predicted_Pa=float(10**p)) for r,p in zip(ks,frozen)]})
# Verify diagnostics from archived chains and calculate held-out calibration.
diagnostics=[];calrows=[];originalmc=json.loads((R/'analysis_source/mcmc_validation.json').read_text());om={r['id']:r for r in originalmc['rows']};mc_errors=[]
for f in ['full']+[f'fold_{i:02}' for i in range(1,17)]:
 arr=np.load(R/'analysis_source/mcmc'/f'{f}.npz')['draws'];di=mg.diagnostics(arr);diagnostics.append(dict(name=f,shape=list(arr.shape),parameters=di))
 if f=='full':continue
 no=int(f.split('_')[1]);held=sorted(set(mg.grp(r) for r in ks))[no-1];tr=[r for r in old+ks if mg.grp(r)!=held];te=[r for r in ks if mg.grp(r)==held];gp=mg.GP(tr);mus=[];vars=[]
 for z in mg.selected(arr):
  mu,_,vv=gp.conditional(z,te);mus.append(mu);vars.append(vv)
 mus=np.array(mus);sds=np.sqrt(vars)
 for j,r in enumerate(te):
  point=10**mus[:,j].mean();mc_errors.append(abs(point/om[r['id']]['G_Pa']-1))
  if r.get('legacy_matches'):continue
  true=np.log10(r['G_observed_Pa']);out=dict(id=r['id'],group=held,observed_log10=float(true),predicted_log10=float(mus[:,j].mean()))
  for level in [.5,.8,.9,.95]:
   alpha=1-level
   def q(p):return brentq(lambda v:np.mean(norm.cdf((v-mus[:,j])/sds[:,j]))-p,float(np.min(mus[:,j]-12*sds[:,j])),float(np.max(mus[:,j]+12*sds[:,j])))
   low,high=q(alpha/2),q(1-alpha/2);score=high-low+2/alpha*max(low-true,0)+2/alpha*max(true-high,0);tag=str(int(level*100));out.update({f'lo{tag}':low,f'hi{tag}':high,f'covered{tag}':int(low<=true<=high),f'interval_score{tag}':score})
  calrows.append(out)
save('mcmc_diagnostics_recomputed.json',diagnostics);csvout('mcmc_calibration_rows.csv',calrows)
cal=[]
for level in [50,80,90,95]:
 tag=str(level);a=np.array([np.mean([r['covered'+tag] for r in calrows if r['group']==gg]) for gg in unique]);cal.append(dict(nominal_percent=level,covered_records=sum(r['covered'+tag] for r in calrows),N=len(calrows),record_coverage=float(np.mean([r['covered'+tag] for r in calrows])),group_weighted_coverage=float(a.mean()),group_coverage_bootstrap95=np.quantile(a[boots].mean(axis=1),[.025,.975]).tolist(),mean_interval_score_log10=float(np.mean([r['interval_score'+tag] for r in calrows]))))
save('mcmc_calibration.json',cal)
summary={'max_relative_GP_prediction_difference':max(reproduced),'max_relative_MCMC_prediction_difference':max(mc_errors),'max_Rhat':max(x['rhat'] for d in diagnostics for x in d['parameters'].values()),'min_ESS':min(x['ess_bulk'] for d in diagnostics for x in d['parameters'].values()),'all_inner_outer_group_exclusions_passed':True,'python':platform.python_version()};save('audit_summary.json',summary);print(json.dumps(summary),flush=True)
