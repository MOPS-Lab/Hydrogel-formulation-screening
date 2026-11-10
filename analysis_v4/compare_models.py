"""Nested composition-blocked validation of KS-only calibration of frozen v2.
No KS target is used to train its outer-fold prediction. Zero-G rows are not
regression targets. All rows of an effective polymer composition share a fold.
"""
from pathlib import Path
import json,math,warnings
import numpy as np
from sklearn.model_selection import GroupKFold,LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from scipy.spatial.distance import cdist
from scipy.linalg import cho_factor,cho_solve
ROOT=Path(__file__).resolve().parent
data=json.loads((ROOT/'ks_records.json').read_text());rows=[r for r in data if r['gel_formed']]
ions=['Ca','Mg','Ba','Sr','Fe','Cu','Al']
groups=np.array([r['group'] for r in rows]);y=np.log10([r['G_observed_Pa'] for r in rows]);base=np.log10([r['G_Pa_full_precision'] for r in rows]);resid=y-base
legacy=json.loads((ROOT/'previous_v2/base_model.json').read_text())
old=[dict(c=[r['c_ALG']/2,r['c_TOCNF']/2,r['c_SCLG']/2],ion=r['ion'],dose_M=r['dose_M']/2,G_observed_Pa=r['G'],source='legacy') for r in legacy['observations'] if r['gel_formed']]

def joint_inputs(rs):
 return np.array([[*r['c'],np.log10(1+r['dose_M']/.005)] for r in rs])

def joint_fit(tr,te,param,return_model=False):
 held_groups={rows[i]['group'] for i in range(len(rows)) if i not in set(tr)};train=[r for r in old if '|'.join(f'{v:.8f}' for v in r['c']) not in held_groups]+[rows[i] for i in tr];test=[rows[i] for i in te];sc=StandardScaler().fit(joint_inputs(train));X=sc.transform(joint_inputs(train));T=sc.transform(joint_inputs(test));it=np.array([r['ion'] for r in train]);isrc=np.array([r.get('source','KS') for r in train]);jt=np.array([r['ion'] for r in test]);jsrc=np.array(['KS']*len(test));length,noise,rho=param
 def kern(A,B,ia,ib,sa,sb):
  t=np.sqrt(5)*cdist(A,B)/length
  return (1+t+t*t/3)*np.exp(-t)*(.3+.7*(ia[:,None]==ib[None,:]))*(rho+(1-rho)*(sa[:,None]==sb[None,:]))
 Y=np.log10([r['G_observed_Pa'] for r in train]);mean=float(Y.mean());K=kern(X,X,it,it,isrc,isrc)+noise*np.eye(len(X));cf=cho_factor(K,lower=True);alpha=cho_solve(cf,Y-mean);p=mean+kern(T,X,jt,it,jsrc,isrc)@alpha
 if return_model:return dict(type='joint_gp',length=length,noise=noise,rho_source=rho,rho_ion=.3,scaler_mean=sc.mean_.tolist(),scaler_scale=sc.scale_.tolist(),X=X.tolist(),ions=it.tolist(),sources=isrc.tolist(),alpha=alpha.tolist(),y_mean=mean,L=np.tril(cf[0]).tolist())
 return p-base[te]

def features(rs,kind):
 out=[]
 for r in rs:
  c=np.array(r['c']);C=c.sum();w=c/C;d=np.log10(1+r['dose_M']/.005);ion=np.array([float(r['ion']==k) for k in ions]);pair=np.array([4*w[0]*w[1],4*w[0]*w[2],4*w[1]*w[2]])
  if kind=='ion':f=ion
  elif kind=='main':f=np.r_[np.log10(C),d,w,ion]
  elif kind=='pair':f=np.r_[np.log10(C),d,w,pair,ion]
  elif kind=='interaction':f=np.r_[np.log10(C),d,w,pair,ion,np.outer(ion,w).ravel()]
  elif kind=='smooth':f=np.r_[np.log10(C),d,w,ion]
  out.append(f)
 return np.array(out)

def fit_predict(kind,param,tr,te,return_model=False):
 if kind=='joint_gp':return joint_fit(tr,te,param,return_model)
 if kind=='base':return np.zeros(len(te))
 if kind=='offset':return np.full(len(te),np.mean(resid[tr]))
 if kind=='target_mean':return np.full(len(te),np.mean(y[tr]))-base[te]
 if kind.startswith('direct_'):
  feat=kind[7:];X=features(rows,'main' if feat=='calibrated' else feat)
  if feat=='calibrated':X=np.column_stack([X,base])
  sc=StandardScaler().fit(X[tr]);m=Ridge(alpha=param).fit(sc.transform(X[tr]),y[tr]);p=m.predict(sc.transform(X[te]))
  if return_model:return dict(type=kind,alpha=param,scaler_mean=sc.mean_.tolist(),scaler_scale=sc.scale_.tolist(),coef=m.coef_.tolist(),intercept=float(m.intercept_))
  return p-base[te]
 if kind.startswith('ridge_'):
  X=features(rows,kind[6:]);sc=StandardScaler().fit(X[tr]);m=Ridge(alpha=param).fit(sc.transform(X[tr]),resid[tr]);p=m.predict(sc.transform(X[te]))
  if return_model:return dict(type=kind,alpha=param,scaler_mean=sc.mean_.tolist(),scaler_scale=sc.scale_.tolist(),coef=m.coef_.tolist(),intercept=float(m.intercept_))
  return p
 if kind=='kernel_residual':
  X=features(rows,'smooth');sc=StandardScaler().fit(X[tr]);Z=sc.transform(X[tr]);T=sc.transform(X[te]);length,alpha=param
  def kernel(a,b):
   t=np.sqrt(5)*cdist(a,b)/length;return (1+t+t*t/3)*np.exp(-t)
  center=float(np.mean(resid[tr]));dual=np.linalg.solve(kernel(Z,Z)+alpha*np.eye(len(tr)),resid[tr]-center);p=center+kernel(T,Z)@dual
  if return_model:return dict(type=kind,length=length,alpha=alpha,scaler_mean=sc.mean_.tolist(),scaler_scale=sc.scale_.tolist(),X=Z.tolist(),dual=dual.tolist(),intercept=center)
  return p
 raise ValueError(kind)

SPECS={'base':[None],'offset':[None],'ridge_ion':[.1,1,10,100],
 'ridge_main':[.1,1,10,100],'ridge_pair':[.1,1,10,100],
 'ridge_interaction':[.1,1,10,100],'kernel_residual':[(l,a) for l in [1,3,10] for a in [.1,1]]}
SPECS.update({'target_mean':[None],**{k:[.1,1,10,100] for k in ['direct_ion','direct_main','direct_pair','direct_calibrated']}})
SPECS['joint_gp']=[(l,n,r) for l in [.5,1.5,4] for n in [.03,.3] for r in [.3,.9]]

def tune(kind,tr):
 candidates=SPECS[kind]
 if len(candidates)==1:return candidates[0]
 n=min(4,len(set(groups[tr])));best=None
 for param in candidates:
  errs=[]
  for itr,ite in GroupKFold(n).split(tr,groups=groups[tr]):
   aa=tr[itr];bb=tr[ite];pred=base[bb]+fit_predict(kind,param,aa,bb);errs.extend((y[bb]-pred)**2)
  score=np.mean(errs)
  if best is None or score<best[0]:best=(score,param)
 return best[1]

def metrics(actual,pred):
 e=np.abs(actual-pred)
 return dict(n=len(actual),RMSE_log10=float(np.sqrt(np.mean((actual-pred)**2))),MAE_log10=float(e.mean()),median_factor=float(10**np.median(e)),geomean_factor=float(10**e.mean()),within_factor2=float(np.mean(e<=math.log10(2))),R2_log=float(r2_score(actual,pred)),MAE_Pa=float(np.mean(np.abs(10**actual-10**pred))))

def main():
 oof={k:np.zeros(len(rows)) for k in SPECS};chosen={k:[] for k in SPECS}
 for fi,(tr,te) in enumerate(LeaveOneGroupOut().split(y,groups=groups),1):
  for kind in SPECS:
   param=tune(kind,tr);oof[kind][te]=base[te]+fit_predict(kind,param,tr,te);chosen[kind].append({'fold':fi,'held_ids':[rows[i]['id'] for i in te],'param':param})
 novel=np.array([not bool(r['legacy_matches']) for r in rows]);metrics_all={k:metrics(y,oof[k]) for k in oof};metrics_novel={k:metrics(y[novel],oof[k][novel]) for k in oof}
 for k,m in metrics_novel.items():print(k,json.dumps(m),flush=True)
 output=dict(validation='Nested leave-one-KS-composition-out; inner grouped 4-fold hyperparameter selection; frozen v2 unchanged; headline excludes 4 same-recipe independent controls; joint GP also removes legacy records matching held-out compositions.',n_groups=len(set(groups)),metrics_all=metrics_all,metrics_novel=metrics_novel,hyperparameters_by_fold=chosen,rows=[dict(id=r['id'],group=r['group'],G_observed_Pa=r['G_observed_Pa'],legacy_overlap=bool(r['legacy_matches']),predictions={k:float(10**oof[k][i]) for k in oof}) for i,r in enumerate(rows)])
 (ROOT/'model_comparison.json').write_text(json.dumps(output,indent=2))
 models={}
 for kind in SPECS:
  if kind not in ('base','offset','target_mean'):
   param=tune(kind,np.arange(len(rows)));models[kind]=fit_predict(kind,param,np.arange(len(rows)),np.arange(len(rows)),return_model=True)
 (ROOT/'fitted_candidates.json').write_text(json.dumps(models,indent=2))
if __name__=='__main__':main()
