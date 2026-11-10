"""Model comparison with nested, composition-blocked KS validation.
All legacy rows with held-out KS compositions are removed at every training split.
Fixed candidate menu, MSE(log10 G) inner selection; no outer target in tuning.
"""
from pathlib import Path
import json,math,itertools
import numpy as np
from scipy.spatial.distance import cdist
from scipy.linalg import cho_factor,cho_solve
from sklearn.model_selection import GroupKFold,LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
R=Path(__file__).resolve().parent
allrecords=json.loads((R/'records.json').read_text());ks=[r for r in allrecords if r['source']=='KS' and r['gel_formed']];old=[r for r in allrecords if r['source']=='legacy' and r['gel_formed']]
ions=['Ca','Mg','Ba','Sr','Fe','Cu','Al','none']
group=lambda r:tuple(round(v,8) for v in r['c'])
g=np.array(['|'.join(f'{v:.8f}' for v in r['c']) for r in ks]);y=np.log10([r['G_observed_Pa'] for r in ks])

def make_features(rs,kind):
 out=[]
 for r in rs:
  c=np.array(r['c']);C=c.sum();w=c/C;dose=np.log10(1+r['dose_M']/.005);ion=np.array([float(r['ion']==i) for i in ions]);source=float(r['source']=='KS');pair=np.array([4*w[0]*w[1],4*w[0]*w[2],4*w[1]*w[2]])
  basic=np.r_[np.log10(C),dose,w[:2]]
  if kind=='raw':f=np.r_[c,dose]
  elif kind=='mixture':f=basic
  elif kind=='ion':f=ion
  elif kind=='main':f=np.r_[basic,ion]
  elif kind=='joint':f=np.r_[basic,ion,source,source*ion]
  elif kind=='pair':f=np.r_[basic,pair,ion,source,source*ion]
  elif kind=='interactions':f=np.r_[basic,pair,ion,np.outer(ion,w[:2]).ravel(),source,source*ion]
  elif kind=='concentration':f=np.r_[basic,ion,source,source*ion,np.log10(C)*w[:2],dose*w[:2]]
  out.append(f)
 return np.array(out)

def trainset(tr):
 # Includes neither outer nor inner held-out compositions, even from legacy.
 allowed={group(ks[i]) for i in tr};held={group(r) for r in ks}-allowed
 return [r for r in old if group(r) not in held]+[ks[i] for i in tr]

def linear_fit(train,kind,alpha,old_weight):
 X=make_features(train,kind);sc=StandardScaler().fit(X);Y=np.log10([r['G_observed_Pa'] for r in train]);weights=np.array([old_weight if r['source']=='legacy' else 1 for r in train]);reg=Ridge(alpha=alpha).fit(sc.transform(X),Y,sample_weight=weights)
 return dict(type='linear',features=kind,scaler_mean=sc.mean_.tolist(),scaler_scale=sc.scale_.tolist(),coef=reg.coef_.tolist(),intercept=float(reg.intercept_))

def gp_fit(train,kind,l,n,rho,mean_model=None):
 X=make_features(train,kind);sc=StandardScaler().fit(X);Z=sc.transform(X);ion=np.array([r['ion'] for r in train]);src=np.array([r['source'] for r in train]);Y=np.log10([r['G_observed_Pa'] for r in train]);center=float(Y.mean()) if mean_model is None else 0.;res=Y-center if mean_model is None else Y-predict(mean_model,train)
 t=np.sqrt(5)*cdist(Z,Z)/l;K=(1+t+t*t/3)*np.exp(-t)*(.3+.7*(ion[:,None]==ion[None,:]))*(rho+(1-rho)*(src[:,None]==src[None,:]))
 alpha=cho_solve(cho_factor(K+n*np.eye(len(Z)),lower=True),res)
 return dict(type='gp',features=kind,scaler_mean=sc.mean_.tolist(),scaler_scale=sc.scale_.tolist(),X=Z.tolist(),ions=ion.tolist(),sources=src.tolist(),alpha=alpha.tolist(),length=l,noise=n,rho_source=rho,rho_ion=.3,intercept=center,mean_model=mean_model)

def predict(m,rs):
 if m['type']=='ensemble':return sum(w*predict(mm,rs) for w,mm in zip(m['weights'],m['models']))
 X=(make_features(rs,m['features'])-m['scaler_mean'])/m['scaler_scale']
 if m['type']=='linear':return X@np.array(m['coef'])+m['intercept']
 ion=np.array([r['ion'] for r in rs]);src=np.array([r['source'] for r in rs]);t=np.sqrt(5)*cdist(X,m['X'])/m['length'];K=(1+t+t*t/3)*np.exp(-t)*(.3+.7*(ion[:,None]==np.array(m['ions'])[None,:]))*(m['rho_source']+(1-m['rho_source'])*(src[:,None]==np.array(m['sources'])[None,:]));out=m['intercept']+K@np.array(m['alpha']);return out if m['mean_model'] is None else out+predict(m['mean_model'],rs)

# Modest, declared grids. More data matter more than a very large search.
SPECS={
 'v4_gp': [('raw',l,n,r) for l in [.5,1.5,4] for n in [.03,.3] for r in [.3,.9]],
 'ks_main': [('main',a,0.) for a in [1,10,100]],
 'joint_ridge': [(f,a,w) for f in ['joint','pair','concentration','interactions'] for a in [1,10,100] for w in [.25,1]],
 'mixture_gp': [('mixture',l,n,r) for l in [1.5,4] for n in [.1,.5] for r in [.3,.9]],
 'trend_gp': [('mixture',l,n,r) for l in [1.5,4] for n in [.1,.5] for r in [.3,.9]],
 'ensemble': [(w,) for w in [.25,.5,.75]]
}
cache={}
def fit(kind,param,tr):
 key=(kind,tuple(param),tuple(sorted(tr)))
 if key in cache:return cache[key]
 train=trainset(tr)
 if kind in ('v4_gp','mixture_gp'):m=gp_fit(train,*param)
 elif kind=='trend_gp':m=gp_fit(train,*param,mean_model=linear_fit(train,'joint',10,.25))
 elif kind=='joint_ridge':m=linear_fit(train,*param)
 elif kind=='ks_main':m=linear_fit([ks[i] for i in tr],*param)
 elif kind=='ensemble':
  # Components use fixed regularization, weight alone selected inside folds.
  m=dict(type='ensemble',weights=[param[0],1-param[0]],models=[gp_fit(train,'raw',4,.3,.9),linear_fit([ks[i] for i in tr],'main',10,1)])
 cache[key]=m;return m

def tune(kind,tr):
 best=None
 for p in SPECS[kind]:
  out=[];err=[]
  for a,b in GroupKFold(min(4,len(set(g[tr])))).split(tr,groups=g[tr]):
   aa=tr[a];bb=tr[b];pr=predict(fit(kind,p,aa),[ks[i] for i in bb]);err.extend((y[bb]-pr)**2);out.extend(zip(bb.tolist(),pr.tolist()))
  score=float(np.mean(err))
  if best is None or score<best[0]:best=(score,p,out)
 return best

def metrics(a,p):
 e=np.abs(a-p);return dict(n=len(a),RMSE_log10=float(np.sqrt(np.mean((a-p)**2))),MAE_log10=float(e.mean()),median_factor=float(10**np.median(e)),within_factor2=float(np.mean(e<=np.log10(2))),R2_log=float(1-np.sum((a-p)**2)/np.sum((a-a.mean())**2)))

def main():
 oof={k:np.zeros(len(ks)) for k in SPECS};oof['nested_selector']=np.zeros(len(ks));models=[];folds=[]
 for no,(tr,te) in enumerate(LeaveOneGroupOut().split(y,groups=g),1):
  tuned={k:tune(k,tr) for k in SPECS}
  chosen=min(tuned,key=lambda k:tuned[k][0]);fold={'fold':no,'held_ids':[ks[i]['id'] for i in te],'chosen':chosen,'params':{k:list(v[1]) for k,v in tuned.items()}}
  for k,(_,p,_) in tuned.items():oof[k][te]=predict(fit(k,p,tr),[ks[i] for i in te])
  oof['nested_selector'][te]=oof[chosen][te];models.append(fit(chosen,tuned[chosen][1],tr));folds.append(fold)
  print('fold',no,'of 16',chosen,flush=True);cache.clear()
 novel=np.array([not r.get('legacy_matches') for r in ks]);out={'description':'Nested leave-one-KS-composition-out, 4-fold inner tuning and algorithm selection; matching legacy compositions excluded; all KS positive records tested; subset excludes 4 near-repeat controls.','metrics_all':{k:metrics(y,v) for k,v in oof.items()},'metrics_novel':{k:metrics(y[novel],v[novel]) for k,v in oof.items()},'folds':folds,'rows':[dict(id=r['id'],group=g[i],family=r['family'],ion=r['ion'],observed=r['G_observed_Pa'],legacy_overlap=bool(r.get('legacy_matches')),predictions={k:float(10**v[i]) for k,v in oof.items()}) for i,r in enumerate(ks)]}
 (R/'comparison.json').write_text(json.dumps(out,indent=2));(R/'fold_models.json').write_text(json.dumps(models))
 tuned={k:tune(k,np.arange(len(ks))) for k in SPECS};chosen=min(tuned,key=lambda k:tuned[k][0]);m=fit(chosen,tuned[chosen][1],np.arange(len(ks)));(R/'full_models.json').write_text(json.dumps({k:fit(k,v[1],np.arange(len(ks))) for k,v in tuned.items()}));(R/'selection.json').write_text(json.dumps({'chosen':chosen,'inner_scores':{k:v[0] for k,v in tuned.items()},'params':{k:list(v[1]) for k,v in tuned.items()}},indent=2));print(json.dumps(out['metrics_novel'],indent=2),flush=True);print('full selected',chosen,flush=True)
if __name__=='__main__':main()
