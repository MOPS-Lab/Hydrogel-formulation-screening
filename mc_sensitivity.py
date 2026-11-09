from pathlib import Path
import sys,json,csv
import numpy as np
from scipy.spatial.distance import cdist
R=Path(__file__).resolve().parent;O=R/'revision_analysis'
sys.path.insert(0,str(R/'analysis_source'));from monte_carlo import draws,RNG
sys.path.insert(0,str(R/'analysis_source/model'));from hydrogel import HydrogelModel
model=HydrogelModel();m=model.m
def evaluate(rr,ion='Ca'):
 X=np.array([[*r['c'],np.log10(1+r['dose']/.005)] for r in rr]);T=(X-m['scaler_mean'])/m['scaler_scale'];t=np.sqrt(5)*cdist(T,m['X'])/m['length'];K=(1+t+t*t/3)*np.exp(-t)*np.array([1 if ion==i else m['rho_ion'] for i in m['ions']])*np.array([1 if i=='KS' else m['rho_source'] for i in m['sources']]);return X,10**(m['y_mean']+K@np.array(m['alpha']))
search=json.loads((O/'search_results.json').read_text());candidates=[x['recipes'][0] for x in search if x['target_Pa'] in [8000,10000,12000,15000]];results=[];geom=[]
for r in candidates:
 recipe=r['recipe'];target=r['mc']['config']['target_Pa'];st,fr,dose=recipe['stocks'],recipe['fractions'],recipe['dose']
 for scenario,factor,corr,seed,n in [('nominal',1,False,20261014,5000),('half_widths',.5,False,20261014,5000),('double_widths',2,False,20261014,5000),('shared_calibration_error',1,True,20261014,5000),('second_seed',1,False,20261101,5000),('larger_MC',1,False,20261014,20000),('polymer_only',1,False,20261014,5000),('volumes_only',1,False,20261014,5000),('crosslinker_only',1,False,20261014,5000)]:
  sw=.1*factor if scenario not in ['volumes_only','crosslinker_only'] else 0;vw=2*factor if scenario not in ['polymer_only','crosslinker_only'] else 0;dw=5*factor if scenario not in ['polymer_only','volumes_only'] else 0
  rr=draws(st,fr,dose,n,seed,sw,vw,dw)
  if corr:
   rng=RNG(seed);rr=[]
   for _ in range(n):
    u=np.array([2*rng.uniform()-1 for i in range(8)]);stocks=np.array(st)+.1*u[0];vol=np.array(fr)*(1+.02*u[3]);vx=1+.02*u[3];total=vol.sum()+vx;dd=dose*(1+.05*u[7])*vx/total;rr.append({'c':(stocks*vol/total).tolist(),'dose':dd})
  X,G=evaluate(rr);lo,med,hi=np.quantile(G,[.05,.5,.95]);out=np.array([not model.inside(x['c'],'Ca',x['dose']) for x in rr]);results.append(dict(target_kPa=target/1000,scenario=scenario,n=n,seed=seed,median_Pa=med,p05_Pa=lo,p95_Pa=hi,relative_halfwidth_percent=100*(hi-lo)/(2*med),target_fraction=np.mean(abs(G/target-1)<=.1),outside_fraction=out.mean()))
  if scenario=='nominal':
   dom=model.meta['domains']['ALG+TOCNF|Ca'];basis=np.array(dom['basis']);D=X-dom['origin'];Z=D@basis;aff=np.linalg.norm(D-Z@basis.T,axis=1)>1e-7;eq=np.array(dom['equations']);fac=Z@eq[:,:-1].T+eq[:,-1];x0=np.array([*r['prediction']['c_percent'],np.log10(1+r['prediction']['dose_M']/.005)]);z0=(x0-dom['origin'])@basis
   geom.append(dict(target_kPa=target/1000,rank=dom['rank'],n_support=dom['n'],affine_violations=int(aff.sum()),hull_violations=int(np.any(fac>1e-7,axis=1).sum()),per_facet_violation_counts=(fac>1e-7).sum(axis=0).tolist(),nominal_signed_facet_values=(z0@eq[:,:-1].T+eq[:,-1]).tolist(),tolerance=1e-7))
   if target==15000:np.savez_compressed(O/'candidate15_MC.npz',G=G,X=X)
for name,data in [('mc_sensitivity.json',results),('support_geometry.json',geom)]: (O/name).write_text(json.dumps(data,indent=2,default=lambda x:x.item()))
print(json.dumps([r for r in results if r['target_kPa']==15],indent=2,default=lambda x:x.item()))
