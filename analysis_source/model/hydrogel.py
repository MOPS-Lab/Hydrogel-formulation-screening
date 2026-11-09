"""Exploratory v4 source-aware Gaussian process for Maxwell-derived G (Pa).
Forward prediction needs only Python standard library. See README.txt for scope.
"""
from pathlib import Path
import math,json,argparse
ROOT=Path(__file__).resolve().parent
class HydrogelModel:
 def __init__(self):
  self.m=json.loads((ROOT/'joint_gp.json').read_text());self.meta=json.loads((ROOT/'metadata.json').read_text());self.diag=json.loads((ROOT/'diagnostics.json').read_text());self.alt=json.loads((ROOT/'alternative.json').read_text())
 def log_mean(self,c,ion,dose,source=None):
  source=source or ('legacy' if ion=='none' else 'KS');m=self.m;x=[*c,math.log10(1+dose/.005)];x=[(v-a)/b for v,a,b in zip(x,m['scaler_mean'],m['scaler_scale'])];s=0.
  for row,a,i,b in zip(m['X'],m['alpha'],m['ions'],m['sources']):
   t=math.sqrt(5*sum((v-u)**2 for v,u in zip(x,row)))/m['length'];s+=(1+t+t*t/3)*math.exp(-t)*(1 if ion==i else m['rho_ion'])*(1 if source==b else m['rho_source'])*a
  return m['y_mean']+s
 def alternative_log(self,c,ion,dose):
  C=sum(c);w=[v/C for v in c];d=math.log10(1+dose/.005);lc=math.log10(C);iv=[float(x==ion) for x in ['Ca','Mg','Ba','Sr','Fe','Cu','Al','none']];src=float(ion!='none');x=[lc,d,w[0],w[1],*iv,src,*[v*src for v in iv],lc*w[0],lc*w[1],d*w[0],d*w[1]];m=self.alt
  return m['intercept']+sum((v-a)/b*k for v,a,b,k in zip(x,m['scaler_mean'],m['scaler_scale'],m['coef']))
 def evidence(self,c,ion,dose,G):
  family='+'.join(p for p,v in zip(['ALG','TOCNF','SCLG'],c) if v>1e-10);count=sum(r['source']=='KS' and r['gel_formed'] and r['ion']==ion and '+'.join(p for p,v in zip(['ALG','TOCNF','SCLG'],r['c']) if v>1e-10)==family for r in self.meta['records']);v=self.diag['by_family_ion'].get(family+'|'+ion);scope='ista družina in ion'
  if not v or v['n']<5 or v['n_compositions']<3:v=self.diag['by_family'].get(family);scope='ista družina, različni ioni'
  if not v or v['n']<5 or v['n_compositions']<3:v=self.diag['all'];scope='celoten nabor novih sestav'
  alt=10**self.alternative_log(c,ion,dose)
  return dict(validation_reference=dict(**v,scope=scope,reference_range_Pa=[G/v['factor80'],G*v['factor80']],definition='Descriptive 80th percentile of held-out errors, not a calibrated prediction interval') if count else None,alternative_G_Pa=alt,model_disagreement_factor=max(G/alt,alt/G),positive_KS_family_ion_count=count)
 def inside(self,c,ion,dose,ks=False):
  family='+'.join(p for p,v in zip(['ALG','TOCNF','SCLG'],c) if v>1e-10);g=self.meta['ks_domains' if ks else 'domains'].get(family+'|'+ion)
  if g is None:return False
  d=[v-z for v,z in zip([*c,math.log10(1+dose/.005)],g['origin'])];rank=g['rank']
  if not rank:return math.hypot(*d)<1e-7
  z=[sum(d[i]*g['basis'][i][j] for i in range(4)) for j in range(rank)]
  if math.hypot(*[d[i]-sum(g['basis'][i][j]*z[j] for j in range(rank)) for i in range(4)])>1e-7:return False
  if rank==1:return g['lower']-1e-7<=z[0]<=g['upper']+1e-7
  return all(sum(e[i]*z[i] for i in range(rank))+e[-1]<=1e-7 for e in g['equations'])
 def predict_effective(self,c,ion='Ca',dose=.045):
  if len(c)!=3 or not all(math.isfinite(v) and v>=0 for v in c) or not 0<sum(c)<=1.5+1e-8 or not math.isfinite(dose) or dose<0:raise ValueError('Invalid composition/dose')
  if ion not in self.meta['support']:raise ValueError('Unknown ion')
  if ion=='none' and dose!=0:raise ValueError('No ion requires zero dose')
  if dose==0:ion='none'
  q=self.meta['support'][ion]
  if not q['dose_min']-1e-9<=dose<=q['dose_max']+1e-9:raise ValueError('Dose outside data range')
  ids=[r['id'] for r in self.meta['records'] if not r['gel_formed'] and r['ion']==ion and abs(r['dose_M']-dose)<1e-8 and all(abs(v-w)<1e-8 for v,w in zip(c,r['c']))]
  inside=self.inside(c,ion,dose)
  common=dict(version=self.meta['version'],c_percent=c,total_percent=sum(c),ion=ion,dose_M=dose,prediction_source='legacy' if ion=='none' else 'KS',inside_combined_input_domain=inside,inside_KS_input_domain=self.inside(c,ion,dose,True),conditional_on_gel_formation=True,gel_probability=None,prediction_interval=None)
  if ids:return dict(**common,G_Pa=None,status='observed_no_gel',known_failure_ids=ids)
  if ion=='none' and not inside:raise ValueError('No-ion extrapolation unsupported')
  G=10**self.log_mean(c,ion,dose);return dict(**common,G_Pa=G,status='conditional_estimate',**self.evidence(c,ion,dose,G))
 def predict(self,stocks=(2,2,2),fractions=(.5,.5,0),ion='Ca',dose=.09):
  if len(stocks)!=3 or len(fractions)!=3 or not all(math.isfinite(v) and .1<=v<=3 for v in stocks) or not all(math.isfinite(v) and 0<=v<=1 for v in fractions) or abs(sum(fractions)-1)>1e-8:raise ValueError('Stocks .1-3%; fractions 0-1, sum 1')
  return self.predict_effective([a*b/2 for a,b in zip(stocks,fractions)],ion,dose/2)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--c',nargs=3,type=float,default=[.5,.5,0],help='final ALG TOCNF SCLG concentrations, percentage points');p.add_argument('--ion',default='Ca');p.add_argument('--dose',type=float,default=.045);a=p.parse_args();print(json.dumps(HydrogelModel().predict_effective(a.c,a.ion,a.dose),indent=2))
