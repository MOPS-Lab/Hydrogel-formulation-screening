"""Reproducible preparation Monte Carlo using physical dispensed volumes.
Nominal polymer mixture : crosslinker solution = 1 : 1. Default ranges are
independent uniform SCENARIOS, not measured tolerances. Model uncertainty is
excluded. See mcmc_gp.py for the separately validated hyperparameter experiment.
"""
from pathlib import Path
import json,math,argparse,sys
import numpy as np
from scipy.spatial.distance import cdist
R=Path(__file__).resolve().parent;sys.path.insert(0,str(R/'model'));from hydrogel import HydrogelModel
class RNG:
 def __init__(self,seed):self.x=seed&0xffffffff
 def uniform(self):
  self.x=(self.x+0x6D2B79F5)&0xffffffff;t=self.x;t=((t^(t>>15))*(t|1))&0xffffffff;t^=(t+(((t^(t>>7))*(t|61))&0xffffffff))&0xffffffff;return ((t^(t>>14))&0xffffffff)/4294967296

def draws(stocks,fractions,dose,n=5000,seed=20261014,stock_halfwidth_pp=.1,volume_halfwidth_percent=2,crosslinker_halfwidth_percent=5):
 random=RNG(seed);out=[]
 for _ in range(n):
  u=np.array([2*random.uniform()-1 for i in range(8)]);ss=np.array(stocks)+stock_halfwidth_pp*u[:3];vol=np.array(fractions)*(1+volume_halfwidth_percent*u[3:6]/100);vx=1+volume_halfwidth_percent*u[6]/100;total=vol.sum()+vx;ionstock=dose*(1+crosslinker_halfwidth_percent*u[7]/100)
  if np.any((np.array(fractions)>0)&(ss<=0)) or total<=0 or ionstock<0:out.append(None)
  else:out.append(dict(c=(ss*vol/total).tolist(),dose=float(ionstock*vx/total),ratio=float(vol.sum()/vx)))
 return out

def simulate(stocks=(2,2,2),fractions=(.5,.5,0),ion='Ca',dose=.09,n=5000,seed=20261014,stock_halfwidth_pp=.1,volume_halfwidth_percent=2,crosslinker_halfwidth_percent=5,target_Pa=15000,tolerance_percent=10):
 if not isinstance(n,int) or not 100<=n<=50000 or not isinstance(seed,int) or not 0<=seed<=4294967295:raise ValueError('Invalid n/seed')
 if not all(math.isfinite(v) for v in [stock_halfwidth_pp,volume_halfwidth_percent,crosslinker_halfwidth_percent,target_Pa,tolerance_percent]) or not 0<=stock_halfwidth_pp<=1 or not 0<=volume_halfwidth_percent<100 or not 0<=crosslinker_halfwidth_percent<100 or target_Pa<=0 or not 0<tolerance_percent<=100:raise ValueError('Invalid scenarios/target')
 model=HydrogelModel();nom=model.predict(stocks,fractions,ion,dose)
 if nom['G_Pa'] is None:return dict(status='observed_no_gel')
 ion=nom['ion'];m=model.m;rs=draws(stocks,fractions,dose,n,seed,stock_halfwidth_pp,volume_halfwidth_percent,crosslinker_halfwidth_percent);fails=[r for r in model.meta['records'] if not r['gel_formed'] and r['ion']==ion];valid=[];invalid=0;no_gel=0
 for r in rs:
  if r is None:invalid+=1;continue
  if any(abs(f['dose_M']-r['dose'])<1e-8 and all(abs(v-w)<1e-8 for v,w in zip(f['c'],r['c'])) for f in fails):no_gel+=1;continue
  valid.append(r)
 if not valid:return dict(status='no_valid_simulations',invalid_n=invalid,known_no_gel_n=no_gel)
 X=np.array([[*r['c'],np.log10(1+r['dose']/.005)] for r in valid]);X=(X-m['scaler_mean'])/m['scaler_scale'];t=np.sqrt(5)*cdist(X,m['X'])/m['length'];K=(1+t+t*t/3)*np.exp(-t)*np.array([1 if ion==a else m['rho_ion'] for a in m['ions']])*np.array([1 if ('legacy' if ion=='none' else 'KS')==a else m['rho_source'] for a in m['sources']]);G=10**(m['y_mean']+K@np.array(m['alpha']));p=float(np.mean(np.abs(G/target_Pa-1)<=tolerance_percent/100));z=1.959963984540054;den=1+z*z/len(G);center=(p+z*z/(2*len(G)))/den;half=z*math.sqrt(p*(1-p)/len(G)+z*z/(4*len(G)**2))/den;out=dict(target_fraction_numerical95=[max(0,center-half),min(1,center+half)],status='ok',n=n,seed=seed,valid_n=len(valid),invalid_n=invalid,known_no_gel_n=no_gel,mean_Pa=float(G.mean()),median_Pa=float(np.median(G)),sd_Pa=float(G.std(ddof=1)),preparation90_Pa=np.quantile(G,[.05,.95]).tolist(),target_fraction_given_valid=p,target_fraction_MC_se=math.sqrt(p*(1-p)/len(G)),outside_training_fraction=float(np.mean([not model.inside(r['c'],ion,r['dose']) for r in valid])),nominal_G_Pa=nom['G_Pa'],conditional_on_gel_formation=True,includes_model_uncertainty=False,scenario={'stock_halfwidth_pp':stock_halfwidth_pp,'volume_halfwidth_percent':volume_halfwidth_percent,'crosslinker_halfwidth_percent':crosslinker_halfwidth_percent});return out
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stocks',nargs=3,type=float,default=[2,2,2]);p.add_argument('--fractions',nargs=3,type=float,default=[.5,.5,0]);p.add_argument('--ion',default='Ca');p.add_argument('--dose',type=float,default=.09,help='crosslinker STOCK mol/l');p.add_argument('--n',type=int,default=5000);p.add_argument('--seed',type=int,default=20261014);p.add_argument('--target',type=float,default=15000);p.add_argument('--tolerance',type=float,default=10);a=p.parse_args();print(json.dumps(simulate(stocks=a.stocks,fractions=a.fractions,ion=a.ion,dose=a.dose,n=a.n,seed=a.seed,target_Pa=a.target,tolerance_percent=a.tolerance),indent=2))
