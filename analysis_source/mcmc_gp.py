"""Hyperparameter MCMC for source-aware GP. Priors fixed before validation.
Posterior: length, observation NOISE VARIANCE, source correlation. GP mean is
training mean (empirical Bayes), amplitude=1 and ion correlation=.3 fixed.
4 adaptive RW Metropolis chains; adaptation only in warmup. Grouped validation
removes all matching legacy polymer compositions. No outer response enters fit.
"""
from pathlib import Path
import json,math,time,argparse
import numpy as np
from scipy.spatial.distance import cdist
from scipy.linalg import cho_factor,cho_solve,solve_triangular
from scipy.special import expit,ndtri
from scipy.stats import rankdata,norm
from scipy.optimize import brentq
from sklearn.preprocessing import StandardScaler
R=Path(__file__).resolve().parent
ALL=json.loads((R/'records.json').read_text());KS=[r for r in ALL if r['source']=='KS' and r['gel_formed']];OLD=[r for r in ALL if r['source']=='legacy' and r['gel_formed']]
def grp(r):return '|'.join(f'{v:.8f}' for v in r['c'])
def features(rs):return np.array([[*r['c'],np.log10(1+r['dose_M']/.005)] for r in rs])
PRIOR_MU=np.array([np.log(2),np.log(.1),np.log(.7/.3)]);PRIOR_SD=np.array([1.,1.2,1.5])
class GP:
 def __init__(self,records):
  self.rs=records;self.sc=StandardScaler().fit(features(records));self.X=self.sc.transform(features(records));self.y=np.log10([r['G_observed_Pa'] for r in records]);self.center=float(self.y.mean());self.y-=self.center;self.dist=cdist(self.X,self.X);self.ions=np.array([r['ion'] for r in records]);self.src=np.array([r['source'] for r in records]);self.ki=.3+.7*(self.ions[:,None]==self.ions[None,:]);self.same=self.src[:,None]==self.src[None,:]
 def factor(self,z):
  l,n=np.exp(z[:2]);rho=expit(z[2]);t=np.sqrt(5)*self.dist/l;K=(1+t+t*t/3)*np.exp(-t)*self.ki*(rho+(1-rho)*self.same);K.flat[::len(K)+1]+=n+1e-10;cf=cho_factor(K,lower=True,check_finite=False);a=cho_solve(cf,self.y,check_finite=False);return cf,a
 def logpost(self,z):
  if np.any(np.abs(z)>20):return -np.inf
  try:
   cf,a=self.factor(z);return float(-.5*self.y@a-np.log(np.diag(cf[0])).sum()-.5*np.sum(((z-PRIOR_MU)/PRIOR_SD)**2))
  except (np.linalg.LinAlgError,ValueError,FloatingPointError):return -np.inf
 def conditional(self,z,rs,export=False):
  cf,a=self.factor(z);l,n=np.exp(z[:2]);rho=expit(z[2]);T=self.sc.transform(features(rs));t=np.sqrt(5)*cdist(T,self.X)/l;ki=.3+.7*(np.array([r['ion'] for r in rs])[:,None]==self.ions[None,:]);ss=np.array([r['source'] for r in rs])[:,None]==self.src[None,:];K=(1+t+t*t/3)*np.exp(-t)*ki*(rho+(1-rho)*ss);mu=self.center+K@a;v=solve_triangular(cf[0],K.T,lower=True,check_finite=False);var=np.maximum(1-np.sum(v*v,axis=0),1e-12)
  if export:return dict(length=float(l),noise=float(n),rho_source=float(rho),alpha=a.tolist(),L_packed=cf[0][np.tril_indices(len(a))].tolist())
  return mu,var,var+n

def split(a):
 n=a.shape[1]//2;return np.concatenate([a[:,:n],a[:,-n:]],axis=0)
def rhat_basic(a):
 n=a.shape[1];W=np.var(a,axis=1,ddof=1).mean();B=n*np.var(a.mean(axis=1),ddof=1);return float(np.sqrt(((n-1)*W/n+B/n)/W))
def ess(a):
 # Multi-chain autocorrelation and initial positive/monotone paired sequence.
 m,n=a.shape;means=a.mean(axis=1);W=np.var(a,axis=1,ddof=1).mean();varplus=(n-1)/n*W+np.var(means,ddof=1);ac=[]
 for x in a:
  y=x-x.mean();fft=np.fft.rfft(y,n=2*n);ac.append(np.fft.irfft(fft*np.conj(fft))[:n]/n)
 ac=np.mean(ac,axis=0);rho=1-(W-ac)/varplus;rho[0]=1.;pairs=[]
 for k in range(0,n-1,2):
  v=float(rho[k]+rho[k+1])
  if v<0:break
  if pairs:v=min(v,pairs[-1])
  pairs.append(v)
 tau=max(1.,-1+2*sum(pairs));return float(min(m*n,m*n/tau))
def diagnostics(a):
 out={}
 for j,name in enumerate(['log_length','log_noise_variance','logit_source_correlation']):
  x=split(a[:,:,j]);z=ndtri((rankdata(x.ravel()).reshape(x.shape)-.375)/(x.size+.25));f=np.abs(x-np.median(x));fz=ndtri((rankdata(f.ravel()).reshape(x.shape)-.375)/(x.size+.25));out[name]=dict(rhat=max(rhat_basic(z),rhat_basic(fz)),ess_bulk=ess(z))
 return out

def sample(gp,seed,warmup=1000,draws=2500):
 chains=[];states=[];accept=[]
 for c in range(4):
  rng=np.random.default_rng(seed+104729*c);z=PRIOR_MU+np.array([[-1.,.8,-1.],[.8,-.8,1.],[-.5,-.5,.5],[1.,.7,-.5]])[c]*PRIOR_SD;lp=gp.logpost(z);hist=[];scale=.5;cov=np.diag([.3,.3,.6])**2;accepted=0;block=0
  for i in range(warmup):
   prop=z+rng.multivariate_normal(np.zeros(3),scale**2*cov);pp=gp.logpost(prop)
   if np.log(rng.uniform())<pp-lp:z=prop;lp=pp;block+=1
   hist.append(z.copy())
   if (i+1)%50==0:
    scale*=np.exp(.5*(block/50-.234));scale=np.clip(scale,.08,5);block=0
    if i>=199:cov=np.cov(np.array(hist[-500:]).T)+np.eye(3)*.005
  L=np.linalg.cholesky(cov)*scale;rr=[]
  for i in range(draws):
   prop=z+L@rng.normal(size=3);pp=gp.logpost(prop)
   if np.log(rng.uniform())<pp-lp:z=prop;lp=pp;accepted+=1
   rr.append(z.copy())
  chains.append(rr);states.append((rng,z,lp,L));accept.append(accepted/draws)
 arr=np.array(chains)
 # Extend stationary chains to resolve concrete convergence/ESS failure.
 for extension in range(3):
  d=diagnostics(arr)
  if all(v['rhat']<1.01 and v['ess_bulk']>=400 for v in d.values()):break
  added=[]
  for c,(rng,z,lp,L) in enumerate(states):
   rr=[]
   for i in range(draws):
    prop=z+L@rng.normal(size=3);pp=gp.logpost(prop)
    if np.log(rng.uniform())<pp-lp:z=prop;lp=pp
    rr.append(z.copy())
   states[c]=(rng,z,lp,L);added.append(rr)
  arr=np.concatenate([arr,np.array(added)],axis=1)
 d=diagnostics(arr);return arr,dict(parameters=d,chains=4,warmup_per_chain=warmup,retained_per_chain=arr.shape[1],acceptance_initial_sampling=accept,passed=all(v['rhat']<1.01 and v['ess_bulk']>=400 for v in d.values()))

def selected(arr,count=128):
 # Equal representation of all chains and evenly spaced retained draws.
 n=count//4;ii=np.linspace(0,arr.shape[1]-1,n,dtype=int);return arr[:,ii,:].reshape(-1,3)
def predictive(gp,draws,rs):
 mus=[];vs=[];obs=[]
 for z in draws:
  a,b,c=gp.conditional(z,rs);mus.append(a);vs.append(b);obs.append(c)
 mus=np.array(mus);vs=np.array(vs);obs=np.array(obs);out=[]
 for j,r in enumerate(rs):
  mm=mus[:,j];vv=obs[:,j];sd=np.sqrt(vv)
  def q(prob):return brentq(lambda v:np.mean(norm.cdf((v-mm)/sd))-prob,float(np.min(mm-12*sd)),float(np.max(mm+12*sd)))
  out.append(dict(id=r['id'],G_Pa=float(10**mm.mean()),G_median_predictive_Pa=float(10**q(.5)),predictive80_Pa=[float(10**q(.1)),float(10**q(.9))],sd_log_hyper_means=float(mm.std(ddof=1)),mean_latent_variance_log=float(vs[:,j].mean()),mean_noise_variance_log=float(np.mean(obs[:,j]-vs[:,j]))))
 return out

def metrics(rows):
 y=np.log10([r['observed'] for r in rows]);p=np.log10([r['G_Pa'] for r in rows]);e=np.abs(y-p);return dict(n=len(rows),RMSE_log10=float(np.sqrt(np.mean((y-p)**2))),median_factor=float(10**np.median(e)),within_factor2=float(np.mean(e<=np.log10(2))),predictive80_empirical_coverage=float(np.mean([r['predictive80_Pa'][0]<=r['observed']<=r['predictive80_Pa'][1] for r in rows])),median_predictive80_upper_lower_ratio=float(np.median([r['predictive80_Pa'][1]/r['predictive80_Pa'][0] for r in rows])))

def run(name,train,test,seed):
 gp=GP(train);cache=R/'mcmc'/f'{name}.npz';dfile=R/'mcmc'/f'{name}_diagnostics.json'
 if cache.exists():arr=np.load(cache)['draws'];d=json.loads(dfile.read_text())
 else:
  arr,d=sample(gp,seed);np.savez_compressed(cache,draws=arr);dfile.write_text(json.dumps(d,indent=2))
 preds=predictive(gp,selected(arr),test);print(name,'n',len(train),'samples',arr.shape,'maxR',round(max(v['rhat'] for v in d['parameters'].values()),4),'passed',d['passed'],flush=True)
 return gp,arr,d,preds

def main():
 started=time.time();train=OLD+KS;gp,arr,d,p=run('full',train,KS,6381)
 # Keep 32 representative hyperparameter draws for optional browser comparison.
 draws=selected(arr,32);bundle=dict(description='Exploratory MCMC; empirical training-mean and fixed amplitude=1, ion correlation=.3',scaler_mean=gp.sc.mean_.tolist(),scaler_scale=gp.sc.scale_.tolist(),X=gp.X.tolist(),ions=gp.ions.tolist(),sources=gp.src.tolist(),y_mean=gp.center,draws=[gp.conditional(z,[KS[0]],True) for z in draws],diagnostics=d,prior={'log_length_normal':[float(PRIOR_MU[0]),float(PRIOR_SD[0])],'log_noise_variance_normal':[float(PRIOR_MU[1]),float(PRIOR_SD[1])],'logit_rho_source_normal':[float(PRIOR_MU[2]),float(PRIOR_SD[2])]},browser_draws=32,validation_draws=128)
 (R/'model/mcmc_model.json').write_text(json.dumps(bundle,separators=(',',':')))
 out=[];folds=[]
 for i,g in enumerate(sorted(set(grp(r) for r in KS)),1):
  tr=[r for r in OLD+KS if grp(r)!=g];te=[r for r in KS if grp(r)==g];_,aa,dd,pp=run(f'fold_{i:02}',tr,te,6381+1000*i)
  for r,v in zip(te,pp):out.append(dict(**v,observed=r['G_observed_Pa'],legacy_overlap=bool(r.get('legacy_matches')),group=g,fold=i))
  folds.append(dict(fold=i,held_ids=[r['id'] for r in te],diagnostics=dd))
  (R/'mcmc/progress.json').write_text(json.dumps(dict(rows=out,folds=folds),indent=2))
 novel=[r for r in out if not r['legacy_overlap']];result=dict(model='MCMC over length/noise/source correlation; no outer-fold tuning',full_diagnostics=d,folds=folds,rows=out,metrics_all=metrics(out),metrics_novel=metrics(novel),all_chains_passed=d['passed'] and all(f['diagnostics']['passed'] for f in folds),seconds=time.time()-started)
 (R/'mcmc_validation.json').write_text(json.dumps(result,indent=2));print('DONE',json.dumps(result['metrics_novel']),flush=True)
if __name__=='__main__':main()
