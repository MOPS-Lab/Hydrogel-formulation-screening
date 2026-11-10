"""Reproducible provisional hydrogel model. Run: python train.py
Input: data/source_records.csv; no LAP/GO data are used for training.
G is a Maxwell-model-derived modulus, conditional on observed gel formation.
"""
from pathlib import Path
import json, warnings, hashlib, sys
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
from scipy.optimize import least_squares
from scipy.special import expit
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import GroupKFold, KFold
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn.exceptions import ConvergenceWarning
import sklearn, scipy, matplotlib

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
IONS=['Ca','Mg','Ba','Sr','Fe','Cu','Al']
POL=['ALG','TOCNF','SCLG']
RAW=['id']+['stock_'+p for p in POL]+['f_'+p for p in POL]+['LAP','GO','pH']+IONS+['G']
KEYS=RAW[1:-1]
NUM=['c_ALG','c_TOCNF','c_SCLG','logdose']
SEED=2026
warnings.filterwarnings('ignore',category=ConvergenceWarning)

def enrich(d):
    d=d.copy()
    for p in POL: d['c_'+p]=d['stock_'+p]*d['f_'+p]
    d['ion']=d[IONS].apply(lambda r: next((k for k,v in r.items() if v>0),'none'),axis=1)
    d['dose_M']=10*d[IONS].sum(axis=1)
    d['logdose']=np.log10(1+d.dose_M/.01)
    d['family']=d[['c_'+p for p in POL]].apply(lambda r:'+'.join(p for p in POL if r['c_'+p]>1e-10),axis=1)
    d['gel_formed']=(d.G>0).astype(int)
    d['G_model_Pa']=d.G.where(d.G>0)
    # Block all doses and preparation routes sharing the same effective polymer input.
    d['recipe_group']=pd.factorize(pd.MultiIndex.from_frame(d[['c_'+p for p in POL]].round(10)))[0]
    return d

def prepare():
    raw=pd.read_csv(ROOT/'data/source_records.csv')
    assert raw.id.notna().all() and raw.id.is_unique
    assert np.allclose(raw[['f_'+p for p in POL]+['LAP','GO']].sum(axis=1),1)
    noadd=raw[(raw.LAP==0)&(raw.GO==0)].copy()
    audit=raw.copy(); audit['decision']=''; audit['retained_id']=np.nan
    audit.loc[~audit.id.isin(noadd.id),'decision']='excluded_LAP_GO'
    noadd['record_ids']=noadd.groupby(KEYS+['G'],dropna=False).id.transform(lambda s:','.join(str(int(v)) for v in s))
    dedup=noadd.drop_duplicates(KEYS+['G']).copy()
    first=noadd.groupby(KEYS+['G'],dropna=False).id.transform('first')
    audit.loc[noadd.index,'retained_id']=first.to_numpy()
    duplicate=noadd.id!=first
    audit.loc[noadd.index[duplicate],'decision']='exact_duplicate'
    conflict=dedup.duplicated(KEYS,keep=False)
    conflicts=dedup[conflict].copy()
    clean=dedup[~conflict].copy()
    audit.loc[audit.id.isin(conflicts.id),'decision']='conflicting_G_quarantined'
    audit.loc[audit.id.isin(clean.id),'decision']='retain_gel'
    audit.loc[audit.id.isin(clean.loc[clean.G==0,'id']),'decision']='retain_no_gel'
    assert (audit.decision!='').all()
    clean=enrich(clean)
    cgroups=[]
    for group, g in conflicts.groupby(KEYS,sort=False,dropna=False):
        vals=g.G.to_numpy(); rel=(vals.max()-vals.min())/vals.min()
        cgroups.append(dict(ids=','.join(g.record_ids),values_Pa='; '.join(f'{v:.9g}' for v in vals),ratio=float(vals.max()/vals.min()),note='Razlika < 0,1 %; možna razlika v zaokroževanju.' if rel<.001 else 'Različni G pri enakih zapisanih vhodih; potrebno preverjanje.',**{k:g.iloc[0][k] for k in KEYS}))
    pd.DataFrame(cgroups).to_csv(OUT/'conflicts.csv',index=False)
    audit.to_csv(OUT/'cleaning_audit.csv',index=False)
    clean.to_csv(OUT/'clean_records.csv',index=False)
    return raw,clean,audit,pd.DataFrame(cgroups)

def make_model(kind='GP'):
    prep=ColumnTransformer([('num',StandardScaler(),NUM),('ion',OneHotEncoder(handle_unknown='ignore',sparse_output=False),['ion'])])
    if kind=='GP':
        estimator=GaussianProcessRegressor(kernel=ConstantKernel(1,(.01,100))*Matern(1,(.03,100),nu=2.5)+WhiteKernel(.03,(1e-5,1)),normalize_y=True,n_restarts_optimizer=2,random_state=SEED)
    else: estimator=ExtraTreesRegressor(n_estimators=400,min_samples_leaf=2,random_state=SEED,n_jobs=1)
    return make_pipeline(prep,estimator)

def predict_gp(model,d):
    return model[-1].predict(model[0].transform(d),return_std=True)

def metrics(y,p):
    err=np.abs(y-p)
    return dict(n=len(y),R2_log=float(r2_score(y,p)),MAE_log=float(err.mean()),median_factor=float(10**np.median(err)),geomean_factor=float(10**err.mean()),MAE_Pa=float(mean_absolute_error(10**y,10**p)),within_factor2=float(np.mean(err<=np.log10(2))))

def empirical_q(scores,coverage=.9):
    # Finite sample order statistic; empirical calibration, no group-exchangeability guarantee.
    a=np.sort(scores); k=min(len(a),int(np.ceil((len(a)+1)*coverage)))
    return float(a[k-1])

def validate(d,regime,kind='GP',nested=False):
    d=d.reset_index(drop=True); y=np.log10(d.G.to_numpy()); p=np.zeros(len(d)); sd=np.full(len(d),np.nan); q=np.full(len(d),np.nan); folds=np.zeros(len(d),int)
    split=GroupKFold(5).split(d,groups=d.recipe_group) if 'blocked' in regime else KFold(5,shuffle=True,random_state=SEED).split(d)
    for fold,(tr,te) in enumerate(split,1):
        m=make_model(kind).fit(d.iloc[tr],y[tr])
        if kind=='GP': p[te],sd[te]=predict_gp(m,d.iloc[te])
        else: p[te]=m.predict(d.iloc[te])
        folds[te]=fold
        if nested:
            dt=d.iloc[tr].reset_index(drop=True); yt=y[tr]; score=[]
            for itr,ite in GroupKFold(4).split(dt,groups=dt.recipe_group):
                mi=make_model().fit(dt.iloc[itr],yt[itr]); mp,ms=predict_gp(mi,dt.iloc[ite]); score.extend(np.abs(yt[ite]-mp)/np.maximum(ms,1e-8))
            q[te]=empirical_q(score)
    res=metrics(y,p); res.update(regime=regime,model=kind)
    if nested:
        res.update(interval_target=.9,interval_coverage=float(np.mean(abs(y-p)<=q*sd)),median_interval_factor=float(np.median(10**(q*sd))))
    tab=pd.DataFrame(dict(id=d.id,record_ids=d.record_ids,ion=d.ion,family=d.family,G_Pa=d.G,logG=y,pred_logG=p,pred_Pa=10**p,std_logG=sd,fold=folds,q90=q))
    if nested: tab['lower90_Pa']=10**(p-q*sd); tab['upper90_Pa']=10**(p+q*sd)
    tab.to_csv(OUT/f'oof_{regime}_{kind}.csv',index=False)
    print(json.dumps(res),flush=True)
    return res,tab

def export_gp(m):
    pre=m[0]; gp=m[-1]; sc=pre.named_transformers_['num']; cats=pre.named_transformers_['ion'].categories_[0].tolist()
    return dict(numeric_features=NUM,scaler_mean=sc.mean_.tolist(),scaler_scale=sc.scale_.tolist(),ion_categories=cats,X=gp.X_train_.tolist(),alpha=gp.alpha_.tolist(),L=gp.L_.tolist(),amplitude=float(gp.kernel_.k1.k1.constant_value),length_scale=float(gp.kernel_.k1.k2.length_scale),noise=float(gp.kernel_.k2.noise_level),y_mean=float(gp._y_train_mean),y_std=float(gp._y_train_std))

def domains(d):
    out={}
    for (family,ion),g in d.groupby(['family','ion']):
        a=g[NUM].to_numpy(); origin=a.mean(axis=0); _,s,v=np.linalg.svd(a-origin,full_matrices=False); rank=int(np.sum(s>1e-8)); basis=v[:rank].T; z=(a-origin)@basis
        item=dict(n=len(g),origin=origin.tolist(),basis=basis.tolist(),rank=rank)
        if rank==1: item.update(lower=float(z.min()),upper=float(z.max()))
        elif rank>1: item['equations']=ConvexHull(z).equations.tolist()
        out[family+'|'+ion]=item
    return out

def equations(d):
    t=d[(d.f_ALG==1)&(d.ion=='Ca')]; p=t.stock_ALG.to_numpy(); c=t.dose_M.to_numpy(); y=np.log(t.G.to_numpy())
    def f(theta,p,c):
        a,n,k,h=theta; return a+n*np.log(p)+np.log(expit(h*(np.log(c)-k)))
    def fit(p,c,y): return least_squares(lambda t:f(t,p,c)-y,[np.log(20000),1.8,np.log(.05),2.5],bounds=([0,0,np.log(.0001),.1],[np.log(1e7),6,np.log(10),10])).x
    pars=fit(p,c,y); cv=np.zeros(len(y))
    for conc in sorted(set(p)):
        te=p==conc; cv[te]=f(fit(p[~te],c[~te],y[~te]),p[te],c[te])
    r=metrics(y/np.log(10),cv/np.log(10)); r['n_samples']=r.pop('n'); r.update(A_Pa=float(np.exp(pars[0])),n=float(pars[1]),K_M=float(np.exp(pars[2])),h=float(pars[3]),formula='G = A*(cALG/1%)^n * (cCa/K)^h / (1+(cCa/K)^h)',validation='Leave one full ALG concentration out (4 folds); pure ALG + Ca only; no gel threshold inferred.')
    return r

def experiment_plan():
    rows=[]
    slices=[('ALG 2 %',2,0,0,1,0,0),('TOCNF 1 %',0,1,0,0,1,0),('ALG/SCLG 75/25',2,0,2,.75,0,.25),('ALG/SCLG 90/10',2,0,2,.9,0,.1)]
    for name,sa,st,ss,fa,ft,fs in slices:
        for dose in [0,.003,.01,.03,.180189620758483]:
            rows.append(dict(phase=1,goal='Meja nastanka gela',series=name,stock_ALG=sa,stock_TOCNF=st,stock_SCLG=ss,f_ALG=fa,f_TOCNF=ft,f_SCLG=fs,ion='Ca' if dose else 'none',dose_M=dose,independent_preparations=3,note='Isti postopek priprave; 0 in zgornja koncentracija sta kontrolna pogoja.'))
    for fa,ft,fs in [(.25,.75,0),(.5,.5,0),(.75,.25,0),(1/3,1/3,1/3),(.6,.2,.2),(.2,.6,.2),(.2,.2,.6)]:
        for dose in [.03,.180189620758483]:
            rows.append(dict(phase=2,goal='Manjkajoče interakcije',series='ALG/TOCNF/SCLG',stock_ALG=2,stock_TOCNF=2,stock_SCLG=2 if fs else 0,f_ALG=fa,f_TOCNF=ft,f_SCLG=fs,ion='Ca',dose_M=dose,independent_preparations=3,note='Predlog za drugo fazo; odmerke prilagoditi po prvi fazi. Napoved G za te kombinacije še ni podprta.'))
    p=pd.DataFrame(rows); p.insert(0,'plan_id',[f'E{i+1:02d}' for i in range(len(p))]); p.to_csv(OUT/'experiment_plan.csv',index=False)

def main():
    raw,clean,audit,conflicts=prepare(); pos=clean[clean.G>0].copy(); ca=pos[pos.ion.isin(['Ca','none'])].copy()
    results=[]
    for kind in ['GP','ExtraTrees']:
        result,tab=validate(pos,'all_random',kind); results.append(result)
    for kind in ['GP','ExtraTrees']:
        result,tab=validate(ca,'calcium_blocked',kind,nested=kind=='GP'); results.append(result)
        if kind=='GP': ca_oof=tab
    gp_ca=make_model().fit(ca,np.log10(ca.G)); gp_all=make_model().fit(pos,np.log10(pos.G))
    q=empirical_q(abs(ca_oof.logG-ca_oof.pred_logG)/ca_oof.std_logG)
    bundle=dict(version='1.0',date='2026-09-29',source_sha256=hashlib.sha256((ROOT/'data/source_records.csv').read_bytes()).hexdigest(),definition='G from generalized Maxwell-model fitting; conditional on observed gel formation; reference DOI 10.1016/j.ijbiomac.2020.05.088.',units='Polymer percentages and fractions as recorded; derived nominal pre-crosslinking component concentration = stock * fraction. Ion dose_M = 10 * recorded mol/100 ml; not internal/free-ion concentration.',scope='Neutral pH, original preparation and fitting protocol. LAP/GO excluded. No trained gel classifier: only four negative observations.',models={'calcium':export_gp(gp_ca),'all_exploratory':export_gp(gp_all)},q90_calcium=q,domains=domains(pos),observations=clean[['id','c_ALG','c_TOCNF','c_SCLG','ion','dose_M','G','gel_formed']].to_dict('records'),metrics=results,hill=equations(clean),counts=dict(source=len(raw),excluded_additives=int((audit.decision=='excluded_LAP_GO').sum()),exact_duplicates=int((audit.decision=='exact_duplicate').sum()),conflict_records=int((audit.decision=='conflicting_G_quarantined').sum()),conflict_recipes=len(conflicts),retained=len(clean),gels=len(pos),no_gel=int((clean.G==0).sum())))
    (OUT/'model.json').write_text(json.dumps(bundle,ensure_ascii=False,separators=(',',':')))
    (OUT/'metrics.json').write_text(json.dumps(results,indent=2))
    experiment_plan()
    (ROOT/'requirements.txt').write_text(f'numpy=={np.__version__}\npandas=={pd.__version__}\nscipy=={scipy.__version__}\nscikit-learn=={sklearn.__version__}\nmatplotlib=={matplotlib.__version__}\n')
    print('Finished',bundle['counts'],'q90',q,flush=True)

if __name__=='__main__': main()
