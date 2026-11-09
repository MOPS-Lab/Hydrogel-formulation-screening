"""Historical cleaning recovered from v1. Not the seven-family benchmark."""
from pathlib import Path
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'revision_analysis'/'reconstructed_data'
OUT.mkdir(parents=True,exist_ok=True)
IONS=['Ca','Mg','Ba','Sr','Fe','Cu','Al']
POL=['ALG','TOCNF','SCLG']
RAW=['id']+['stock_'+p for p in POL]+['f_'+p for p in POL]+['LAP','GO','pH']+IONS+['G']
KEYS=RAW[1:-1]
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
    raw=pd.read_csv(ROOT/'data/raw/historical_source_records.csv')
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
