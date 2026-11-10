"""Checks input identity, leakage exclusions, deployed GP and reported scores."""
from pathlib import Path
import json,csv,sys,math
import numpy as np
R=Path(__file__).resolve().parent
sys.path.insert(0,str(R/'analysis_source/model'));from hydrogel import HydrogelModel
sys.path.insert(0,str(R/'analysis_v5'));import train_v5 as v5
def load(p):return json.loads((R/p).read_text())
rs=load('data/curated/records.json');byid={r['id']:r for r in rs}
assert len(rs)==205 and len(byid)==205 and sum(bool(r['gel_formed']) for r in rs)==198
assert rs==load('analysis_source/records.json')==load('analysis_v5/records.json')
assert len([r for r in rs if r['source']=='legacy' and r['gel_formed'] and r['ion'] in ['Fe','Al']])==12
assert len({tuple(r['c']) for r in rs if r['source']=='legacy' and r['gel_formed'] and r['ion'] in ['Mg','Ba','Sr','Cu']})==1
group=lambda rid:tuple(round(v,8) for v in byid[rid]['c'])
folds=load('revision_analysis/fold_manifest.json');tested=[]
for f in folds:
    train=set(f['training_ids']);test=set(f['held_ids']);tested+=list(test)
    assert train.isdisjoint(test)
    assert {group(i) for i in train}.isdisjoint({group(i) for i in test})
    for inner in f['inner']:
        a=set(inner['training_ids']);b=set(inner['validation_ids'])
        assert a<=train and b<=train and a.isdisjoint(b)
        assert {group(i) for i in a}.isdisjoint({group(i) for i in b|test})
assert len(tested)==34 and len(set(tested))==34 and len(folds)==16
refs=load('tests/reference_results.json')['revision_analysis/new_series_metrics.json']
pred=list(csv.DictReader((R/'revision_analysis/new_series_predictions.csv').open()))
y=np.log10([float(r['observed_Pa']) for r in pred]);assert len(pred)==30
max_score_error=0.
for name,ref in refs.items():
    p=np.log10([float(r[name]) for r in pred]);score=v5.metrics(y,p)
    for k in ['RMSE_log10','MAE_log10','median_factor','within_factor2','R2_log']:
        err=abs(score[k]-ref[k]);max_score_error=max(max_score_error,err);assert err<1e-9,(name,k,err)
m=HydrogelModel();idx=np.arange(len(v5.ks));_,param,_=v5.tune('v4_gp',idx);fit=v5.fit('v4_gp',param,idx)
assert list(param)==['raw',4,.3,.9]
refit=v5.predict(fit,v5.old+v5.ks)
deployed=np.array([m.log_mean(r['c'],r['ion'],r['dose_M'],source=r['source']) for r in v5.old+v5.ks])
max_fit_error=float(np.max(abs(refit-deployed)));assert max_fit_error<1e-10
for r in rs:
    if not r['gel_formed']:assert m.predict_effective(r['c'],r['ion'],r['dose_M'])['G_Pa'] is None
for ion in ['Ca','Fe','Al']:
    a=m.predict([1.7]*3,[.65,.35,0],ion,.06);b=m.predict_effective([.5525,.2975,0],ion,.03)
    assert math.isclose(a['G_Pa'],b['G_Pa'],rel_tol=1e-12)
cal=load('revision_analysis/mcmc_calibration.json');assert next(x for x in cal if x['nominal_percent']==80)['covered_records']==18
di=load('revision_analysis/mcmc_diagnostics_recomputed.json')
assert len(di)==17 and all(v['rhat']<1.01 and v['ess_bulk']>=400 for f in di for v in f['parameters'].values())
search=load('revision_analysis/search_results.json');cand=next(x for x in search if x['target_Pa']==15000)['recipes'][0]
sys.path.insert(0,str(R/'analysis_source'));from monte_carlo import simulate
q=cand['recipe'];sim=simulate(stocks=q['stocks'],fractions=q['fractions'],ion=q['ion'],dose=q['dose'],target_Pa=15000,n=5000,seed=20261014)
assert np.allclose(sim['preparation90_Pa'],cand['mc']['preparation90_Pa'],rtol=1e-11)
assert sim['outside_training_fraction']==.8138
out={'records_verified':205,'positive_records':198,'known_failures':7,'outer_folds_checked':16,'inner_folds_checked':64,'headline_records':30,'max_metric_difference_from_20261101_reference':max_score_error,'max_full_refit_log10_difference':max_fit_error,'Monte_Carlo_cross_language_quantiles_match':True,'candidate15_outside_support_fraction':sim['outside_training_fraction'],'MCMC_80pct_coverage':'18/30','scope':'Available analyses verified; original historical seven-family benchmark not rerun; no raw rheology refit.'}
(R/'results').mkdir(exist_ok=True);(R/'results/verification.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
