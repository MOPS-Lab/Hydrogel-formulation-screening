"""Export machine-readable manuscript tables from saved analysis outputs."""
from pathlib import Path
import json,csv
from collections import Counter
from openpyxl import load_workbook
R=Path(__file__).resolve().parent;A=R/'revision_analysis';O=R/'results/tables';O.mkdir(parents=True,exist_ok=True)
def load(p):return json.loads((R/p).read_text())
def save(name,rows):
    with (O/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
rs=load('data/curated/records.json')
hist=[r for r in rs if r['source']=='legacy' and r['gel_formed'] and r['ion'] not in ['Fe','Al']]
counts=Counter(r['ion'] for r in hist)
save('table_1_historical_B_ion_counts.csv',[dict(ion=i,n=counts[i]) for i in ['Ca','none','Mg','Ba','Sr','Cu']])
book=load_workbook(R/'data/raw/historical_benchmark_report.xlsx',data_only=True)
for n,sheet in zip([2,3],list(book)[:2]):
    rows=list(sheet.values);hdr=rows[25];data=[{hdr[i]:row[i] for i in range(11)} for row in rows[26:] if row[1] is not None]
    save(f'table_{n}_historical_reported_metrics.csv',data)
metrics=load('revision_analysis/new_series_metrics.json');mm=load('analysis_source/mcmc_validation.json')['metrics_novel']
rows=[]
for model,m in metrics.items():rows.append(dict(model=model,n=m['n'],RMSE_log10=m['RMSE_log10'],median_factor=m['median_factor'],within_factor2=m['within_factor2'],R2_log=m['R2_log'],status='reported in Table 4' if model not in ['RBF_SVR_added','hist_gradient_boosting_added'] else 'reported in abstract; absent from current Table 4'))
pred=load('analysis_source/mcmc_validation.json')['rows'];import numpy as np
keep=[x for x in pred if not x['legacy_overlap']];y=np.log10([x['observed'] for x in keep]);p=np.log10([x['G_Pa'] for x in keep]);r2=1-np.sum((y-p)**2)/np.sum((y-y.mean())**2)
rows.append(dict(model='GP_MCMC',n=mm['n'],RMSE_log10=mm['RMSE_log10'],median_factor=mm['median_factor'],within_factor2=mm['within_factor2'],R2_log=float(r2),status='reported in Table 4'))
save('table_4_extended_new_series_metrics.csv',rows)
save('table_5_independent_repeat_controls.csv',load('analysis_v5/diagnostics.json')['replicates'])
rows=[]
for search in load('revision_analysis/search_results.json'):
    if search['target_Pa'] not in [8000,10000,12000,15000]:continue
    r=search['recipes'][0];p=r['prediction'];q=r['recipe'];mc=r['mc']
    rows.append(dict(target_Pa=search['target_Pa'],stock_ALG_pct=q['stocks'][0],stock_TOCNF_pct=q['stocks'][1],stock_SCLG_pct=q['stocks'][2],f_ALG=q['fractions'][0],f_TOCNF=q['fractions'][1],f_SCLG=q['fractions'][2],ion=q['ion'],crosslinker_stock_M=q['dose'],ion_final_M=p['dose_M'],c_ALG_final_pct=p['c_percent'][0],c_TOCNF_final_pct=p['c_percent'][1],c_SCLG_final_pct=p['c_percent'][2],nominal_G_Pa=p['G_Pa'],MC_p05_Pa=mc['preparation90_Pa'][0],MC_p95_Pa=mc['preparation90_Pa'][1],scenario_target_fraction=mc['target_fraction_given_valid'],outside_support_fraction=mc['outside_training_fraction'],nominal_inside_combined_support=p['inside_combined_input_domain'],nominal_inside_new_support=p['inside_KS_input_domain'],MC_N=mc['valid_n'],seed=mc['config']['seed'],experimentally_tested=False))
save('table_6_candidates_full_precision.csv',rows)
rows=[]
for fold in load('revision_analysis/fold_manifest.json'):
    for role,key in [('train','training_ids'),('test','held_ids')]:
        for rid in fold[key]:rows.append(dict(outer_fold=fold['fold'],role=role,record_id=rid))
save('outer_fold_membership.csv',rows)
print('Exported manuscript tables and outer-fold membership.')
