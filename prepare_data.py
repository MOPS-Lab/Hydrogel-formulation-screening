"""Reconstruct the numerical inputs from recovered tabular sources and verify them.

No rheology fitting is performed. The imported Maxwell-derived G is the target.
"""
from pathlib import Path
import json,csv,math
import numpy as np
import pandas as pd
from openpyxl import load_workbook
from historical_preprocessing import prepare
R=Path(__file__).resolve().parent
def main():
    raw,clean,audit,conflicts=prepare()
    records=json.loads((R/'data/curated/records.json').read_text())
    byid={r['id']:r for r in records}
    for r in clean.to_dict('records'):
        saved=byid['OLD_'+str(float(r['id']))]
        assert np.allclose(saved['c'],[r['c_'+p]/2 for p in ['ALG','TOCNF','SCLG']],rtol=0,atol=1e-12)
        assert math.isclose(saved['dose_M'],r['dose_M']/2,rel_tol=1e-12,abs_tol=1e-14)
        assert saved['G_observed_Pa']==r['G'] and saved['ion']==r['ion']
    w=load_workbook(R/'data/raw/new_series_original.xlsx',data_only=True)
    tab=list(w.active.values);provenance=json.loads((R/'data/curated/new_series_provenance.json').read_text())
    def num(x):return float(str(x).replace(',','.')) if x is not None else 0.
    transformed=[]
    for p in provenance:
        line=tab[p['row']-1];rid='KS_'+str(int(num(line[0])))
        stock=[num(v) for v in line[1:4]];fractions=[num(v) for v in line[4:7]];correction=''
        if rid=='KS_23':fractions=[.57,.43,0.];correction='0.43 assigned to TOCNF; confirmed correction in v4'
        if rid=='KS_6':fractions=[v/sum(fractions) for v in fractions];correction='0.33/0.33/0.33 normalised to equal thirds'
        c=[s*f/2 for s,f in zip(stock,fractions)]
        ioncols=[num(v) for v in line[7:14]];active=[i for i,v in enumerate(ioncols) if v>0]
        assert len(active)<=1,'Mixtures require a separate encoding'
        ion=['Ca','Mg','Ba','Sr','Fe','Cu','Al'][active[0]] if active else 'none'
        dose=10*sum(ioncols)/2;observed=num(line[15]);saved=byid[rid]
        assert np.allclose(c,saved['c'],rtol=0,atol=1e-12)
        assert math.isclose(dose,saved['dose_M'],rel_tol=1e-12,abs_tol=1e-14)
        assert observed==saved['G_observed_Pa'] and ion==saved['ion']
        transformed.append(dict(id=rid,source_row=p['row'],response_cell='P'+str(p['row']),stock_ALG_pct=stock[0],stock_TOCNF_pct=stock[1],stock_SCLG_pct=stock[2],f_ALG=fractions[0],f_TOCNF=fractions[1],f_SCLG=fractions[2],ion=ion,crosslinker_stock_M=2*dose,nominal_final_M=dose,G_Pa=observed,correction=correction))
    assert len(clean)==168 and len(transformed)==37
    out=R/'revision_analysis/reconstructed_data';out.mkdir(exist_ok=True,parents=True)
    with (out/'new_series_inputs.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=transformed[0]);writer.writeheader();writer.writerows(transformed)
    summary={'raw_historical_rows':len(raw),'historical_retained':len(clean),'new_retained':len(transformed),'all_205_records_verified':True,'historical_cleaning_decisions':audit.decision.value_counts().to_dict(),'new_response_column':'P (tt)','ion_conversion':'recorded mol/100mL * 10 / 2 = nominal final mol/L; independent of valence'}
    (out/'verification.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
