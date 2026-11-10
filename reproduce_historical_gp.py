"""Recover the Ca/no-ion GP benchmark only. Other original families are missing."""
from pathlib import Path
import sys,json
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parent
sys.path.insert(0,str(R/'analysis_historical_gp'))
import train_original as original
_,clean,_,_=original.prepare()
ca=clean[(clean.G>0)&clean.ion.isin(['Ca','none'])].copy()
metrics,tab=original.validate(ca,'calcium_blocked','GP',nested=True)
reference=pd.read_csv(R/'analysis_historical_gp/archived_results/oof_calcium_blocked_GP.csv').set_index('id')
current=tab.set_index('id').loc[reference.index]
assert np.array_equal(current.fold,reference.fold)
maxdiff=float(np.max(abs(current.pred_logG-reference.pred_logG)))
rmse=float(np.sqrt(np.mean((current.logG-current.pred_logG)**2)))
assert abs(rmse-.134773550357283)<1e-5,rmse
out={'n':len(ca),'groups':int(ca.recipe_group.nunique()),'RMSE_log10':rmse,'max_abs_log10_prediction_difference_from_archived':maxdiff,'exact_original_fold_membership_matches':True,'metrics':metrics,'scope':'Historical subset A Gaussian process only; not the other six A families or three B models.'}
(R/'results/historical_GP_reproduction.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
