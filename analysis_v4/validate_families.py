import json
import numpy as np
from compare_models import *
family=np.array([r['family'] for r in rows]);out=[]
for fam in sorted(set(family)):
 tr=np.flatnonzero(family!=fam);te=np.flatnonzero(family==fam)
 p=tune('joint_gp',tr);pred=base[te]+fit_predict('joint_gp',p,tr,te)
 for i,v in zip(te,pred):out.append(dict(id=rows[i]['id'],family=fam,observed=rows[i]['G_observed_Pa'],G_Pa=float(10**v),legacy_overlap=bool(rows[i]['legacy_matches'])))
novel=[r for r in out if not r['legacy_overlap']]
result=dict(description='All KS samples of one polymer family held out; legacy retained except compositions matching held-out KS records. Diagnostic only, not entirely unseen family when legacy has it.',metrics=metrics(np.log10([r['observed'] for r in novel]),np.log10([r['G_Pa'] for r in novel])),by_family={f:metrics(np.log10([r['observed'] for r in novel if r['family']==f]),np.log10([r['G_Pa'] for r in novel if r['family']==f])) for f in sorted(set(r['family'] for r in novel)) if len([r for r in novel if r['family']==f])>1},rows=out)
(ROOT/'family_validation.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
