"""Build the English offline application; original v6 numerical engines unchanged."""
from pathlib import Path
import json
R=Path(__file__).resolve().parent;M=R/'analysis_source/model'
def load(p):return json.loads(p.read_text())
core='\n'.join((M/n).read_text() for n in ['hydrogel.js','monte_carlo.js','mcmc_predict.js'])
bundle={k:load(M/n) for k,n in [('gp','joint_gp.json'),('meta','metadata.json'),('diagnostics','diagnostics.json'),('alternative','alternative.json'),('mcmc','mcmc_model.json')]}
bundle['metrics']=load(R/'revision_analysis/new_series_metrics.json');bundle['csv']=(R/'data/curated/records.csv').read_text()
html=(R/'app/template.html').read_text().replace('/*__CORE__*/',core).replace('/*__BUNDLE__*/',json.dumps(bundle,separators=(',',':'),ensure_ascii=False).replace('</','<\\/'))
(R/'index.html').write_text(html)
print('Built index.html:',len(html.encode()),'bytes')
