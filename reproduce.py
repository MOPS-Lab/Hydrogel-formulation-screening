"""Portable entry point; analyses write to revision_analysis/ and results/."""
from pathlib import Path
import argparse,subprocess,sys,os,json,platform,shutil
R=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--mode',choices=['check','all','mcmc-fresh','figures'],default='check')
a=p.parse_args()
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
def run(script,node=False):
    print('Running',script,flush=True)
    subprocess.run(['node' if node else sys.executable,str(R/script)],cwd=R,env=env,check=True)
if a.mode=='all':
    run('prepare_data.py')
    run('reproduce_historical_gp.py')
    run('analysis_v5/train_v5.py')
    run('audit_revision.py')
    run('run_search.cjs',True)
    run('mc_sensitivity.py')
elif a.mode=='mcmc-fresh':
    dest=R/'results/mcmc_fresh'
    if dest.exists():raise SystemExit('results/mcmc_fresh already exists. Rename it before a new uncached sampling run.')
    (dest/'mcmc').mkdir(parents=True);(dest/'model').mkdir()
    for f in ['mcmc_gp.py','records.json']:shutil.copy2(R/'analysis_source'/f,dest/f)
    run('results/mcmc_fresh/mcmc_gp.py')
    raise SystemExit(0)
if a.mode in ['all','figures']:run('revision_figures.py')
if a.mode in ['all','check']:
    if a.mode=='check':run('prepare_data.py')
    run('export_tables.py');run('verify_release.py')
    if shutil.which('node'):run('tests/check_app.cjs',True)
    else:print('Node unavailable: JavaScript verification skipped. Browser app itself needs no Node installation.')
    import numpy,scipy,sklearn,pandas,matplotlib,openpyxl
    versions={'python':platform.python_version(),**{m.__name__:m.__version__ for m in [numpy,scipy,sklearn,pandas,matplotlib,openpyxl]}}
    if shutil.which('node'):versions['node']=subprocess.check_output(['node','--version'],text=True).strip()
    (R/'results/execution_environment.json').write_text(json.dumps(versions,indent=2))
