"""Run all finite checks and regenerate the complete reproducibility data.

Run: python verification/verify_all.py
All output paths are relative to this file's parent project, not the cwd.
"""
from contextlib import redirect_stdout
from datetime import datetime,timezone
from pathlib import Path
import hashlib,json,platform,time
import networkx,sympy
from diagram_checks import ROOT,generate_cores,check_local
from marking_counts import check_counts
from long_checks import run as check_long
from metric_checks import run as check_metrics
from global_checks import check_faces,check_algebra

def main():
    (ROOT/'tables').mkdir(exist_ok=True)
    (ROOT/'logs').mkdir(exist_ok=True)
    start=time.perf_counter()
    results={}
    with (ROOT/'logs/verification_run.log').open('w',encoding='utf-8') as log:
        with redirect_stdout(log):
            cores=generate_cores()
            (ROOT/'tables/minimal_core_catalog.json').write_text(json.dumps([
                dict(name=name,n=G[0],edges=G[1]) for name,G in cores],indent=2)+'\n')
            results['local']=check_local(cores)
            results['marking_counts']=check_counts(cores)
            results['type_IV']=check_long()
            results['metric']=check_metrics()
            results['small_three_faces']=check_faces()
            results['algebra']=check_algebra()
    summary=dict(status='passed',utc=datetime.now(timezone.utc).isoformat(),
                 elapsed_seconds=round(time.perf_counter()-start,3),python=platform.python_version(),
                 networkx=networkx.__version__,sympy=sympy.__version__,results=results,
                 scope='Exact finite computations and symbolic identities; geometric reductions and classification completeness are proved or cited in the manuscript.')
    (ROOT/'logs/verification_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    manifest={}
    files=[*sorted((ROOT/'verification').glob('*.py')),*sorted((ROOT/'verification').glob('*.json')),
           *sorted((ROOT/'tables').glob('*.json')),*sorted((ROOT/'tables').glob('*.csv')),
           ROOT/'requirements.txt']
    for path in files:
        manifest[path.relative_to(ROOT).as_posix()]=hashlib.sha256(path.read_bytes()).hexdigest()
    (ROOT/'logs/SHA256SUMS.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
