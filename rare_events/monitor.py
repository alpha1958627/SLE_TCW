"""Read only one explicit pilot job and its exact directory; no restart/submission."""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import subprocess
import time
from . import core


def command(args):
    p=subprocess.run(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    return dict(returncode=p.returncode,stdout=p.stdout.strip(),stderr=p.stderr.strip())


def snapshot(root,job):
    root=Path(root).resolve();campaign=core.load(root/'campaign.json');assert campaign['job']==job
    runs={}
    for mode in ('ams','dmc','qdmc','direct'):
        base=root/mode
        if not (base/'manifest.json').exists():continue
        manifest=core.load(base/'manifest.json');assert int(manifest['slurm_job'])==job
        for name,sha in manifest['source_sha256'].items():
            assert hashlib.sha256((core.ROOT.parent/name).read_bytes()).hexdigest()==sha,name+' changed'
        assert manifest['build']==core.build_identity()
        calls=[]
        for request in sorted(base.rglob('request.json')):
            case=request.parent;spec=core.load(request)
            exe=core.load(case/'execution.json') if (case/'execution.json').exists() else None
            summary=core.load(case/'summary.json') if (case/'summary.json').exists() else None
            rows=[]
            if (case/'trajectory.tsv').exists():
                for line in (case/'trajectory.tsv').read_text().splitlines()[1:]:
                    fields=line.split('\t')
                    if len(fields)==7:rows.append([float(v) for v in fields])
            calls.append(dict(case=case.relative_to(root).as_posix(),request=spec,execution=exe,summary=summary,
                observations=rows,active_elapsed_core_hours=(time.time()-request.stat().st_mtime)*spec['ranks']/3600 if not exe else 0.))
        def read(name):return core.load(base/name) if (base/name).exists() else None
        state=read('state.json');compact=None
        if state:
            if mode in ('dmc','qdmc'):
                compact=dict(stage=state['stage'],logz=state['logz'],statistics=state.get('statistics'),
                    log_compensation=[state['logz']-p['V'] for p in state['particles']])
            elif mode=='ams':compact=dict(round=state['round'],log_weight=state['log_weight'],events=sum(p['tau'] is not None for p in state['particles']))
            else:compact=dict(events=sum(p['tau'] is not None for p in state['particles']))
        runs[mode]=dict(manifest=manifest,status=read('status.json'),state=compact,result=read('result.json'),
                       selections=[core.load(p) for p in sorted(base.glob('selection_*.json'))],calls=calls)
    return dict(job=job,utc=datetime.datetime.utcnow().isoformat()+'Z',root=str(root),runs=runs,
        campaign_status=core.load(root/'campaign_status.json') if (root/'campaign_status.json').exists() else None,
        squeue=command(['squeue','-j',str(job),'-h','-o','%i|%T|%M|%C|%R']),
        sacct=command(['sacct','-j',str(job),'-X','-n','-P','--format=JobIDRaw,State,ElapsedRaw,AllocCPUS,ExitCode']))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workdir',required=True);p.add_argument('--job',type=int,required=True)
    a=p.parse_args();print(json.dumps(snapshot(a.workdir,a.job),indent=2))
