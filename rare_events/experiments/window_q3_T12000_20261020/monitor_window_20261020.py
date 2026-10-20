"""One read-only snapshot of the four exact independently submitted window jobs."""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import sys
import time

def snapshot(repo):
    repo=Path(repo).resolve();sys.path.insert(0,str(repo))
    from rare_events import core
    from rare_events.monitor import command
    folder=Path(__file__).resolve().parent;ledger=core.load(folder/'submission.json');plan=core.load(folder/'plan.json')
    assert ledger['status']=='submitted' and set(ledger['jobs'])==set(plan['configs'])
    build=core.build_identity();runs={};lifecycle={}
    for mode,job in ledger['jobs'].items():
        work=repo/'experiments'/plan['experiment']/mode/('job_%d'%job);base=work/'production'
        def read(path):return core.load(path) if path.exists() else None
        lifecycle[mode]=dict(job=job,workdir=str(work),experiment_manifest=read(work/'experiment_manifest.json'),startup_gate=read(work/'startup_gate.json'))
        manifest=read(base/'manifest.json')
        if not manifest:continue
        assert int(manifest['slurm_job'])==job and manifest['config']==plan['configs'][mode] and manifest['build']==build
        for name,digest in manifest['source_sha256'].items():assert hashlib.sha256((repo/name).read_bytes()).hexdigest()==digest
        calls=[]
        for request in sorted(base.rglob('request.json')):
            case=request.parent;spec=core.load(request);exe=read(case/'execution.json');summary=read(case/'summary.json')
            rows=[];trace=case/'trajectory.tsv'
            if trace.exists():
                for line in trace.read_text().splitlines()[1:]:
                    fields=line.split('\t')
                    if len(fields)==7:rows.append([float(v) for v in fields])
            calls.append(dict(case=mode+'/'+case.relative_to(base).as_posix(),request=spec,execution=exe,summary=summary,
                observations=rows,active_elapsed_core_hours=(time.time()-request.stat().st_mtime)*spec['ranks']/3600 if not exe else 0.))
        state=read(base/'state.json');compact=None
        if state:
            if mode in ('dmc','qdmc'):compact=dict(stage=state['stage'],logz=state['logz'],statistics=state.get('statistics'))
            elif mode=='ams':compact=dict(round=state['round'],log_weight=state['log_weight'],events=sum(p['tau'] is not None for p in state['particles']))
            else:compact=dict(events=sum(p['tau'] is not None for p in state['particles']))
        runs[mode]=dict(manifest=manifest,status=read(base/'status.json'),state=compact,result=read(base/'result.json'),
                       selections=[core.load(p) for p in sorted(base.glob('selection_*.json'))],calls=calls)
    ids=','.join(str(v) for v in ledger['jobs'].values())
    return dict(job=ledger['jobs'],utc=datetime.datetime.utcnow().isoformat()+'Z',experiment=plan['experiment'],lifecycle=lifecycle,runs=runs,
                squeue=command(['squeue','-j',ids,'-h','-o','%i|%T|%M|%C|%R']),
                sacct=command(['sacct','-j',ids,'-X','-n','-P','--format=JobIDRaw,State,ElapsedRaw,AllocCPUS,ExitCode,Start,End']))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repo',required=True);a=p.parse_args();print(json.dumps(snapshot(a.repo),indent=2))
