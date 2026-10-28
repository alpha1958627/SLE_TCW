"""One source-checked snapshot of the exact q17 candidate jobs, no solver actions."""
from pathlib import Path
import datetime,hashlib,json,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'source'))
from rare_events import core
def command(args):
    p=subprocess.run(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    return dict(returncode=p.returncode,stdout=p.stdout.strip(),stderr=p.stderr.strip())
def main():
    ledger=core.load(ROOT/'submission_production.json');assert ledger['status']=='submitted';build=core.build_identity();runs={};life={}
    for mode,job in ledger['jobs'].items():
        work=ROOT/'runs'/mode/('job_%d'%job);base=work/'production'
        def read(p):return core.load(p) if p.exists() else None
        life[mode]=dict(job=job,workdir=str(work),experiment_manifest=read(work/'experiment_manifest.json'),startup_gate=read(work/'startup_gate.json'))
        manifest=read(base/'manifest.json')
        if not manifest:continue
        assert int(manifest['slurm_job'])==job and manifest['build']==build
        for name,digest in manifest['source_sha256'].items():assert hashlib.sha256((ROOT/'source'/name).read_bytes()).hexdigest()==digest
        calls=[]
        for p in sorted(base.rglob('request.json')):
            spec=core.load(p);exe=read(p.parent/'execution.json');trace=p.parent/'trajectory.tsv';rows=[]
            if trace.exists():
                for line in trace.read_text().splitlines()[1:]:
                    v=line.split('\t')
                    if len(v)==7:rows.append([float(x) for x in v])
            calls.append(dict(case=mode+'/'+p.parent.relative_to(base).as_posix(),request=spec,execution=exe,summary=read(p.parent/'summary.json'),observations=rows,
                              active_elapsed_core_hours=(time.time()-p.stat().st_mtime)*spec['ranks']/3600 if not exe else 0.))
        state=read(base/'state.json');compact=None
        if state:
            compact=dict(round=state['round'],log_weight=state['log_weight'],events=sum(p['tau'] is not None for p in state['particles'])) if mode=='ams' else dict(stage=state['stage'],logz=state['logz'],statistics=state.get('statistics'))
        runs[mode]=dict(manifest=manifest,status=read(base/'status.json'),state=compact,result=read(base/'result.json'),selections=[core.load(p) for p in sorted(base.glob('selection_*.json'))],calls=calls)
    ids=','.join(str(j) for j in ledger['jobs'].values())
    return dict(job=ledger['jobs'],utc=datetime.datetime.utcnow().isoformat()+'Z',experiment='width_q17_20261021',lifecycle=life,runs=runs,
       squeue=command(['squeue','-j',ids,'-h','-o','%i|%T|%M|%C|%R']),sacct=command(['sacct','-j',ids,'-X','-n','-P','--format=JobIDRaw,State,ElapsedRaw,AllocCPUS,ExitCode,Start,End']))
if __name__=='__main__':print(json.dumps(main(),indent=2))
