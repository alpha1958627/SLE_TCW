"""One independent finite-window method per Slurm job. Frozen physical engine unchanged."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import sys

def main(repo,plan_path,mode):
    repo=Path(repo).resolve();plan_path=Path(plan_path).resolve()
    sys.path.insert(0,str(repo))
    from rare_events import core
    from rare_events.test_estimators import self_check
    from rare_events.ams import run as ams
    from rare_events.dmc import run as dmc
    from rare_events.direct import run as direct
    plan=core.load(plan_path);cfg=plan['configs'][mode]
    build=core.build_identity();gate=core.load(repo/plan['prior_gate']/'numerical_gate.json')
    assert gate['status']=='PASS' and gate['build']==build and build['engine_signature']==plan['expected_engine_signature']
    prior=repo/'runs'/('job_%d'%plan['prior_job']);status=core.load(prior/'campaign_status.json')
    assert status['phase']=='pilot_review_required' and not status['errors']
    assert status['results']['ams']['status']=='budget_unresolved'
    assert all(status['results'][m]['status']=='complete' and status['results'][m]['P']==0 for m in ('dmc','qdmc','direct'))
    for m in plan['configs']:
        manifest=core.load(prior/m/'manifest.json')
        assert manifest['build']==build
        for name,digest in manifest['source_sha256'].items():
            assert hashlib.sha256((repo/name).read_bytes()).hexdigest()==digest,name
    old=[core.seed(core.load(prior/m/'manifest.json')['config']['namespace'],'initial:%d'%i) for m in plan['configs'] for i in range(12)]
    new=[core.seed(c['namespace'],'initial:%d'%i) for c in plan['configs'].values() for i in range(12)]
    assert len(set(new))==48 and not set(new).intersection(old)
    assert core.check_config(cfg)==4000000 and cfg['score']=='q3' and cfg['dt']==.003 and cfg['ratio']==1.01 and cfg['slip_ref']==.2
    job=int(os.environ['SLURM_JOB_ID']);base=repo/'experiments'/plan['experiment']/mode
    base.mkdir(parents=True,exist_ok=True)
    claim=dict(job=job,mode=mode,namespace=cfg['namespace'],plan_sha256=hashlib.sha256(plan_path.read_bytes()).hexdigest())
    # Claim each scientific repetition once. An uncertain attempt is preserved, never replayed automatically.
    with (base/'claim.json').open('x') as handle:json.dump(claim,handle,indent=2)
    work=base/('job_%d'%job);work.mkdir()
    core.save(work/'experiment_manifest.json',dict(claim,config=cfg,build=build,
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),prior_job=plan['prior_job']))
    core.save(work/'estimator_gate.json',self_check())
    function=ams if mode=='ams' else direct if mode=='direct' else dmc
    smoke=dict(cfg,N=4,workers=3,horizon=3.,namespace=cfg['namespace']+':runtime-job%d'%job)
    if mode=='ams':smoke.update(max_rounds=1)
    elif mode in ('dmc','qdmc'):smoke.update(block=1.5)
    checked=function(work/'startup_gate',smoke)
    assert checked['status'] in ('complete','budget_unresolved') and checked['MFPT'] is None
    core.save(work/'startup_gate.json',dict(status='PASS',result=checked,scope='same-mode current job actual MPI4/controller startup only'))
    result=function(work/'production',cfg)
    core.save(work/'result.json',result)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repo',required=True);p.add_argument('--plan',required=True);p.add_argument('--mode',choices=('ams','dmc','qdmc','direct'),required=True)
    a=p.parse_args();main(a.repo,a.plan,a.mode)
