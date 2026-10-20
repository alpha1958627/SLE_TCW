"""Preflight and submit each distinct scientific method exactly once, recording uncertain attempts."""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys

def main(repo):
    repo=Path(repo).resolve();folder=Path(__file__).resolve().parent
    sys.path.insert(0,str(repo));from rare_events import core
    plan=core.load(folder/'plan.json');build=core.build_identity()
    assert build['engine_signature']==plan['expected_engine_signature']
    prior=core.load(repo/'runs'/('job_%d'%plan['prior_job'])/'campaign_status.json')
    assert prior['phase']=='pilot_review_required' and not prior['errors']
    assert prior['results']['ams']['status']=='budget_unresolved'
    seeds=[core.seed(c['namespace'],'initial:%d'%i) for c in plan['configs'].values() for i in range(12)]
    assert len(set(seeds))==48
    old=[core.seed(core.load(repo/'runs'/('job_%d'%plan['prior_job'])/m/'manifest.json')['config']['namespace'],'initial:%d'%i) for m in plan['configs'] for i in range(12)]
    assert not set(seeds).intersection(old)
    for c in plan['configs'].values():assert core.check_config(c)==4000000
    ledger_path=folder/'submission.json'
    with ledger_path.open('x') as handle:json.dump(dict(status='preflight',jobs={},attempts=[],unique_new_initial_seeds=48),handle)
    ledger=core.load(ledger_path)
    common=['sbatch','--parsable','--output',str(folder/'%x-%j.out'),'--error',str(folder/'%x-%j.err')]
    for mode in ('ams','dmc','qdmc','direct'):
        command=common+['--job-name','q3W12k-'+mode,'--export','ALL,MODE='+mode,str(folder/'window_20261020.slurm')]
        test=subprocess.run(command[:-1]+['--test-only',command[-1]],cwd=str(repo),stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
        ledger['attempts'].append(dict(mode=mode,test_only=True,returncode=test.returncode,stdout=test.stdout,stderr=test.stderr));core.save(ledger_path,ledger)
        assert test.returncode==0,'Preflight failed; no automatic fallback or retry here'
    ledger['status']='submitting';core.save(ledger_path,ledger)
    for mode in ('ams','dmc','qdmc','direct'):
        command=common+['--job-name','q3W12k-'+mode,'--export','ALL,MODE='+mode,str(folder/'window_20261020.slurm')]
        ledger['uncertain_mode']=mode;core.save(ledger_path,ledger)
        out=subprocess.run(command,cwd=str(repo),stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
        ledger['attempts'].append(dict(mode=mode,test_only=False,returncode=out.returncode,stdout=out.stdout,stderr=out.stderr));core.save(ledger_path,ledger)
        assert out.returncode==0,'Unknown submission state; inspect ledger/sacct, never replay blindly'
        job=int(out.stdout.strip().split(';')[0]);ledger['jobs'][mode]=job;ledger.pop('uncertain_mode');core.save(ledger_path,ledger)
    ledger['status']='submitted';core.save(ledger_path,ledger);print(json.dumps(ledger))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repo',required=True);a=p.parse_args();main(a.repo)
