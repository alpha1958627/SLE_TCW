"""Claim and submit each new candidate once. Test-only IDs are never actual jobs."""
from pathlib import Path
import hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
def run(command):
    p=subprocess.run(command,cwd=str(ROOT),stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    return dict(command=command,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr)
def main():
    sys.path.insert(0,str(ROOT/'source'));from rare_events import core
    sys.path.insert(0,str(ROOT/'admin'));from check_d3_oracle import check
    gate=core.load(ROOT/'source/validation/gate_30902/score_gate.json');assert gate['status']=='PASS'
    assert gate['build']==core.build_identity() and max(gate['baseline_field_errors'].values())==0
    plan=core.load(ROOT/'plan.json');assert plan['status']=='prepared_not_submitted'
    baseline=ROOT.parents[1];seen=[]
    # Deliberately separate production namespaces from both prior independent ensembles.
    for mode in ('ams','dmc','qdmc','direct'):
        for ns in ('rupture-q3-T12000-20261020-v1:'+mode+':R0','rupture-fed7bb3-r101-no-repulsion-20261019-v1:job29264:'+mode+':R0'):
            seen += [core.seed(ns,'initial:%d'%i) for i in range(12)]
    fresh=[core.seed(plan['namespace_prefix']+':'+m+':R0','initial:%d'%i) for m in plan['modes'] for i in range(12)]
    assert len(set(fresh))==36 and not set(fresh).intersection(seen)
    core.save(ROOT/'d3_oracle_sipf.json',check())
    # Check the completed baseline finite window; leave the baseline AMS controller untouched.
    for mode,job in (('dmc',29485),('qdmc',29486),('direct',29487)):
        r=core.load(baseline/'experiments/window_q3_T12000_20261020'/mode/('job_%d'%job)/'production/result.json')
        assert r['status']=='complete' and r['horizon']==12000 and r['P']==0 and r['MFPT'] is None
    with (ROOT/'production_submit_claim.json').open('x') as f:json.dump(dict(plan_sha256=hashlib.sha256((ROOT/'plan.json').read_bytes()).hexdigest()),f)
    ledger=dict(status='submitting',jobs={},attempts=[],scope='independent q17 candidate; baseline q3 unchanged')
    def save():(ROOT/'submission_production.json').write_text(json.dumps(ledger,indent=2)+'\n')
    save()
    for mode in plan['modes']:
        pre=run(['sbatch','--test-only','--job-name=q17-'+mode,'--export=ALL,MODE='+mode,'admin/candidate.slurm'])
        ledger['attempts'].append(dict(mode=mode,test_only=True,result=pre));save();assert pre['returncode']==0,pre
        actual=run(['sbatch','--parsable','--job-name=q17-'+mode,'--export=ALL,MODE='+mode,'admin/candidate.slurm'])
        ledger['attempts'].append(dict(mode=mode,test_only=False,result=actual));save();assert actual['returncode']==0,actual
        ledger['jobs'][mode]=int(actual['stdout'].strip().split(';')[0]);save()
    ledger['status']='submitted';ledger['initial_squeue']=run(['squeue','-j',','.join(str(j) for j in ledger['jobs'].values()),'-h','-o','%i|%T|%M|%C|%R']);save()
    print(json.dumps(ledger,indent=2))
if __name__=='__main__':main()
