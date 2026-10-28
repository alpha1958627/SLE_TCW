"""Independent q17 candidate, actual controller startup first, fixed physical kernel."""
from pathlib import Path
import argparse,hashlib,json,os,sys
def main(root,mode):
    root=Path(root).resolve();source=root/'source';sys.path.insert(0,str(source))
    from rare_events import core
    from rare_events.ams import run as ams
    from rare_events.dmc import run as dmc
    plan=core.load(root/'plan.json');build=core.build_identity();gate=core.load(source/'validation/gate_30902/score_gate.json')
    assert gate['status']=='PASS' and gate['build']==build and max(gate['baseline_field_errors'].values())==0
    assert hashlib.sha256((source/'rare_events/physics.h').read_bytes()).hexdigest()=='1f108bce15a106991904c0d4744ef142b4fd6fc6bbf50e01796635aaa1d46290'
    assert core.load(root/'d3_oracle_sipf.json')['status']=='PASS'
    cfg=dict(plan['physics'],N=plan['N'],horizon=plan['horizon'],score='q17',mode=mode,
             namespace=plan['namespace_prefix']+':'+mode+':R0',**plan['modes'][mode])
    # Event/floor/grid are fixed inside the frozen engine, not alternate config knobs.
    for key in ('grid','event','floor'):cfg.pop(key)
    assert core.check_config(cfg)==4000000 and cfg['ratio']==1.01 and cfg['slip_ref']==.2 and cfg['repulsion'] is False
    job=int(os.environ['SLURM_JOB_ID']);work=root/'runs'/mode/('job_%d'%job);work.mkdir(parents=True)
    helper=root/'admin';manifest=dict(job=job,mode=mode,config=cfg,build=build,plan_sha256=hashlib.sha256((root/'plan.json').read_bytes()).hexdigest(),
        helpers={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in helper.iterdir() if p.is_file()},
        gate_job=30902,scope='independent q17/parameter feasibility test, not isolated score/framework causality')
    core.save(work/'experiment_manifest.json',manifest)
    function=ams if mode=='ams' else dmc
    smoke=dict(cfg,horizon=3.,namespace=cfg['namespace']+':runtime-job%d'%job)
    if mode=='ams':smoke.update(max_rounds=2) # Actual N12/D3, including multi-child parallel replacement.
    else:smoke.update(N=4,block=1.5)
    result=function(work/'startup_gate',smoke)
    assert result['status'] in ('complete','budget_unresolved') and result['MFPT'] is None
    core.save(work/'startup_gate.json',dict(status='PASS',result=result,
        scope='current actual MPI4/controller startup; AMS N12D3, DMC/QDMC actual chosen strength; not rupture validation'))
    result=function(work/'production',cfg);core.save(work/'result.json',result)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--mode',choices=('ams','dmc','qdmc'),required=True)
    a=p.parse_args();main(a.root,a.mode)
