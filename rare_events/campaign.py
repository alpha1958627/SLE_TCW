"""Bounded first pilot, four independent groups sharing exactly the same physical kernel.

The pilot never claims an efficiency winner or a complete MFPT. Later experiments
must use fresh namespaces, fixed comparison windows and independent whole runs.
"""
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import argparse
import os
import traceback
from . import core
from .direct import run as direct
from rare_events.ams import run as ams
from rare_events.dmc import run as dmc


def main(root,gate):
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=True)
    verified=core.load(Path(gate)/'numerical_gate.json')
    assert verified['status']=='PASS' and verified['build']==core.build_identity()
    job=int(os.environ['SLURM_JOB_ID'])
    configs={}
    for mode in ('ams','dmc','qdmc','direct'):
        cfg=dict(mode=mode,N=12,ratio=1.01,slip_ref=.2,dt=.003,horizon=3000.,
            ranks=4,workers=3,score='q3',repulsion=False,
            namespace='rupture-fed7bb3-r101-no-repulsion-20261019-v1:job%d:%s:R0'%(job,mode))
        if mode=='ams':cfg.update(D=1,max_rounds=16)
        if mode in ('dmc','qdmc'):cfg.update(block=150.,strength=4. if mode=='dmc' else .5)
        configs[mode]=cfg
    seeds=[core.seed(c['namespace'],'initial:%d'%i) for c in configs.values() for i in range(c['N'])]
    assert len(set(seeds))==48
    metadata=dict(job=job,phase='runtime_gate',configs=configs,gate=str(Path(gate).resolve()),
                  independent_initial_seeds=48,physics='no_repulsion',
                  window_selection='initial fixed T3000, not a proven sufficient physical window')
    core.save(root/'campaign.json',metadata)
    from .gates import controller_smoke
    controller_smoke(root/'_runtime_gate',configs['ams'])
    metadata['phase']='running_bounded_pilot';core.save(root/'campaign.json',metadata)
    results={};errors={}
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(ams if m=='ams' else direct if m=='direct' else dmc,root/m,c):m for m,c in configs.items()}
        for future in as_completed(futures):
            mode=futures[future]
            try:results[mode]=future.result()
            except Exception:errors[mode]=traceback.format_exc()
            core.save(root/'campaign_status.json',dict(job=job,phase='running',results=results,errors=errors,
                whole_runs_complete=sum(r['status']=='complete' for r in results.values()),MFPT=None))
    core.save(root/'campaign_status.json',dict(job=job,phase='pilot_review_required',results=results,errors=errors,
        MFPT=None,efficiency_identified=False,
        next='inspect event evidence, field structure, families, cutoff and cost before new independent window/score tests'))
    if errors:raise RuntimeError('Failed groups retained: '+','.join(errors))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workdir',required=True);p.add_argument('--gate',required=True)
    a=p.parse_args();main(a.workdir,a.gate)
