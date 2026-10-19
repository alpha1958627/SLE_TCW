"""Shared CLI and immutable run manifest. Local runs require an explicitly new namespace."""
from pathlib import Path
import argparse
import hashlib
import os
from . import core


def arguments(mode):
    p=argparse.ArgumentParser(description=mode.upper()+' finite-window rupture sampling')
    p.add_argument('--workdir',required=True)
    p.add_argument('--namespace',required=True)
    p.add_argument('--N',type=int,default=12)
    p.add_argument('--ratio',type=float,default=1.01)
    p.add_argument('--slip',type=float,default=.2)
    p.add_argument('--dt',type=float,default=.003)
    p.add_argument('--horizon',type=float,default=12000.)
    p.add_argument('--ranks',type=int,choices=(1,4,16),default=4)
    p.add_argument('--workers',type=int,default=4)
    p.add_argument('--resume',action='store_true')
    if mode=='ams':
        p.add_argument('--D',type=int,default=1)
        p.add_argument('--max-rounds',type=int,default=32)
    elif mode=='dmc':
        p.add_argument('--mode',choices=('dmc','qdmc'),default='dmc')
        p.add_argument('--block',type=float,default=150.)
        p.add_argument('--strength',type=float,default=None)
    a=p.parse_args()
    config=dict(N=a.N,ratio=a.ratio,slip_ref=a.slip,dt=a.dt,horizon=a.horizon,
                ranks=a.ranks,workers=a.workers,repulsion=False,score='q3',namespace=a.namespace)
    if mode=='ams':config.update(mode='ams',D=a.D,max_rounds=a.max_rounds)
    elif mode=='dmc':config.update(mode=a.mode,block=a.block,strength=a.strength if a.strength is not None else (6. if a.mode=='dmc' else 1.))
    else:config.update(mode='direct')
    core.check_config(config)
    assert 1<=config['workers']<=config['N']
    return Path(a.workdir).resolve(),config,a.resume


def begin(workdir,config,resume):
    import fcntl
    workdir.mkdir(parents=True,exist_ok=True)
    lock=(workdir/'controller.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    files=list(core.ROOT.glob('*.py'))+list((core.ROOT.parent/'AMS').glob('*.py'))+list((core.ROOT.parent/'DMC').glob('*.py'))
    manifest=dict(config=config,build=core.build_identity(),
        source_sha256={p.relative_to(core.ROOT.parent).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        slurm_job=os.environ.get('SLURM_JOB_ID'),model='no_repulsion',checkpoint_abi='RESTFE26',
        time_reference_factor=config['ratio'],old_states_allowed=False)
    path=workdir/'manifest.json'
    if path.exists():
        assert resume,'Existing run requires explicit --resume after confirming the old process is stopped'
        old=core.load(path)
        for key in ('config','build','source_sha256'):assert old[key]==manifest[key], 'Run identity changed: '+key
    else:core.save(path,manifest)
    return lock
