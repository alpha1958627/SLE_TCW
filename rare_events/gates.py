"""Run physical and controller integration gates inside one >=16-core Slurm job."""
from pathlib import Path
import argparse
import math
import struct
import subprocess
import time
from . import core
from .test_estimators import self_check
from AMS.run import run as ams
from DMC.run import run as dmc
from .direct import run as direct


def controller_smoke(root,config):
    root=Path(root);results={}
    short=dict(config,N=4,workers=min(3,config['workers']),horizon=3.,mode='ams',D=1,max_rounds=2,
               namespace='new-engine-runtime-smoke-ams:'+root.name)
    ar=ams(root/'ams',short)
    assert ar['status'] in ('complete','budget_unresolved');results['ams']=ar
    for mode in ('dmc','qdmc'):
        cfg=dict(config,N=4,workers=min(3,config['workers']),horizon=3.,block=1.5,mode=mode,strength=2.,
                 namespace='new-engine-runtime-smoke-'+mode+':'+root.name)
        result=dmc(root/mode,cfg)
        assert result['status']=='complete' and result['horizon']==3. and result['MFPT'] is None
        results[mode]=result
    cfg=dict(config,N=4,workers=min(3,config['workers']),horizon=3.,mode='direct',
             namespace='new-engine-runtime-smoke-direct:'+root.name)
    result=direct(root/'direct',cfg)
    assert result['status']=='complete' and result['horizon']==3. and result['MFPT'] is None
    results['direct']=result
    report=dict(status='PASS',results=results,build=core.build_identity(),
                scope='actual current four-controller smoke; distinct from pilot statistics')
    core.save(root/'runtime_gate.json',report)
    return report


def main(root):
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=True)
    core.save(root/'estimator_gate.json',self_check())
    config=dict(N=4,ratio=1.01,slip_ref=.2,dt=.003,horizon=6.,ranks=4,workers=4,
                repulsion=False,score='q3',namespace='rupture-reference-gates-20261019')
    build=core.build_identity();errors={};seed1=2026101901
    def execute(name,seed=seed1,binary='mpi',ranks=4,steps=2000,restore=None,records=False,cfg=None):
        return core.execute(root/name,seed,steps,cfg or config,restore,records,binary,ranks)
    def values(out):return core.checkpoint(out['checkpoint'])['field']
    def difference(a,b):return max(abs(x-y) for x,y in zip(values(a),values(b)))
    for initial_seed in (seed1,seed1+1):
        prefix='seed%d_'%initial_seed
        reference=execute(prefix+'serial',initial_seed,'serial',1)
        # Full native reference uses integration i=0..1999: exactly 2000 kicks.
        native=root/(prefix+'native');native.mkdir(exist_ok=True)
        command=[str(core.ROOT/'build/reference_serial'),str(initial_seed),'.003','1999',str(.2/1.01),'1.','1.01']
        with (native/'stdout.log').open('w') as out,(native/'stderr.log').open('w') as err:
            code=subprocess.call(command,cwd=str(native),stdout=out,stderr=err)
        assert code==0
        raw=(native/'stop.hfield').read_bytes();assert struct.unpack_from('<i',raw)[0]==64
        field=struct.unpack_from('<4096d',raw,4)
        errors[prefix+'native_full_program']=max(abs(x-y) for x,y in zip(field,values(reference)))
        mpi4=None
        for rank in (1,4,16):
            out=execute(prefix+'mpi%d'%rank,initial_seed,ranks=rank)
            errors[prefix+'serial_mpi%d'%rank]=difference(reference,out)
            if rank==4:mpi4=out
        observed=execute(prefix+'observed',initial_seed,records=True)
        errors[prefix+'observer_noninterference']=difference(mpi4,observed)
        errors[prefix+'repeat4']=difference(mpi4,execute(prefix+'repeat',initial_seed))
        half=execute(prefix+'first_half',initial_seed,steps=1000)
        resumed=execute(prefix+'resumed',initial_seed,restore=half['checkpoint'])
        errors[prefix+'restart']=difference(mpi4,resumed)
        assert max(abs(out['summary']['mass_relative_drift']) for out in (reference,mpi4,resumed,observed))<1e-9
        branch=execute(prefix+'new_future',initial_seed+100,restore=half['checkpoint'])
        assert difference(branch,resumed)>1e-8
    assert max(errors.values())<1e-10,errors
    assert all(value==0 for key,value in errors.items() if 'repeat' in key or 'observer' in key)
    # A genuine half-dt execution: fixed .0015 must not ramp back to .003.
    fine=dict(config,dt=.0015)
    out=execute('fixed_half_dt',steps=2000,cfg=fine)
    assert out['summary']['time']==3. and core.checkpoint(out['checkpoint'])['dt']==.0015
    # Foreign checkpoints must fail before any physical update.
    foreign=root/'foreign.chk';foreign.write_bytes(b'AMS3'+b'\0'*128)
    rejected=root/'reject_foreign';rejected.mkdir(exist_ok=True)
    with (rejected/'stdout.log').open('w') as out,(rejected/'stderr.log').open('w') as err:
        code=subprocess.call([str(core.ROOT/'build/serial'),'--seed','1','--restore',str(foreign)],
                             cwd=str(rejected),stdout=out,stderr=err)
    assert code!=0 and not (rejected/'final.chk').exists()
    short=dict(config,horizon=3.,mode='ams',D=1,max_rounds=2,namespace='new-engine-controller-gate-20261019-ams')
    ar=ams(root/'controller_ams',short)
    assert ar['status'] in ('complete','budget_unresolved')
    for mode in ('dmc','qdmc'):
        cfg=dict(config,horizon=3.,block=1.5,mode=mode,strength=2.,namespace='new-engine-controller-gate-20261019-'+mode)
        result=dmc(root/('controller_'+mode),cfg)
        assert result['status']=='complete' and result['horizon']==3. and result['MFPT'] is None
    report=dict(status='PASS',build=build,errors=errors,half_dt_executed=True,foreign_checkpoint_rejected=True,
                reference='full film_rupture program, no repulsion, fixed dt, multigrid and Jacobi; 2000 updates',
                scope='two seeds 2000 steps, true serial/MPI1/4/16, observer/repeat/restart, real controller smoke; not rupture/MFPT/near-event/grid/noise acceptance')
    core.save(root/'numerical_gate.json',report)
    print(report)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--workdir',required=True)
    main(parser.parse_args().workdir)
