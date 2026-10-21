"""Actual new-variant serial/MPI/reference/restart/controller and q17 score gates.
No production campaign is started. Python3.6/stdlib; gate seeds are paired tests.
"""
from pathlib import Path
import argparse,hashlib,json,math,os,struct,subprocess
from rare_events import core,gates
H=struct.Struct('<8s64sIIQQ11d')
def direct_score(field):
    lo=min(field)
    if lo<=.01:return lo
    # Deliberately direct 2-D summation, independent of C's separable observer.
    box=min(math.fsum(field[((y+dy)%64)*64+(x+dx)%64] for dy in range(-8,9) for dx in range(-8,9))/289
            for y in range(64) for x in range(64))
    return lo+max(0.,min(1.,(lo-.01)/.19))*(max(box,lo)-lo)
def main(workdir,baseline):
    work=Path(workdir).resolve();baseline=Path(baseline).resolve();gates.main(work)
    build=core.build_identity();cfg=dict(N=4,ratio=1.01,slip_ref=.2,dt=.003,horizon=6.,ranks=4,workers=3,
                                       repulsion=False,score='q17',namespace='q17-numerical-gate-20261021')
    maximum=0.;score_cases=[]
    cases=[]
    for low in (.005,.009,.01,.010001,.19,.2,.200001,.7):
        field=[1.]*4096
        # Pit straddles both periodic boundaries.
        for y in (63,0,1):
            for x in (63,0,1):field[y*64+x]=low
        cases.append(('pit_'+str(low),field))
    cases += [('uniform',[1.]*4096),('wave',[1+.05*math.cos(2*math.pi*x/64)+.02*math.sin(2*math.pi*y/64) for y in range(64) for x in range(64)])]
    for label,field in cases:
        folder=work/'synthetic_inputs';folder.mkdir(exist_ok=True);path=folder/(label+'.chk')
        mass=math.fsum(field)*(100/1.01/64)**2
        header=H.pack(b'RESTFW17',build['engine_signature'].encode(),64,1,0,2026102101,.003,1.01,.2,.01/1.01**2,
                      .0013159/1.01**2,100/1.01,.01,.005,mass,1.,0.)
        path.write_bytes(header+struct.pack('<4096d',*field))
        out=core.execute(work/'score_static'/label,2026102101,0,cfg,restore=str(path),binary='serial',ranks=1)
        q=out['summary']['q17'];expected=direct_score(field);error=abs(q-expected);maximum=max(maximum,error)
        assert error<1e-12 and q>=min(field)-1e-14
        assert (q<=.01)==(min(field)<=.01)
        score_cases.append(dict(case=label,q17=q,expected=expected,error=error,rupture=out['summary']['reason']=='rupture'))
    production_fields=[];baseline_field_errors={}
    for seed in (2026101901,2026101902):
        p=work/('seed%d_mpi4'%seed)/'final.chk';data=core.checkpoint(p)
        expected=direct_score(data['field']);end=core.load(p.parent/'summary.json');error=abs(expected-end['q17']);maximum=max(maximum,error)
        assert error<1e-12;production_fields.append(dict(seed=seed,error=error))
        old=(baseline/('validation/gate_29241/seed%d_serial/final.chk'%seed)).read_bytes()
        assert H.unpack_from(old)[0]==b'RESTFE26'
        prior=struct.unpack_from('<4096d',old,H.size)
        now=core.checkpoint(work/('seed%d_serial'%seed)/'final.chk')['field']
        baseline_field_errors[str(seed)]=max(abs(x-y) for x,y in zip(prior,now))
        assert baseline_field_errors[str(seed)]<1e-12
    # Explicitly reject a same-physics q3 checkpoint. Score history is incompatible.
    old=baseline/'validation/gate_29241/seed2026101901_mpi4/final.chk';assert old.is_file()
    rejected=work/'reject_q3';rejected.mkdir()
    with (rejected/'stdout.log').open('w') as out,(rejected/'stderr.log').open('w') as err:
        code=subprocess.call([str(core.ROOT/'build/serial'),'--seed','1','--restore',str(old)],cwd=str(rejected),stdout=out,stderr=err)
    assert code!=0 and not (rejected/'final.chk').exists()
    # Physical reference C is identical to the already validated baseline reference.
    assert (core.ROOT/'physics.h').read_bytes()==(baseline/'rare_events/physics.h').read_bytes()
    assert (core.ROOT/'build/reference.c').read_bytes()==(baseline/'rare_events/build/reference.c').read_bytes()
    report=dict(status='PASS',job=int(os.environ['SLURM_JOB_ID']),build=build,score='q17',checkpoint_abi='RESTFW17',
      maximum_score_error=maximum,synthetic_cases=score_cases,actual_fields=production_fields,baseline_field_errors=baseline_field_errors,q3_checkpoint_rejected=True,
      scope='paired two-seed 2000-step physical-reference/true-serial/MPI1/4/16/repeat/restart/controller gate plus static score endpoints; not rupture/MFPT/near-event convergence')
    core.save(work/'score_gate.json',report);print(json.dumps(report,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workdir',required=True);p.add_argument('--baseline',required=True);a=p.parse_args();main(a.workdir,a.baseline)
