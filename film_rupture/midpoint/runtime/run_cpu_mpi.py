"""Direct-MC midpoint solver: one independent trajectory per MPI rank.

Uniform periodic FV, FP64, same-path rejection, strict checkpoint fingerprints.
This is ensemble parallelism; it does not spatially decompose a single film.
"""
import os
for _name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
import argparse,csv,datetime,fcntl,hashlib,json,math,sys,time,traceback
from pathlib import Path
import numpy as np
from mpi4py import MPI
from accelerated_runner import Trajectory,fingerprint
from native_midpoint_candidate import build

ROOT=Path(__file__).resolve().parent
COMM=MPI.COMM_WORLD

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def atomic_json(p,value):
    temp=p.with_name(p.name+'.tmp');temp.write_text(json.dumps(value,indent=2,sort_keys=True));os.replace(str(temp),str(p))

def raw_field(p,s,protocol_sha):
    power=abs(np.fft.fft2(s.h-s.h.mean())/s.film.n**2)**2
    temp=p.with_name(p.name+'.tmp')
    with temp.open('wb') as f:
        np.savez_compressed(f,h=s.h,power=power,t=s.time,ticks=s.ticks,protocol_SHA256=protocol_sha)
    os.replace(str(temp),str(p))
    return sha(p)

def identity(s,case,args):
    assert s.film.height==case['height'] and s.noise.seed==case['seed'],'Checkpoint case identity mismatch'
    assert s.film.n==args.n and s.cap_ticks==s.noise.ticks(args.dt_max),'Checkpoint grid/dt mismatch'
    assert s.noise.base_dt==.064 and s.noise.depth==18 and s.tiered,'Checkpoint protocol mismatch'

def solve(case,args,library,protocol_sha):
    folder=args.out/case['id'];folder.mkdir(exist_ok=True)
    checkpoint=folder/'checkpoint.npz';trace=folder/'trace.csv'
    s=None;start=time.perf_counter();first_steps=0
    try:
        if checkpoint.exists():
            s=Trajectory.load(checkpoint,library=library)
        else:
            s=Trajectory(case['height'],args.n,case['seed'],dt_cap=args.dt_max,base_dt=.064,depth=18,
                         kernel=args.kernel,library=library)
        identity(s,case,args)
        first_steps=s.step_count
        # A killed process may have appended trace rows after its last checkpoint.
        if trace.exists():
            with trace.open() as f:old=list(csv.DictReader(f))
            keep=[r for r in old if int(r['ticks'])<=s.ticks]
            with trace.open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=['ticks','t','step','hmin','hmax','mean','dt','residual']);w.writeheader();w.writerows(keep)
        last_logged=None
        if trace.exists() and trace.stat().st_size:
            with trace.open() as f:rr=list(csv.DictReader(f))
            if rr:last_logged=int(rr[-1]['ticks'])
        def log():
            nonlocal last_logged
            if last_logged==s.ticks:return
            exists=trace.exists() and trace.stat().st_size>0
            with trace.open('a',newline='') as f:
                w=csv.writer(f)
                if not exists:w.writerow(['ticks','t','step','hmin','hmax','mean','dt','residual'])
                dt=0 if s.last_interval is None else s.last_interval[1]-s.last_interval[0]
                w.writerow([s.ticks,s.time,s.step_count,float(s.h.min()),float(s.h.max()),float(s.h.mean()),dt,s.ledger['max_residual']])
            last_logged=s.ticks
        def persist(state):
            assert np.isfinite(s.h).all() and s.h.min()>0,'Invalid accepted field'
            assert abs(float(s.h.mean())/s.initial_mean-1)<1e-10,'Independent mass check failed'
            assert s.ledger['floor_cells']==0 and s.ledger['max_residual']<=1e-10
            log();s.save(checkpoint)
            field_sha=raw_field(folder/'last_field.npz',s,protocol_sha)
            result=dict(case=case,state=state,at=now(),rank=COMM.rank,node=MPI.Get_processor_name(),
                resumed_from_checkpoint=bool(first_steps),
                CPU_affinity=sorted(os.sched_getaffinity(0)),t=s.time,ticks=s.ticks,steps=s.step_count,
                steps_this_run=s.step_count-first_steps,hmin=float(s.h.min()),hmax=float(s.h.max()),
                mean=float(s.h.mean()),independent_volume=float(s.h.sum()*s.film.dx**2),
                relative_mass_error=abs(float(s.h.mean())/s.initial_mean-1),event=s.event,ledger=s.ledger,
                wall_seconds_this_run=time.perf_counter()-start,protocol_SHA256=protocol_sha,
                checkpoint_SHA256=sha(checkpoint),field_SHA256=field_sha,fingerprint=s._fingerprint)
            atomic_json(folder/'summary.json',result)
            return result
        target=s.noise.ticks(args.until) if args.until>0 else None
        interval=s.noise.ticks(args.frame_interval)
        assert interval>0 and abs(interval*s.noise.quantum-args.frame_interval)<1e-12
        if target is not None:assert target>=s.ticks,'Target precedes checkpoint'
        log()
        while s.status=='RUNNING' and (target is None or s.ticks<target):
            if time.perf_counter()-start>=args.wall_seconds or (args.steps_this_run and s.step_count-first_steps>=args.steps_this_run):
                return persist('PAUSED_WITH_CHECKPOINT')
            next_frame=(s.ticks//interval+1)*interval
            stop=next_frame if target is None else min(next_frame,target)
            s.advance(stop_ticks=stop)
            if s.step_count%args.trace_every==0:log()
            if s.ticks==next_frame:
                raw_field(folder/('frame_%020d.npz'%s.ticks),s,protocol_sha)
            if s.step_count%args.checkpoint_every==0 or s.ticks==next_frame:persist(s.status)
        state=s.status if s.status!='RUNNING' else 'CENSORED_AT_TARGET'
        return persist(state)
    except Exception:
        # I/O/provenance errors must not be mislabeled as physical rupture or
        # exhausted nonlinear solves. Preserve the existing checkpoint as is.
        result=dict(case=case,state='EXECUTION_ERROR',at=now(),rank=COMM.rank,error=traceback.format_exc())
        atomic_json(folder/'execution_error.json',result)
        return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--heights',type=float,nargs='+',default=[1.2])
    p.add_argument('--seeds',type=int,nargs='+',default=[2026102201])
    p.add_argument('--n',type=int,default=64);p.add_argument('--dt-max',type=float,default=.3)
    p.add_argument('--until',type=float,default=300.,help='Fixed code-time target; 0 means event-only')
    p.add_argument('--wall-seconds',type=float,default=3600.)
    p.add_argument('--frame-interval',type=float,default=300.)
    p.add_argument('--checkpoint-every',type=int,default=1000)
    p.add_argument('--trace-every',type=int,default=100)
    p.add_argument('--steps-this-run',type=int,default=0,help='Pause after this many accepted steps; 0 disables')
    p.add_argument('--kernel',choices=['native','reference'],default='native')
    p.add_argument('--resume',action='store_true',help='Explicitly resume this output protocol')
    p.add_argument('--cpu-base',type=int,help='Optional Linux CPU id of rank0, rank r uses base+r')
    a=p.parse_args();a.out=a.out.resolve()
    assert a.n>=4 and a.n%2==0 and math.isfinite(a.dt_max) and a.dt_max>0
    assert math.isfinite(a.until) and a.until>=0 and math.isfinite(a.wall_seconds) and a.wall_seconds>0
    assert math.isfinite(a.frame_interval) and a.frame_interval>0
    assert a.checkpoint_every>0 and a.trace_every>0 and a.steps_this_run>=0
    assert len(set(a.heights))==len(a.heights) and all(math.isfinite(h) and h>0 for h in a.heights)
    assert len(set(a.seeds))==len(a.seeds) and all(0<=s<2**64 for s in a.seeds)
    if a.cpu_base is not None:os.sched_setaffinity(0,{a.cpu_base+COMM.rank})
    cases=[dict(id='case_%04d'%i,height=h,seed=s) for i,(h,s) in enumerate((h,s) for h in a.heights for s in a.seeds)]
    guard=None;setup=None
    if COMM.rank==0:
        try:
            a.out.mkdir(parents=True,exist_ok=True)
            guard=(a.out/'.run.lock').open('a');fcntl.flock(guard,fcntl.LOCK_EX|fcntl.LOCK_NB)
            if a.resume:assert (a.out/'protocol.json').exists(),'No identified protocol to resume'
            else:assert not (a.out/'protocol.json').exists(),'Output exists: use --resume explicitly'
            library=None
            if a.kernel=='native':
                compiled=ROOT/'compiled'
                if not compiled.exists():build(compiled)
                library=str(compiled/'kernel.so')
            core=fingerprint(a.kernel,library)
            protocol=dict(cases=cases,n=a.n,dt_max=a.dt_max,until=a.until,frame_interval=a.frame_interval,
                checkpoint_every=a.checkpoint_every,trace_every=a.trace_every,kernel=a.kernel,
                fingerprint=core,driver_SHA256=sha(__file__),base_dt=.064,depth=18,
                scheme='Conservative stochastic midpoint; direct MC; ensemble MPI')
            dest=a.out/'protocol.json'
            if dest.exists():assert json.loads(dest.read_text())==protocol,'Existing output protocol differs'
            else:
                assert not list(a.out.glob('case_*')),'Refuse to adopt unidentified output cases'
                atomic_json(dest,protocol)
            setup=dict(ok=True,library=library,protocol_SHA256=sha(dest))
        except Exception:setup=dict(ok=False,error=traceback.format_exc())
    setup=COMM.bcast(setup,root=0)
    if not setup['ok']:
        if COMM.rank==0:print(setup['error'],file=sys.stderr,flush=True)
        return 2
    COMM.Barrier();start=MPI.Wtime()
    results=[solve(case,a,setup['library'],setup['protocol_SHA256']) for i,case in enumerate(cases) if i%COMM.size==COMM.rank]
    gathered=COMM.gather(results,root=0);elapsed=COMM.reduce(MPI.Wtime()-start,op=MPI.MAX,root=0)
    code=None
    if COMM.rank==0:
        results=sorted(sum(gathered,[]),key=lambda r:r['case']['id']);assert len(results)==len(cases)
        counts={state:sum(r['state']==state for r in results) for state in sorted({r['state'] for r in results})}
        summary=dict(at=now(),MPI_ranks=COMM.size,identities=len(cases),counts=counts,cases=results,
            batch_wall_seconds=elapsed,protocol_SHA256=setup['protocol_SHA256'],MPI=MPI.Get_library_version(),
            implementation_scope='Independent trajectories; no spatial MPI decomposition')
        atomic_json(a.out/'campaign_summary.json',summary)
        print(json.dumps(dict(out=str(a.out),ranks=COMM.size,counts=counts,batch_wall_seconds=elapsed)),flush=True)
        code=1 if any(r['state'] in ('NUMERICAL_FAILURE','EXECUTION_ERROR') for r in results) else 0
    code=COMM.bcast(code,root=0)
    if guard is not None:guard.close()
    return code

if __name__=='__main__':sys.exit(main())
