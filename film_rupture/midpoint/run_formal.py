"""Formal direct MC, event-terminated; no artificial finite observation horizon.
Only this campaign's fingerprinted checkpoint may be resumed explicitly.
"""
import hashlib,json,os,signal,time
from pathlib import Path
import numpy as np
from basilisk_midpoint_candidate import Trajectory

ROOT=Path(__file__).resolve().parent
task=int(os.environ['SLURM_ARRAY_TASK_ID'])
campaign=json.loads((ROOT/'FORMAL_SELECTION.json').read_text())
assert campaign['PASS'] and json.loads((ROOT/'TCW_PASS.json').read_text())['PASS']
case=campaign['cases'][task];assert case['task']==task
folder=ROOT/'formal'/('task_%02d'%task)
resume=(folder/'state.npz').exists()
s=Trajectory.load(folder/'state.npz') if resume else Trajectory(case['height'],64,case['seed'],dt_cap=.3)
rank=s.film.rank
assert s.film.ranks==case['mpi_ranks']
stop=[False]
for sig in (signal.SIGTERM,signal.SIGINT):signal.signal(sig,lambda signum,frame:stop.__setitem__(0,True))
start=time.perf_counter();initial_time=s.time;initial_steps=s.step_count
if rank==0:
    folder.mkdir(parents=True,exist_ok=True)
    protocol=dict(case=case,cap=.3,n=64,hstar=0,ls_ref=.2,phi_ref=.01,initial_H=1,
        observation_interval=300.,event_resolution_nominal=1e-5,no_fixed_horizon=True,
        source_and_library_fingerprint=s.metadata()['fingerprint'],mpi_ranks=s.film.ranks,
        driver_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    protocol_path=folder/'protocol.json'
    if resume:assert json.loads(protocol_path.read_text())==protocol
    else:protocol_path.write_text(json.dumps(protocol,indent=2))
    trace=(folder/'observations.jsonl').open('a',buffering=1)

def atomic_json(path,data):
    temp=path.with_name(path.name+'.tmp');temp.write_text(json.dumps(data,indent=2));os.replace(str(temp),str(path))

def persist(outcome,field=False):
    if rank!=0:return
    s.save(folder/'state.npz')
    summary=dict(case=case,outcome=outcome,time=s.time,steps=s.step_count,event=s.event,
        hmin=float(s.h.min()),hmax=float(s.h.max()),mean=float(s.h.mean()),
        independent_volume=float(s.h.sum()*s.film.dx**2),initial_volume=s.film.length**2,
        finite=bool(np.isfinite(s.h).all()),ledger=s.ledger,mpi_ranks=s.film.ranks,
        node=os.environ.get('SLURMD_NODENAME'),job=os.environ.get('SLURM_JOB_ID'),
        wall_seconds_this_segment=time.perf_counter()-start,
        code_time_this_segment=s.time-initial_time,steps_this_segment=s.step_count-initial_steps,
        fingerprint=s.metadata()['fingerprint'],no_fixed_horizon=True,MFPT_validated=False)
    if field:
        power=abs(np.fft.fft2(s.h-s.h.mean())/64**2)**2
        name='field_t%016d.npz'%s.ticks
        np.savez_compressed(folder/name,h=s.h,time=s.time,power=power)
        summary['last_field']=name
        trace.write(json.dumps(dict(time=s.time,steps=s.step_count,hmin=summary['hmin'],hmax=summary['hmax'],
            mean=summary['mean'],volume=summary['independent_volume'],ledger=s.ledger,field=name))+'\n')
        print('OBSERVATION',json.dumps(summary),flush=True)
    atomic_json(folder/'summary.json',summary)

if s.status=='RUPTURE_BRACKET':
    persist(s.status,True)
else:
    persist('RUNNING',True)
    spacing=s.noise.ticks(300.)
    while s.status=='RUNNING':
        boundary=(s.ticks//spacing+1)*spacing
        s.advance(stop_ticks=boundary)
        if s.ticks==boundary or s.status!='RUNNING':persist(s.status,True)
        elif s.step_count%1000==0:persist('RUNNING')
        # All ranks consult the stop flag at the same accepted-step boundary.
        if s.step_count%1000==0 and s.film.library.bm_any(int(stop[0])):break
    outcome=s.status if s.status!='RUNNING' else 'PAUSED_WITH_CHECKPOINT'
    persist(outcome,True)
    if rank==0:print('CASE_DONE',outcome,flush=True)
    if outcome=='NUMERICAL_FAILURE':raise SystemExit(2)
    if outcome=='PAUSED_WITH_CHECKPOINT':raise SystemExit(3)
