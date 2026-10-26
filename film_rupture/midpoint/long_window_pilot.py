"""Fixed-horizon flat-film test; events, numerical failures and censoring stay separate."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import signal
import time
import numpy as np
from adaptive_runner import Trajectory,source_hashes,runtime


def atomic_json(path,data):
    temporary=path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(data,indent=2),encoding='utf-8')
    os.replace(str(temporary),str(path))


def observation(state):
    film=state.film
    power=abs(np.fft.fft2(state.h-state.h.mean())/film.n**2)**2
    k=np.fft.fftfreq(film.n)*film.n
    radius=np.sqrt(k[:,None]**2+k[None,:]**2)
    mask=(radius>24)&(radius<=46)
    expected=np.zeros_like(film.omega)
    nz=abs(film.omega)>1e-18
    d=film.phi*(1+film.slip)*film.lam/film.length**2
    expected[nz]=-d[nz]/film.omega[nz]*np.expm1(-2*film.omega[nz]*state.time)
    expected[~nz]=2*d[~nz]*state.time
    denominator=float(expected[mask].sum())
    return dict(time=state.time,steps=state.step_count,hmin=float(state.h.min()),hmax=float(state.h.max()),
        rms=float(np.std(state.h)),energy=film.energy(state.h),mean=float(state.h.mean()),
        volume=float(state.h.sum()*film.dx**2),
        high_band_power=float(power[mask].sum()),high_band_finite_time_linear_TCW=denominator,
        high_band_ratio=float(power[mask].sum()/denominator) if denominator else None,
        rejections=state.ledger['rejections'],max_residual=state.ledger['max_residual'])


def run_case(case,height,n,seed,cap,horizon,interval,wall_seconds,resume=False,stop_after_steps=None,trajectory_options=None):
    options=dict(trajectory_options or {})
    State=Trajectory
    if options:
        from accelerated_runner import Trajectory as State
    checkpoint=case/'state.npz'
    state=(State.load(checkpoint) if resume else State(height,n,seed,dt_cap=cap,**options))
    identity=state.metadata()
    protocol=dict(height=height,n=n,seed=seed,cap=cap,horizon=horizon,observation_interval=interval,
        tiered=True,event_resolution=1e-5,flat_initial_mean=1.,sources=identity['sources'],runtime=identity['runtime'],
        driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        interpretation='Fixed observation boundaries are part of the numerical protocol; not a MFPT estimate')
    if options:protocol['trajectory_options']=options
    trace=case/'observations.csv'
    if resume:
        assert json.loads((case/'protocol.json').read_text())==protocol,'Resume protocol changed'
        with trace.open(newline='',encoding='utf-8') as stream:
            rows=[dict(r) for r in csv.DictReader(stream) if float(r['time'])<=state.time+1e-12]
    else:
        case.mkdir(parents=True,exist_ok=False)
        atomic_json(case/'protocol.json',protocol)
        rows=[]
    end=state.noise.ticks(horizon);spacing=state.noise.ticks(interval)
    assert spacing>0 and end>0 and end%spacing==0
    stopped=[False]
    def stop(signum,frame):stopped[0]=True
    for sig in (signal.SIGTERM,signal.SIGINT):signal.signal(sig,stop)
    start=time.perf_counter();start_steps=state.step_count;start_time=state.time
    def observe():
        row=observation(state)
        if not rows or abs(float(rows[-1]['time'])-state.time)>1e-12:rows.append(row)
        with trace.open('w',newline='',encoding='utf-8') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(row));writer.writeheader();writer.writerows(rows)
        np.savez_compressed(case/('field_t%012d.npz'%state.ticks),h=state.h,
                            time=state.time,power=abs(np.fft.fft2(state.h-state.h.mean())/n**2)**2)
        print(json.dumps(dict(height=height,seed=seed,cap=cap,observation=row)),flush=True)
    def persist(outcome):
        state.save(checkpoint)
        identity=state.metadata()
        row=dict(height=height,n=n,seed=seed,cap=cap,outcome=outcome,time=state.time,
            event=state.event,steps=state.step_count,hmin=float(state.h.min()),finite=bool(np.isfinite(state.h).all()),
            independent_volume=float(state.h.sum()*state.film.dx**2),initial_volume=state.initial_mean*state.film.length**2,
            ledger=state.ledger,wall_seconds_this_segment=time.perf_counter()-start,
            steps_this_segment=state.step_count-start_steps,code_time_this_segment=state.time-start_time,
            sources=identity['sources'],runtime=identity['runtime'],kernel=identity.get('kernel','reference'),
            brownian_base_dt=state.noise.base_dt,brownian_depth=state.noise.depth,
            native_fingerprint=identity.get('fingerprint'),MFPT_validated=False)
        atomic_json(case/'summary.json',row)
        return row
    if not rows:observe()
    persist('RUNNING')
    while state.status=='RUNNING' and state.ticks<end:
        if stopped[0] or time.perf_counter()-start>=wall_seconds:break
        if stop_after_steps is not None and state.step_count-start_steps>=stop_after_steps:break
        boundary=min(end,(state.ticks//spacing+1)*spacing)
        state.advance(stop_ticks=boundary)
        if state.ticks==boundary or state.status!='RUNNING':observe();persist(state.status)
        elif state.step_count%1000==0:persist('RUNNING')
    if state.status!='RUNNING':outcome=state.status
    elif state.ticks==end:outcome='RIGHT_CENSORED_AT_FIXED_HORIZON'
    else:outcome='PAUSED_WITH_CHECKPOINT'
    if not rows or abs(float(rows[-1]['time'])-state.time)>1e-12:observe()
    result=persist(outcome)
    print('CASE_DONE',json.dumps(result),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--height',type=float,required=True)
    p.add_argument('--seed',type=int,required=True)
    p.add_argument('--caps',nargs='+',type=float,default=[.01,.03])
    p.add_argument('--n',type=int,default=64)
    p.add_argument('--horizon',type=float,default=300.)
    p.add_argument('--interval',type=float,default=3.)
    p.add_argument('--wall-seconds',type=float,default=6000.)
    p.add_argument('--resume',action='store_true')
    p.add_argument('--native-library',type=Path)
    p.add_argument('--brownian-base',type=float,default=.064)
    p.add_argument('--brownian-depth',type=int,default=18)
    args=p.parse_args()
    results=[]
    native=dict(kernel='native',library=str(args.native_library.resolve()),base_dt=args.brownian_base,depth=args.brownian_depth) if args.native_library else None
    for cap in args.caps:
        results.append(run_case(args.out/('cap%g'%cap),args.height,args.n,args.seed,cap,
            args.horizon,args.interval,args.wall_seconds,resume=args.resume,trajectory_options=native))
    atomic_json(args.out/'pair_summary.json',results)
    raise SystemExit(2 if any(r['outcome']=='NUMERICAL_FAILURE' for r in results) else
                     (3 if any(r['outcome']=='PAUSED_WITH_CHECKPOINT' for r in results) else 0))
