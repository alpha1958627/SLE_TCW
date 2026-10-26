"""Same-path implementation gates and same-node allocation benchmarks."""
import argparse,hashlib,json,os,time
from pathlib import Path
import numpy as np
from basilisk_midpoint_candidate import Trajectory as BasiliskTrajectory
from accelerated_runner import Trajectory as NativeTrajectory
from midpoint_candidate import Film as ReferenceFilm

ROOT=Path(__file__).resolve().parent
HEIGHTS=[.6,.7,.8,.9,1.,1.1,1.2,1.3]

def save(path,data):
    path.parent.mkdir(exist_ok=True,parents=True)
    path.write_text(json.dumps(data,indent=2))

def gate():
    results=[]
    for height in HEIGHTS:
        for cap in [.03,.3]:
            b=BasiliskTrajectory(height,64,2026102201,dt_cap=cap)
            n=NativeTrajectory(height,64,2026102201,dt_cap=cap,base_dt=.064,depth=18,library=ROOT/'compiled/kernel.so')
            maximum=0.;restart=None
            for step in range(60):
                b.advance();n.advance()
                maximum=max(maximum,float(np.max(abs(b.h-n.h))))
                assert b.status==n.status=='RUNNING' and b.ticks==n.ticks
                if step==20:
                    if b.film.rank==0:b.save(ROOT/'admission/restart.npz')
                    b.film.library.bm_sync()
                    restart=BasiliskTrajectory.load(ROOT/'admission/restart.npz')
                elif step>20:
                    restart.advance();assert np.array_equal(restart.h,b.h),'Restart is not bitwise reproducible'
            assert maximum<1e-8 and b.ledger['max_residual']<=1e-10 and abs(b.h.mean()-1)<1e-12
            row=dict(height=height,cap=cap,max_field_difference=maximum,restart_bitwise=True,
                residual=b.ledger['max_residual'],volume_relative_error=abs(float(b.h.mean())-1),mpi_ranks=b.film.ranks)
            results.append(row)
            if b.film.rank==0:print('GATE',json.dumps(row),flush=True)
    if b.film.rank==0:save(ROOT/'ADMISSION_PASS.json',dict(PASS=True,rows=results))

def mpi_gate():
    b=BasiliskTrajectory(1.05,64,2026102201,dt_cap=.3)
    start=time.perf_counter()
    for step in range(100):b.advance()
    elapsed=time.perf_counter()-start
    if b.film.rank==0:
        folder=ROOT/'mpi_gates';folder.mkdir(exist_ok=True)
        np.savez_compressed(folder/('ranks%d.npz'%b.film.ranks),h=b.h,time=b.time)
        save(folder/('ranks%d.json'%b.film.ranks),dict(ranks=b.film.ranks,time=b.time,wall=elapsed,
            steps=b.step_count,hmin=float(b.h.min()),mean=float(b.h.mean()),max_residual=b.ledger['max_residual']))

def benchmark(kind,repeat,height=1.05):
    if kind=='basilisk':s=BasiliskTrajectory(height,64,2026102201,dt_cap=.3)
    else:s=NativeTrajectory(height,64,2026102201,dt_cap=.3,base_dt=.064,depth=18,library=ROOT/'compiled/kernel.so')
    rank=s.film.rank if kind=='basilisk' else int(os.environ.get('SLURM_PROCID','0'))
    ranks=s.film.ranks if kind=='basilisk' else 1
    # Identical warmup and a fixed code-time window, no asynchronous field output.
    warm=s.noise.ticks(3.)
    while s.ticks<warm:s.advance(stop_ticks=warm)
    start=time.perf_counter();start_steps=s.step_count;end=s.noise.ticks(33.)
    while s.ticks<end:s.advance(stop_ticks=end)
    elapsed=time.perf_counter()-start
    if rank==0:
        result=dict(solver=kind,mpi_ranks=ranks,repeat=repeat,node=os.environ.get('SLURMD_NODENAME'),
            allocation=os.environ.get('SLURM_JOB_ID'),time_start=3.,time_end=33.,time_advanced=30.,
            steps=s.step_count-start_steps,wall_seconds=elapsed,speed=30./elapsed,
            height=height,n=64,cap=.3,seed=2026102201,max_residual=s.ledger['max_residual'],hmin=float(s.h.min()),
            CPU_affinity=sorted(os.sched_getaffinity(0)))
        save(ROOT/'benchmark'/('%s_h%g_r%d_repeat%d.json'%(kind,height,ranks,repeat)),result)
        print('BENCHMARK',json.dumps(result),flush=True)

if __name__=='__main__':
    (ROOT/'admission').mkdir(exist_ok=True)
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['gate','mpi','benchmark']);p.add_argument('--solver',default='basilisk');p.add_argument('--repeat',type=int,default=0);p.add_argument('--height',type=float,default=1.05)
    a=p.parse_args()
    if a.mode=='gate':gate()
    elif a.mode=='mpi':mpi_gate()
    else:benchmark(a.solver,a.repeat,a.height)
