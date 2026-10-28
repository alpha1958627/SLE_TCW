"""CLI for resident batched GPU direct-MC midpoint trajectories."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
import argparse,datetime,fcntl,hashlib,json,math,time,traceback
from pathlib import Path
import numpy as np
from batch_runner import Batch,fingerprint

def atomic(p,value):
    temp=p.with_name(p.name+'.tmp');temp.write_text(json.dumps(value,indent=2));os.replace(str(temp),str(p))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--heights',nargs='+',type=float,default=[1.2])
    parser.add_argument('--seeds',nargs='+',type=int,required=True)
    parser.add_argument('--n',type=int,default=64);parser.add_argument('--dt-max',type=float,default=.3)
    parser.add_argument('--until',type=float,default=300.,help='Fixed code-time target; 0 means event-only')
    parser.add_argument('--wall-seconds',type=float,default=3600.)
    parser.add_argument('--frame-interval',type=float,default=300.)
    parser.add_argument('--checkpoint-rounds',type=int,default=1000)
    parser.add_argument('--noise-workers',type=int,default=1)
    parser.add_argument('--resume',action='store_true')
    a=parser.parse_args();heights=a.heights*len(a.seeds) if len(a.heights)==1 else a.heights
    assert len(heights)==len(a.seeds) and a.checkpoint_rounds>0
    assert math.isfinite(a.frame_interval) and a.frame_interval>0 and math.isfinite(a.wall_seconds) and a.wall_seconds>0
    assert math.isfinite(a.until) and a.until>=0
    a.out.mkdir(parents=True,exist_ok=True)
    with (a.out/'.run.lock').open('a') as guard:
        fcntl.flock(guard,fcntl.LOCK_EX|fcntl.LOCK_NB)
        config=dict(heights=heights,seeds=a.seeds,n=a.n,dt_max=a.dt_max,until=a.until,
            frame_interval=a.frame_interval,checkpoint_rounds=a.checkpoint_rounds,fingerprint=fingerprint(),
            driver_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        protocol=a.out/'protocol.json';checkpoint=a.out/'checkpoint.npz'
        if a.resume:
            assert json.loads(protocol.read_text())==config,'Resume protocol mismatch'
            batch=Batch.load(checkpoint,expected_seeds=a.seeds,noise_workers=a.noise_workers)
        else:
            assert not protocol.exists(),'Output exists: explicit --resume required'
            batch=Batch(heights,a.seeds,a.n,a.dt_max,noise_workers=a.noise_workers)
            atomic(protocol,config)
        start=time.perf_counter();interval=batch.s[0].noise.ticks(a.frame_interval)
        assert interval>0 and abs(interval*batch.s[0].noise.quantum-a.frame_interval)<1e-12
        try:
            while any(s.status=='RUNNING' and (a.until==0 or s.time<a.until) for s in batch.s):
                remaining=a.wall_seconds-(time.perf_counter()-start)
                if remaining<=0:break
                next_ticks=min((s.ticks//interval+1)*interval for s in batch.s if s.status=='RUNNING' and (a.until==0 or s.time<a.until))
                stop=next_ticks*batch.s[0].noise.quantum
                if a.until>0:stop=min(stop,a.until)
                batch.run(stop,wall_seconds=remaining,max_rounds=a.checkpoint_rounds)
                h=batch.save(checkpoint)
                times=np.array([s.time for s in batch.s]);power=abs(np.fft.fft2(h-h.mean(axis=(-2,-1))[:,None,None],axes=(-2,-1))/a.n**2)**2
                # A ruptured lane retains its event time; never relabel it as
                # this frame's requested target time.
                temp=a.out/'last_field.tmp'
                with temp.open('wb') as f:np.savez_compressed(f,h=h,power=power,t=times,seeds=a.seeds,heights=heights)
                os.replace(str(temp),str(a.out/'last_field.npz'))
                if any(s.ticks==next_ticks for s in batch.s):
                    (a.out/('frame_%020d.npz'%next_ticks)).write_bytes((a.out/'last_field.npz').read_bytes())
                report=batch.summary(a.until)
                report.update(at=datetime.datetime.now(datetime.timezone.utc).isoformat(),wall_seconds=time.perf_counter()-start,
                    checkpoint_SHA256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
                    field_SHA256=hashlib.sha256((a.out/'last_field.npz').read_bytes()).hexdigest())
                atomic(a.out/'summary.json',report)
                print(json.dumps(dict(times=times.tolist(),states=[r['state'] for r in report['records']],wall=report['wall_seconds'])),flush=True)
            return 1 if any(s.status=='NUMERICAL_FAILURE' for s in batch.s) else 0
        finally:batch.close()

if __name__=='__main__':
    try:raise SystemExit(main())
    except Exception:
        traceback.print_exc();raise SystemExit(2)
