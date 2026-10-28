"""Independent adaptive trajectories sharing batched resident GPU solves."""
import hashlib,json,math,os,platform,time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import scipy
import cupy as cp
from accelerated_runner import Trajectory
from batched_film import BatchedFilm

ROOT=Path(__file__).resolve().parent

def fingerprint():
    names=['batch_runner.py','batched_film.py','operators.cu','midpoint_candidate.py',
           'adaptive_runner.py','accelerated_runner.py','native_midpoint_candidate.py','native_midpoint_kernel.c']
    return dict(sources={n:hashlib.sha256((ROOT.parent/({'operators.cu':'sle_smid_fv_cuda.cu','native_midpoint_kernel.c':'sle_smid_fv_c99.c'}[n]) if n.endswith(('.c','.cu')) else ROOT/n).read_bytes()).hexdigest() for n in names},
        python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,cupy=cp.__version__,
        cuda_runtime=cp.cuda.runtime.runtimeGetVersion(),cuda_driver=cp.cuda.runtime.driverGetVersion(),
        GPU=str(cp.cuda.runtime.getDeviceProperties(0)['name']),precision='FP64',fmad=False)

class Batch:
    def __init__(self,heights,seeds,n=64,dt_cap=.3,h=None,noise_workers=1,restart=8):
        assert len(heights)>0 and len(heights)==len(seeds) and len(set(seeds))==len(seeds)
        assert all(math.isfinite(v) and v>0 for v in heights)
        assert all(isinstance(s,(int,np.integer)) and 0<=s<2**64 for s in seeds)
        assert n>=4 and n%2==0 and 1<=noise_workers and 1<=restart<=25
        self.B=len(seeds);self.n=n;self.restart=restart;self.dt_cap=dt_cap
        initial=np.ones((self.B,n,n)) if h is None else np.asarray(h,dtype=np.float64)
        assert initial.shape==(self.B,n,n) and np.isfinite(initial).all() and initial.min()>.01
        self.s=[Trajectory(float(height),n,int(seed),dt_cap=dt_cap,h=initial[i],base_dt=.064,depth=18,kernel='reference') for i,(height,seed) in enumerate(zip(heights,seeds))]
        self.film=BatchedFilm(heights,n,restart=restart);self.h=cp.asarray(initial)
        self.minimum=initial.min(axis=(-2,-1));self._fingerprint=fingerprint()
        self.pool=ThreadPoolExecutor(max_workers=noise_workers) if noise_workers>1 else None
        self.noise_workers=noise_workers
        self.stats=dict(noise_seconds=0.,GPU_solve_seconds=0.,noise_upload_bytes=0,field_downloads=0,batch_calls=0)

    def close(self):
        if self.pool is not None:self.pool.shutdown();self.pool=None

    def fields(self):
        self.stats['field_downloads']+=1
        return self.h.get()

    def advance(self,stop=None):
        live=np.array([s.status=='RUNNING' and (stop is None or s.time<stop) for s in self.s])
        if not live.any():return False
        maximum=np.zeros(self.B,dtype=np.int64);planned=maximum.copy();proposed=maximum.copy()
        boundary=np.zeros(self.B,dtype=bool);rejected=np.zeros(self.B,dtype=np.int64)
        for b,s in enumerate(self.s):
            if not live[b]:continue
            cap=.003 if self.minimum[b]<.1 else (.01 if self.minimum[b]<.3 else self.dt_cap)
            maximum[b]=min(s.cap_ticks,s.noise.ticks(cap));planned[b]=min(maximum[b],s.proposal)
            proposed[b]=planned[b] if stop is None else min(planned[b],s.noise.ticks(stop)-s.ticks)
            boundary[b]=proposed[b]<planned[b]
        pending=live.copy()
        while pending.any():
            dt=np.zeros(self.B);noise=np.zeros((self.B,2,self.n,self.n))
            ids=list(np.flatnonzero(pending));start=time.perf_counter()
            def generate(b):
                s=self.s[b]
                return b,s.noise.increment(s.ticks,s.ticks+int(proposed[b]))
            generated=map(generate,ids) if self.pool is None else self.pool.map(generate,ids)
            for b,dw in generated:
                noise[b]=dw;dt[b]=int(proposed[b])*self.s[b].noise.quantum
            self.stats['noise_seconds']+=time.perf_counter()-start
            noise_device=cp.asarray(noise);self.stats['noise_upload_bytes']+=noise.nbytes
            start=time.perf_counter();candidate,diag=self.film.step(self.h,dt,noise_device,enabled=pending)
            self.stats['GPU_solve_seconds']+=time.perf_counter()-start;self.stats['batch_calls']+=1
            accepted=np.zeros(self.B,dtype=bool)
            for b in ids:
                s=self.s[b];s.ledger['attempts']+=1
                why=None
                if diag[b,2] or not diag[b,6]:why='GPU nonlinear solve failed or nonfinite endpoint'
                elif diag[b,3]<=0:why='Nonpositive endpoint; no repair'
                elif diag[b,3]<=.01 and proposed[b]>s.event_resolution:
                    why='Refine event bracket on the same Brownian path';s.ledger['event_refinements']+=1
                if why is not None:
                    rejected[b]+=1;s.ledger['rejections']+=1;s.ledger['last_error']=why
                    if proposed[b]<=1:s.status='NUMERICAL_FAILURE';pending[b]=False
                    else:proposed[b]=max(1,int(proposed[b])//2)
                    continue
                accepted[b]=True;pending[b]=False
                oldtime=s.time;s.ticks+=int(proposed[b]);s.step_count+=1;s.noise.commit(s.ticks)
                s.last_interval=[oldtime,s.time];self.minimum[b]=diag[b,3]
                s.ledger['max_residual']=max(s.ledger['max_residual'],float(diag[b,0]))
                s.ledger['max_newton']=max(s.ledger['max_newton'],int(diag[b,1]))
                s.ledger['max_relative_volume_error']=max(s.ledger['max_relative_volume_error'],abs(float(diag[b,5])/s.initial_mean-1))
                value=float(dt[b]);previous=s.ledger['min_accepted_dt']
                s.ledger['min_accepted_dt']=value if previous is None else min(previous,value)
                s.ledger['max_accepted_dt']=max(s.ledger['max_accepted_dt'],value)
                if diag[b,3]<=.01:
                    s.status='RUPTURE_BRACKET';s.event=dict(lower=oldtime,upper=s.time,hmin=float(diag[b,3]),step=s.step_count)
                s.proposal=int(planned[b] if boundary[b] and rejected[b]==0 else proposed[b])
                s.streak=0 if rejected[b] else s.streak+1
                if s.streak>=20:
                    s.proposal=min(int(maximum[b]),max(s.proposal+1,5*s.proposal//4));s.streak=0
            self.h=cp.where(cp.asarray(accepted)[:,None,None],candidate,self.h)
        return True

    def run(self,until,wall_seconds=3600,max_rounds=0):
        assert math.isfinite(until) and until>=0 and math.isfinite(wall_seconds) and wall_seconds>0 and max_rounds>=0
        if until>0:assert self.s[0].noise.ticks(until)>0 and abs(self.s[0].noise.ticks(until)*self.s[0].noise.quantum-until)<1e-12,'Unrepresentable end time'
        start=time.perf_counter();rounds=0
        while any(s.status=='RUNNING' and (until==0 or s.time<until) for s in self.s):
            if time.perf_counter()-start>=wall_seconds or (max_rounds and rounds>=max_rounds):break
            self.advance(stop=until if until>0 else None);rounds+=1
        cp.cuda.get_current_stream().synchronize()
        return time.perf_counter()-start

    def summary(self,until):
        records=[]
        for b,s in enumerate(self.s):
            state=s.status
            if state=='RUNNING':state='CENSORED_AT_TARGET' if until>0 and s.time>=until else 'PAUSED_WITH_CHECKPOINT'
            records.append(dict(height=s.film.height,seed=s.noise.seed,state=state,t=s.time,ticks=s.ticks,
                step=s.step_count,hmin=float(self.minimum[b]),event=s.event,ledger=s.ledger))
        return dict(implementation='Resident batched GPU midpoint',B=self.B,n=self.n,dt_cap=self.dt_cap,
            records=records,stats=self.stats,diagnostic_host_bytes=self.film.diagnostic_bytes,
            control_scalar_synchronizations=self.film.synchronizations,noise_workers=self.noise_workers,
            fingerprint=self._fingerprint)

    def save(self,path):
        path=Path(path);assert self._fingerprint==fingerprint(),'GPU code/runtime changed'
        h=self.fields();keys=[];values=[]
        for b,s in enumerate(self.s):
            for key in sorted(s.noise.cache):keys.append([b]+list(key));values.append(s.noise.cache[key])
        values=np.stack(values) if values else np.empty((0,2,self.n,self.n))
        meta=dict(schema='resident_GPU_batch_v1',fingerprint=self._fingerprint,n=self.n,dt_cap=self.dt_cap,
            restart=self.restart,states=[s.metadata() for s in self.s],keys=keys,stats=self.stats)
        payload=json.dumps(meta,sort_keys=True,separators=(',',':'))
        digest=hashlib.sha256(payload.encode()+h.tobytes()+values.tobytes()).hexdigest()
        temp=path.with_name(path.name+'.tmp')
        with temp.open('wb') as f:np.savez_compressed(f,h=h,noise_values=values,metadata=payload,sha256=digest)
        os.replace(str(temp),str(path))
        return h

    @classmethod
    def load(cls,path,expected_seeds=None,noise_workers=1):
        with np.load(path,allow_pickle=False) as z:
            h=z['h'].copy();values=z['noise_values'].copy();payload=str(z['metadata'].item());digest=str(z['sha256'].item())
        assert hashlib.sha256(payload.encode()+h.tobytes()+values.tobytes()).hexdigest()==digest,'Checkpoint SHA mismatch'
        meta=json.loads(payload);assert meta['schema']=='resident_GPU_batch_v1' and meta['fingerprint']==fingerprint(),'Checkpoint code/runtime mismatch'
        seeds=[d['seed'] for d in meta['states']]
        if expected_seeds is not None:assert list(expected_seeds)==seeds,'Checkpoint seed identities mismatch'
        # Constructor requires an above-threshold initial state. Restore exact
        # event fields afterwards; no physical field is repaired or evolved.
        initial=h if h.min()>.01 else np.ones_like(h)
        obj=cls([d['height'] for d in meta['states']],seeds,meta['n'],meta['dt_cap'],h=initial,
                noise_workers=noise_workers,restart=meta['restart'])
        obj.h=cp.asarray(h);obj.minimum=h.min(axis=(-2,-1));obj.stats=meta['stats']
        for s,d in zip(obj.s,meta['states']):
            assert s._fingerprint==d['fingerprint'] and s.noise.base_dt==d['base_dt'] and s.noise.depth==d['depth']
            for name in ('ticks','proposal','streak','step_count','event_resolution','status','event','initial_mean','ledger','last_interval'):
                setattr(s,name,d[name])
            s.noise.cursor=d['noise_cursor'];assert s.ticks==s.noise.cursor;s.noise.cache={}
        assert values.shape==(len(meta['keys']),2,obj.n,obj.n)
        for key,value in zip(meta['keys'],values):obj.s[key[0]].noise.cache[tuple(key[1:])]=value.copy()
        return obj
