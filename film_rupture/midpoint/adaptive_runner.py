"""Adaptive midpoint trajectory with replayable Brownian bridges and checkpoints.

Rejected steps reuse the same Brownian path. No floor/volume repair is applied.
This is an implementation gate, not a proof of unbiased first-passage statistics.
"""
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import numpy as np
import scipy
from midpoint_candidate import Film


class BrownianTape:
    def __init__(self,seed,shape,base_dt=.001,depth=12):
        assert 0<=int(seed)<2**64 and 1<=depth<=20
        self.seed=int(seed);self.shape=tuple(shape);self.base_dt=float(base_dt);self.depth=int(depth)
        self.width=1<<self.depth
        self.quantum=self.base_dt/self.width
        self.cursor=0
        self.cache={}

    def ticks(self,time):
        return int(round(float(time)/self.quantum))

    def _normal(self,cell,level,node,kind):
        assert 0<=cell<2**128
        entropy=[self.seed&0xffffffff,self.seed>>32]+[(cell>>(32*i))&0xffffffff for i in range(4)]+[level,node,kind]
        return np.random.default_rng(np.random.SeedSequence(entropy)).standard_normal(self.shape)

    def node(self,cell,level,index):
        key=(cell,level,index)
        if key in self.cache:
            return self.cache[key]
        if level==0:
            self.cache[key]=math.sqrt(self.base_dt)*self._normal(cell,0,0,1)
        else:
            parent=self.node(cell,level-1,index//2)
            left=(cell,level,2*(index//2));right=(cell,level,2*(index//2)+1)
            if left in self.cache:
                self.cache[right]=parent-self.cache[left]
            elif right in self.cache:
                self.cache[left]=parent-self.cache[right]
            else:
                duration=self.base_dt/(1<<(level-1))
                self.cache[left]=.5*parent+math.sqrt(duration/4)*self._normal(cell,level-1,index//2,2)
                self.cache[right]=parent-self.cache[left]
        return self.cache[key]

    def _inside(self,cell,a,b,level=0,index=0,start=0,width=None):
        width=self.width if width is None else width
        if a==start and b==start+width:
            return self.node(cell,level,index)
        half=width//2
        assert half>0
        middle=start+half
        if b<=middle:
            return self._inside(cell,a,b,level+1,2*index,start,half)
        if a>=middle:
            return self._inside(cell,a,b,level+1,2*index+1,middle,half)
        return (self._inside(cell,a,middle,level+1,2*index,start,half)+
                self._inside(cell,middle,b,level+1,2*index+1,middle,half))

    def increment(self,start,end):
        assert self.cursor<=start<end
        result=np.zeros(self.shape)
        for cell in range(start//self.width,(end-1)//self.width+1):
            a=max(start-cell*self.width,0);b=min(end-cell*self.width,self.width)
            result+=self._inside(cell,a,b)
        return result

    def commit(self,end):
        assert end>=self.cursor
        self.cursor=int(end)
        self.cache={key:value for key,value in self.cache.items() if (key[0]+1)*self.width>self.cursor}


def source_hashes():
    root=Path(__file__).resolve().parent
    return {name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ('adaptive_runner.py','midpoint_candidate.py')}


def runtime():
    return dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)


class Trajectory:
    def __init__(self,height,n,seed,dt_cap=.03,h=None,tiered=True,event_resolution=1e-5,depth=12):
        self.film=Film(height,n)
        self.h=np.ones((n,n)) if h is None else np.array(h,dtype=float,copy=True)
        assert self.h.shape==(n,n) and np.all(np.isfinite(self.h)) and self.h.min()>.01
        self.noise=BrownianTape(seed,(2,n,n),depth=depth)
        self.cap_ticks=self.noise.ticks(dt_cap)
        assert self.cap_ticks>0 and abs(self.cap_ticks*self.noise.quantum-dt_cap)<1e-12
        self.proposal=self.cap_ticks;self.streak=0;self.step_count=0;self.ticks=0
        self.tiered=bool(tiered)
        self.event_resolution=max(1,self.noise.ticks(event_resolution))
        self.status='RUNNING';self.event=None
        self.initial_mean=float(self.h.mean())
        self.ledger=dict(attempts=0,rejections=0,event_refinements=0,max_residual=0.,max_newton=0,
                         max_relative_volume_error=0.,min_accepted_dt=None,max_accepted_dt=0.,
                         floor_cells=0,floor_added_volume=0.,last_error=None)
        self.last_interval=None

    @property
    def time(self):
        return self.ticks*self.noise.quantum

    def advance(self,stop_ticks=None):
        assert self.status=='RUNNING'
        maximum=self.cap_ticks
        if self.tiered:
            cap=.003 if self.h.min()<.1 else (.01 if self.h.min()<.3 else self.cap_ticks*self.noise.quantum)
            maximum=min(maximum,self.noise.ticks(cap))
        planned=min(maximum,self.proposal)
        proposed=planned
        if stop_ticks is not None:
            assert stop_ticks>self.ticks
            proposed=min(proposed,int(stop_ticks)-self.ticks)
        planned_boundary_clip=proposed<planned
        rejected=0
        while True:
            dt=proposed*self.noise.quantum
            dw=self.noise.increment(self.ticks,self.ticks+proposed)
            self.ledger['attempts']+=1
            why=None
            try:
                new,diagnostic=self.film.step(self.h,dt,dw)
                if not np.all(np.isfinite(new)) or new.min()<=0:
                    why='Nonpositive/nonfinite endpoint; no floor repair'
                elif new.min()<=.01 and proposed>self.event_resolution:
                    why='Refine event bracket on the same Brownian path'
                    self.ledger['event_refinements']+=1
            except (ValueError,RuntimeError,FloatingPointError) as exc:
                why=str(exc)
            if why is None:
                break
            self.ledger['rejections']+=1
            self.ledger['last_error']=why
            rejected+=1
            if proposed<=1:
                self.status='NUMERICAL_FAILURE'
                return dict(state=self.status,time=self.time,error=why,rejected=rejected)
            proposed=max(1,proposed//2)
        old_time=self.time
        self.h=new;self.ticks+=proposed;self.step_count+=1
        self.noise.commit(self.ticks)
        self.last_interval=[old_time,self.time]
        self.ledger['max_residual']=max(self.ledger['max_residual'],diagnostic['residual'])
        self.ledger['max_newton']=max(self.ledger['max_newton'],diagnostic['newton'])
        self.ledger['max_relative_volume_error']=max(self.ledger['max_relative_volume_error'],abs(float(new.mean())/self.initial_mean-1))
        self.ledger['min_accepted_dt']=dt if self.ledger['min_accepted_dt'] is None else min(self.ledger['min_accepted_dt'],dt)
        self.ledger['max_accepted_dt']=max(self.ledger['max_accepted_dt'],dt)
        if new.min()<=.01:
            self.status='RUPTURE_BRACKET'
            self.event=dict(lower=old_time,upper=self.time,hmin=float(new.min()),step=self.step_count)
        self.proposal=planned if planned_boundary_clip and rejected==0 else proposed
        self.streak=0 if rejected else self.streak+1
        if self.streak>=20:
            self.proposal=min(maximum,max(self.proposal+1,5*self.proposal//4));self.streak=0
        return dict(state=self.status,time=self.time,dt=dt,hmin=float(new.min()),
                    residual=diagnostic['residual'],newton=diagnostic['newton'],rejected=rejected)

    def run_steps(self,count):
        records=[]
        for _ in range(count):
            if self.status!='RUNNING':break
            records.append(self.advance())
        return records

    def fork(self,seed):
        child=copy.deepcopy(self)
        child.noise.seed=int(seed)
        assert 0<=child.noise.seed<2**64
        return child

    def metadata(self):
        return dict(schema=1,height=self.film.height,n=self.film.n,A=self.film.a,phi=self.film.phi,slip=self.film.slip,
            pressure=self.film.pressure,seed=self.noise.seed,base_dt=self.noise.base_dt,depth=self.noise.depth,
            ticks=self.ticks,noise_cursor=self.noise.cursor,cap_ticks=self.cap_ticks,proposal=self.proposal,
            streak=self.streak,step_count=self.step_count,tiered=self.tiered,event_resolution=self.event_resolution,
            status=self.status,event=self.event,initial_mean=self.initial_mean,ledger=self.ledger,
            last_interval=self.last_interval,sources=source_hashes(),runtime=runtime())

    def save(self,path):
        path=Path(path)
        keys=sorted(self.noise.cache)
        values=np.stack([self.noise.cache[k] for k in keys]) if keys else np.empty((0,*self.noise.shape))
        metadata=self.metadata();metadata['noise_keys']=keys
        payload=json.dumps(metadata,sort_keys=True,separators=(',',':'))
        digest=hashlib.sha256(payload.encode()+self.h.tobytes()+values.tobytes()).hexdigest()
        temporary=path.with_name(path.name+'.tmp')
        with temporary.open('wb') as f:
            np.savez_compressed(f,h=self.h,noise_values=values,metadata=payload,sha256=digest)
        os.replace(str(temporary),str(path))

    @classmethod
    def load(cls,path):
        with np.load(path,allow_pickle=False) as archive:
            h=archive['h'].copy();values=archive['noise_values'].copy();payload=str(archive['metadata'].item())
            expected=str(archive['sha256'].item())
        assert hashlib.sha256(payload.encode()+h.tobytes()+values.tobytes()).hexdigest()==expected,'Checkpoint SHA mismatch'
        data=json.loads(payload)
        assert data['schema']==1 and data['sources']==source_hashes(),'Checkpoint source mismatch'
        assert data['runtime']==runtime(),'Checkpoint runtime mismatch; cross-version replay not implicitly accepted'
        obj=cls(data['height'],data['n'],data['seed'],dt_cap=data['cap_ticks']*(data['base_dt']/(1<<data['depth'])),
                h=h if h.min()>.01 else np.ones_like(h),tiered=data['tiered'],depth=data['depth'])
        assert obj.film.a==data['A'] and obj.film.phi==data['phi'] and obj.film.slip==data['slip'] and obj.film.pressure==data['pressure']
        obj.h=h
        for name in ('ticks','proposal','streak','step_count','event_resolution','status','event','initial_mean','ledger','last_interval'):
            setattr(obj,name,data[name])
        obj.noise.cursor=data['noise_cursor']
        obj.noise.cache={tuple(key):value.copy() for key,value in zip(data['noise_keys'],values)}
        assert obj.ticks==obj.noise.cursor
        return obj
