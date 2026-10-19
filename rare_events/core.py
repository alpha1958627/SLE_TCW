"""Shared process identity, checkpoint validation and compensated resampling.

Python 3.6; standard library only. A partial estimate is never an MFPT.
"""
from pathlib import Path
from collections import Counter
import hashlib
import json
import math
import os
import random
import struct
import subprocess
import time

ROOT = Path(__file__).resolve().parent
HEADER = struct.Struct('<8s64sIIQQ11d')


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = Path(str(path)+'.tmp')
    temp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
    os.replace(str(temp), str(path))


def load(path):
    return json.loads(Path(path).read_text())


def seed(namespace, label):
    return int(hashlib.sha256((namespace+':'+label).encode()).hexdigest()[:16],16) or 1


def check_config(config):
    assert config['N']>=2 and config['ratio']>0 and config['slip_ref']>=0
    assert config['dt']>0 and config['horizon']>0 and config['ranks'] in (1,4,16)
    assert config['repulsion'] is False and config['score']=='q3'
    n=round(config['horizon']/config['dt'])
    assert abs(n*config['dt']-config['horizon'])<1e-8 and n<=2147483647
    assert config['namespace'], 'A new campaign must use a new independent seed namespace'
    return int(n)


def build_identity():
    data=load(ROOT/'build/build.json')
    for name,digest in data['files'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest, name+' changed'
    return data


def checkpoint(path, config=None, signature=None):
    raw=Path(path).read_bytes()
    assert len(raw)==HEADER.size+4096*8
    header=HEADER.unpack_from(raw)
    magic,tag,n,version,step,source_seed=header[:6]
    assert magic==b'RESTFE26' and n==64 and version==1
    if signature:assert tag.decode()==signature
    dt,ratio,slip,phi,A,L,event,floor,mass0,bestq,clamp=header[6:]
    assert dt>0 and ratio>0 and event==.01 and floor==.005
    if config:
        assert (dt,ratio,slip)==(config['dt'],config['ratio'],config['slip_ref'])
        assert abs(phi/(.01/ratio**2)-1)<1e-14 and abs(A/(.0013159/ratio**2)-1)<1e-14 and L==100/ratio
    field=struct.unpack_from('<4096d',raw,HEADER.size)
    assert all(math.isfinite(v) and v>=floor for v in field)
    return dict(step=step,time=step*dt,seed=source_seed,field=field,dt=dt,ratio=ratio,
                mass0=mass0,best_q3=bestq,clamp_mass=clamp,signature=tag.decode())


def execute(case, particle_seed, end_step, config, restore=None, records=False, binary='mpi', ranks=None):
    case=Path(case);case.mkdir(parents=True,exist_ok=True)
    identity=build_identity()
    spec=dict(seed=particle_seed,end_step=end_step,config=config,restore=restore,
              records=records,binary=binary,ranks=ranks or config['ranks'],engine_signature=identity['engine_signature'])
    request=case/'request.json'
    if request.exists():assert load(request)==spec, 'Case identity changed'
    else:save(request,spec)
    if not (case/'execution.json').exists():
        assert not (case/'stdout.log').exists(), 'Unconfirmed attempt: preserve files; do not replay automatically'
        command=[] if binary=='serial' else ['srun','--input=none','--exclusive','--mem=256M','--mpi=pmix',
                                            '-n',str(spec['ranks']),'-c','1']
        command += [str(ROOT/'build'/binary),'--seed',str(particle_seed),'--end-step',str(end_step),
                    '--dt',str(config['dt']),'--ratio',str(config['ratio']),'--slip',str(config['slip_ref'])]
        if restore:command+=['--restore',str(restore)]
        if records:command+=['--records']
        began=time.time()
        with (case/'stdout.log').open('w') as out,(case/'stderr.log').open('w') as err:
            result=subprocess.run(command,cwd=str(case),stdout=out,stderr=err)
        save(case/'execution.json',dict(command=command,returncode=result.returncode,wall_seconds=time.time()-began))
    exe=load(case/'execution.json')
    assert exe['returncode']==0, 'Solver failure retained at '+str(case)
    summary=load(case/'summary.json')
    assert summary['seed']==particle_seed and summary['reason'] in ('horizon','rupture')
    data=checkpoint(case/'final.chk',config,identity['engine_signature'])
    assert data['step']==summary['step']<=end_step and abs(data['time']-summary['time'])<1e-8
    assert abs(min(data['field'])-summary['hmin'])<1e-12
    assert abs(summary['mass_relative_drift'])<=1e-5
    return dict(case=str(case),checkpoint=str(case/'final.chk'),summary=summary,execution=exe)


def inverse_normal(u):
    assert 0<u<1
    low,high=-9.,9.
    for unused in range(64):
        middle=(low+high)/2
        if .5*math.erfc(-middle/math.sqrt(2))<u:low=middle
        else:high=middle
    return (low+high)/2


def potential(scores,previous,mode,strength):
    assert mode in ('dmc','qdmc') and len(scores)==len(previous)>0 and strength>=0
    assert all(math.isfinite(v) for v in scores+previous)
    if mode=='dmc':return [strength*s for s in scores]
    shift=max(-v for v in previous)
    w=[math.exp(-v-shift) for v in previous];total=math.fsum(w);w=[v/total for v in w]
    result=[]
    for s in scores:
        u=math.fsum(v*((x<s)+.5*(x==s)) for v,x in zip(w,scores))
        # Bound the potential quantile only; never clip or renormalize a probability estimate.
        result.append(strength*inverse_normal(max(1e-12,min(1-1e-12,u))))
    return result


def selection(previous,current):
    delta=[b-a for a,b in zip(previous,current)];shift=max(delta)
    g=[math.exp(v-shift) for v in delta];total=math.fsum(g)
    return [v/total for v in g],shift+math.log(total/len(g))


def resample(weights,offset):
    n=len(weights);assert n>0 and 0<=offset<1/n and abs(math.fsum(weights)-1)<1e-12
    parents=[];j=0;cum=weights[0]
    for i in range(n):
        while j<n-1 and offset+i/n>=cum:j+=1;cum+=weights[j]
        parents.append(j)
    return parents


def statistics(particles,logz,horizon,kind='dmc'):
    n=len(particles)
    weights=[math.exp(logz-(p['V'] if kind!='ams' else 0))/n for p in particles]
    probability=math.fsum(w for p,w in zip(particles,weights) if p['tau'] is not None)
    deficit=math.fsum(w*(horizon-p['tau']) for p,w in zip(particles,weights) if p['tau'] is not None)
    assert all(p['tau'] is None or 0<=p['tau']<=horizon for p in particles)
    result=dict(horizon=horizon,P=probability,RMST=horizon-deficit,MFPT=None,
                log_normalizer=logz,
                out_of_physical_range=probability>1 or horizon-deficit<0,
                uncertainty='independent whole runs required; clones are dependent',
                time_reference_factor=None,tail_certified=False)
    result['path_weight_sum' if kind=='ams' else 'random_mass']=math.fsum(weights)
    return result


def cdf(particles,weights,grid):
    return [math.fsum(w for p,w in zip(particles,weights) if p['tau'] is not None and p['tau']<=t) for t in grid]


def finalize(particles,logz,config,kind='dmc'):
    result=statistics(particles,logz,config['horizon'],kind)
    result['time_reference_factor']=config['ratio']
    result['RMST_reference']=config['ratio']*result['RMST']
    grid=[config['horizon']*j/16 for j in range(17)]
    weights=[math.exp(logz-(p['V'] if kind!='ams' else 0))/len(particles) for p in particles]
    result.update(status='complete',CDF_times=grid,CDF=cdf(particles,weights,grid))
    return result
