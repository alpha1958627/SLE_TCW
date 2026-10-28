"""Independent finite-chain enumeration plus continuous AMS analytic probability."""
from itertools import product
from pathlib import Path
import math
import random
import statistics
from . import core
from rare_events.ams import cutoff,first_crossing


def offsets(weights):
    n=len(weights);edges={0.,1/n};total=0.
    for weight in weights[:-1]:
        total+=weight
        for j in range(n):
            x=total-j/n
            if 0<x<1/n:edges.add(x)
    edges=sorted(edges)
    result=[(n*(b-a),core.resample(weights,(a+b)/2)) for a,b in zip(edges,edges[1:]) if b-a>1e-15]
    assert abs(sum(p for p,parents in result)-1)<1e-13
    return result


def exact_chain(mode):
    answers=[]
    def advance(population,logz,k,prob):
        if k==3:
            stat=core.statistics(population,logz,2.)
            answers.append((prob,stat,sum(p['tau'] is not None for p in population)/2))
            return
        hazard=.02 if k==1 else .03
        choices=[]
        for p in population:
            if p['tau'] is not None:choices.append([(1.,dict(p))])
            else:
                a=dict(p);b=dict(p);b['tau']=float(k);b['x']=1.
                choices.append([(1-hazard,a),(hazard,b)])
        for combination in product(*choices):
            moved=[p for q,p in combination];density=math.prod([q for q,p in combination]) if hasattr(math,'prod') else combination[0][0]*combination[1][0]
            old=[p['V'] for p in population]
            new=core.potential([p['x'] for p in moved],old,mode,2*k/2)
            weights,increment=core.selection(old,new)
            for p_offset,parents in offsets(weights):
                children=[]
                for parent in parents:
                    child=dict(moved[parent]);child['V']=new[parent];children.append(child)
                advance(children,logz+increment,k+1,prob*density*p_offset)
    advance([dict(tau=None,x=0.,V=0.) for unused in range(2)],0.,1,1.)
    values={key:math.fsum(p*stat[key] for p,stat,naive in answers) for key in ('P','RMST','random_mass')}
    values['uncompensated_P']=math.fsum(p*naive for p,stat,naive in answers)
    expected=dict(P=.02+.98*.03,RMST=1.98,random_mass=1.)
    assert abs(math.fsum(p for p,stat,naive in answers)-1)<1e-12
    for key,truth in expected.items():assert abs(values[key]-truth)<1e-12,(mode,key,values)
    assert abs(values['uncompensated_P']-expected['P'])>1e-4
    return values


def self_check():
    assert cutoff([.1,.1,.2,.4],1)==(.1,[0,1])
    assert first_crossing([dict(q3=.8),dict(q3=.7)],.2)[0]==1
    assert core.cdf([dict(tau=2.),dict(tau=None)],[.25,.25],[0.,2.,10.])==[0.,.25,.25]
    assert core.statistics([dict(tau=2.,V=0.),dict(tau=None,V=0.)],math.log(.5),10.)['RMST']==8.
    chains={mode:exact_chain(mode) for mode in ('dmc','qdmc')}
    # Analytic conditional mutation: X~U(0,1), event X>=.95, p=.05.
    rng=random.Random(2026101907);estimates=[];N=12
    for unused in range(2000):
        values=[rng.random() for j in range(N)];weight=1.
        while True:
            level,dead=cutoff(values,1)
            if level>=.95:
                estimates.append(weight*sum(x>=.95 for x in values)/N);break
            assert len(dead)==1
            weight*=1-1/N
            values[dead[0]]=level+(1-level)*rng.random()
    mean=statistics.mean(estimates);sem=statistics.stdev(estimates)/math.sqrt(len(estimates))
    assert abs(mean-.05)<4*sem,(mean,sem)
    return dict(status='PASS',finite_chain=chains,chain_truth=dict(P=.0494,RMST=1.98,mass=1.),
                ams_uniform_truth=.05,ams_mean=mean,ams_sem=sem,ams_independent_repetitions=2000,
                scope='probability accounting only; physical solver validation is separate')


if __name__=='__main__':
    import json
    print(json.dumps(self_check(),indent=2))
