"""AMS: strict first crossing, full-state cloning, actual tie kill count, bounded budget."""
from concurrent.futures import ThreadPoolExecutor
from collections import Counter
from pathlib import Path
import math
import random
from rare_events import core,cli


def cutoff(scores,D):
    assert 1<=D<len(scores) and all(math.isfinite(v) for v in scores)
    level=sorted(scores)[D-1]
    return level,[i for i,s in enumerate(scores) if s<=level]


def first_crossing(records,level):
    for j,r in enumerate(records):
        if 1-r['q3']>level:return j,r
    raise AssertionError('Surviving path has no strict crossing checkpoint')


def propagate(spec,config):
    case,initial_seed,end_step,restore,prefix,ancestor=spec
    out=core.execute(case,initial_seed,end_step,config,restore=restore,records=True)
    records=list(prefix)
    for line in (Path(case)/'records.tsv').read_text().splitlines():
        step,t,q,name=line.split('\t')
        records.append(dict(step=int(step),time=float(t),q3=float(q),checkpoint=str(Path(case)/name)))
    assert all(a['step']<b['step'] and a['q3']>b['q3'] for a,b in zip(records,records[1:]))
    summary=out['summary']
    assert not records or abs(records[-1]['q3']-summary['best_q3'])<1e-12
    return dict(case=str(case),score=1-summary['best_q3'],tau=summary['time'] if summary['reason']=='rupture' else None,
                records=records,checkpoint=out['checkpoint'],ancestor=ancestor)


def run(workdir,config,resume=False):
    lock=cli.begin(workdir,config,resume)
    try:
        N,D=config['N'],config['D'];assert 1<=D<N and config['max_rounds']>=0
        end=core.check_config(config);namespace=config['namespace']
        state_path=workdir/'state.json'
        state=core.load(state_path) if state_path.exists() else None
        if state is None:
            core.save(workdir/'status.json',dict(phase='initial',requested=N))
            specs=[(workdir/'initial'/('p%03d'%i),core.seed(namespace,'initial:%d'%i),end,None,[],i) for i in range(N)]
            with ThreadPoolExecutor(max_workers=config['workers']) as pool:
                particles=list(pool.map(lambda x:propagate(x,config),specs))
            state=dict(round=0,log_weight=0.,particles=particles)
            core.save(state_path,state);core.save(workdir/'state_000.json',state)
        if (workdir/'result.json').exists():return core.load(workdir/'result.json')
        while True:
            particles=state['particles'];level,dead=cutoff([p['score'] for p in particles],D)
            if level>=.99:
                result=core.finalize(particles,state['log_weight'],config,'ams');break
            if len(dead)==N:
                result=dict(status='complete',extinction=True,P=0.,RMST=config['horizon'],MFPT=None,
                            log_weight_before_extinction=state['log_weight'],tail_certified=False,
                            CDF_times=[config['horizon']*j/16 for j in range(17)],CDF=[0.]*17)
                break
            if state['round']>=config['max_rounds']:
                result=dict(status='budget_unresolved',rounds=state['round'],log_weight=state['log_weight'],
                            P=None,RMST=None,MFPT=None,horizon=config['horizon'])
                break
            survivors=[i for i in range(N) if i not in dead]
            rng=random.Random(core.seed(namespace,'parents:%d'%state['round']))
            specs=[];parents=[]
            for i in dead:
                parent=rng.choice(survivors);p=particles[parent];j,record=first_crossing(p['records'],level)
                parents.append(parent)
                specs.append((workdir/('round_%04d'%(state['round']+1))/('p%03d'%i),
                    core.seed(namespace,'branch:%d:%d'%(state['round']+1,i)),end,record['checkpoint'],
                    p['records'][:j+1],p['ancestor']))
            core.save(workdir/'status.json',dict(phase='splitting',round=state['round'],cut=level,
                      replacing=len(dead),log_weight=state['log_weight'],P=None,MFPT=None))
            with ThreadPoolExecutor(max_workers=config['workers']) as pool:
                children=list(pool.map(lambda x:propagate(x,config),specs))
            for i,p in zip(dead,children):particles[i]=p
            state=dict(round=state['round']+1,log_weight=state['log_weight']+math.log1p(-len(dead)/N),particles=particles)
            core.save(workdir/('selection_%04d.json'%state['round']),dict(round=state['round'],cut=level,
                dead=dead,parents=parents,log_weight=state['log_weight'],families=dict(Counter(p['ancestor'] for p in particles))))
            core.save(workdir/('state_%03d.json'%state['round']),state);core.save(state_path,state)
        result.update(config=config,rounds=state['round'])
        core.save(workdir/'result.json',result);core.save(workdir/'status.json',dict(phase='finished',result=result))
        return result
    finally:lock.close()


if __name__=='__main__':
    root,config,resume=cli.arguments('ams')
    print(run(root,config,resume))
