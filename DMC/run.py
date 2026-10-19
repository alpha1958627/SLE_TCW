"""DMC and QDMC on the same physical engine and q3 score as AMS."""
from concurrent.futures import ThreadPoolExecutor
from collections import Counter
import math
import random
from rare_events import core,cli


def run(workdir,config,resume=False):
    lock=cli.begin(workdir,config,resume)
    try:
        core.check_config(config)
        N=config['N'];block=config['block'];stages=round(config['horizon']/block)
        assert block>0 and abs(stages*block-config['horizon'])<1e-8 and config['strength']>=0
        state_path=workdir/'state.json';namespace=config['namespace']
        state=core.load(state_path) if state_path.exists() else dict(stage=0,logz=0.,particles=[
            dict(seed=core.seed(namespace,'initial:%d'%i),checkpoint=None,ancestor=i,tau=None,V=0.,score=0.) for i in range(N)])
        if not state_path.exists():core.save(state_path,state)
        if (workdir/'result.json').exists():return core.load(workdir/'result.json')
        for k in range(state['stage']+1,int(stages)+1):
            now=k*block;end=round(now/config['dt']);assert abs(end*config['dt']-now)<1e-8
            pending=[(i,p) for i,p in enumerate(state['particles']) if p['tau'] is None]
            core.save(workdir/'status.json',dict(phase='propagating',completed_stage=state['stage'],
                horizon_completed=state['stage']*block,target_time=now,planned_horizon=config['horizon']))
            def execute(item):
                i,p=item
                return core.execute(workdir/('stage_%04d'%k)/('p%03d'%i),p['seed'],int(end),config,restore=p['checkpoint'])
            with ThreadPoolExecutor(max_workers=config['workers']) as pool:outputs=list(pool.map(execute,pending))
            propagated=[dict(p) for p in state['particles']]
            for (i,p),out in zip(pending,outputs):
                summary=out['summary'];propagated[i].update(checkpoint=out['checkpoint'],score=1-summary['q3'])
                if summary['reason']=='rupture':propagated[i]['tau']=summary['time']
            old=[p['V'] for p in state['particles']]
            new=core.potential([p['score'] for p in propagated],old,config['mode'],config['strength']*now/config['horizon'])
            weights,increment=core.selection(old,new)
            parents=core.resample(weights,random.Random(core.seed(namespace,'selection:%d'%k)).random()/N)
            counts=Counter(parents);children=[]
            for slot,parent in enumerate(parents):
                p=dict(propagated[parent]);p['V']=new[parent]
                if counts[parent]>1:p['seed']=core.seed(namespace,'branch:%d:%d'%(k,slot))
                children.append(p)
            assert len({p['seed'] for p in children})==N
            logz=state['logz']+increment;stats=core.statistics(children,logz,now)
            stats['time_reference_factor']=config['ratio']
            families=dict(Counter(p['ancestor'] for p in children))
            ledger=dict(stage=k,time=now,parents=parents,previous_V=old,new_V=new,selection_weights=weights,
                        log_mean_increment=increment,selection_ESS=1/math.fsum(w*w for w in weights),
                        families=families,propagated=propagated,statistics=stats)
            core.save(workdir/('selection_%04d.json'%k),ledger)
            state=dict(stage=k,logz=logz,particles=children,statistics=stats)
            core.save(state_path,state)
            core.save(workdir/'status.json',dict(phase='selected',stage=k,time=now,statistics=stats,
                      families=families,selection_ESS=ledger['selection_ESS']))
        result=core.finalize(state['particles'],state['logz'],config)
        result.update(config=config,stages=int(stages))
        core.save(workdir/'result.json',result);core.save(workdir/'status.json',dict(phase='finished',result=result))
        return result
    finally:lock.close()


if __name__=='__main__':
    root,config,resume=cli.arguments('dmc')
    print(run(root,config,resume))
