"""Independent direct-MC finite-window baseline; no splitting or resampling."""
from concurrent.futures import ThreadPoolExecutor
import math
from . import core,cli


def run(workdir,config,resume=False):
    lock=cli.begin(workdir,config,resume)
    try:
        if (workdir/'result.json').exists():return core.load(workdir/'result.json')
        end=core.check_config(config)
        core.save(workdir/'status.json',dict(phase='independent_paths',N=config['N'],horizon=config['horizon']))
        def execute(i):
            return core.execute(workdir/('p%03d'%i),core.seed(config['namespace'],'initial:%d'%i),end,config)
        with ThreadPoolExecutor(max_workers=config['workers']) as pool:outputs=list(pool.map(execute,range(config['N'])))
        particles=[dict(tau=x['summary']['time'] if x['summary']['reason']=='rupture' else None,ancestor=i)
                   for i,x in enumerate(outputs)]
        result=core.finalize(particles,0.,config,'ams')
        result.update(config=config,censored=sum(p['tau'] is None for p in particles),
            uncertainty='fixed-N independent paths; zero observed SEM under complete censoring does not certify precision')
        if result['censored']==config['N']:
            result['fixed_N_one_sided95_event_upper']=1-.05**(1/config['N'])
            result['fixed_N_one_sided95_RMST_lower']=config['horizon']*.05**(1/config['N'])
        core.save(workdir/'state.json',dict(particles=particles))
        core.save(workdir/'result.json',result);core.save(workdir/'status.json',dict(phase='finished',result=result))
        return result
    finally:lock.close()


if __name__=='__main__':
    root,config,resume=cli.arguments('direct')
    print(run(root,config,resume))
