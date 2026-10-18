"""Read existing server logs and rebuild the same local dashboard; no uploads."""
import csv
import io
import json
import math
import os
import shlex
from pathlib import Path
import subprocess
import sys
import tarfile
from monitor_volume import read_log

root=Path(__file__).resolve().parent
app=Path(os.environ.get('DASHBOARD_APP', str(root/'dashboard')))
ledger=json.loads((root/'deployment.json').read_text())
remote=ledger['base']
(root/'audits').mkdir(exist_ok=True)
errors=[]
for site, settings in ledger['sites'].items():
    host, python = settings['host'], settings.get('python', 'python3')
    ssh=['ssh','-o','BatchMode=yes','-o','ConnectTimeout=15',host]
    result=subprocess.run(ssh+[shlex.quote(python)+' '+shlex.quote(remote+'/monitor_volume.py')],capture_output=True,text=True,timeout=120)
    if result.returncode:
        errors.append(dict(site=site,error=result.stderr[-2500:]))
        continue
    data=json.loads(result.stdout)
    specs=[r for r in ledger['tasks'] if r['server']==site]
    jobs=','.join(settings['jobs'])
    assert all(job.isdigit() for job in settings['jobs']), 'Set numeric Slurm job IDs in deployment.json'
    state_result=subprocess.run(ssh+['squeue -r -j '+jobs+" -h -o '%i|%T|%M|%S|%e|%R'; sacct -X -j "+jobs+' -P -n --format=JobID,State,ExitCode,Elapsed,Start,End,NodeList,AllocCPUS'],capture_output=True,text=True,timeout=60)
    if state_result.returncode:
        errors.append(dict(site=site,error=state_result.stderr[-2500:]))
        continue
    states={}
    for line in state_result.stdout.splitlines():
        p=line.split('|')
        if len(p)>=2:
            states[p[0]]=p[1].strip()
    names=[r['version']+'/h0_'+r['h0_ratio']+'/run_0'+r['run']+'/log_restart' for r in specs]
    archive=subprocess.run(ssh+['tar -czf - -C '+shlex.quote(remote)+' '+' '.join(shlex.quote(name) for name in names)],capture_output=True,timeout=120)
    if archive.returncode:
        errors.append(dict(site=site,error=archive.stderr.decode(errors='replace')[-2500:]))
        continue
    local_dir=root/'audits'/'raw_logs'/site
    local_dir.mkdir(parents=True,exist_ok=True)
    with tarfile.open(fileobj=io.BytesIO(archive.stdout),mode='r:gz') as tf:
        members=tf.getmembers()
        assert {m.name for m in members}==set(names)
        for spec,name in zip(specs,names):
            member=tf.getmember(name)
            assert member.isfile() and not member.issym()
            target=local_dir/('task_'+spec['task']+'.log')
            target.write_bytes(tf.extractfile(member).read())
    kept=[]
    for spec in specs:
        rec=next(r for r in data['runs'] if r['task']==spec['task'])
        points,_=read_log(local_dir/('task_'+spec['task']+'.log'))
        if len(points)>2500:
            width=int(math.ceil(len(points)/1000.))
            retained={}
            for i in range(0,len(points),width):
                group=points[i:i+width]
                for item in [group[0],group[-1],min(group,key=lambda r:r['hmin']),max(group,key=lambda r:r['hmin'])]:
                    retained[item['step']]=item
            points=[retained[k] for k in sorted(retained)]
            rec['trace_sampling']='bucket_first_last_min_max_hmin'
        else:
            rec['trace_sampling']='every_available_log_row'
        rec['trace']=points
        # Metric snapshot and raw trajectory reads have separate explicit endpoints.
        rec['scheduler_state']=states.get(spec['job_id'],'UNKNOWN')
        rec['retained_job_id']=spec['job_id']
        kept.append(rec)
    data['runs']=kept
    data['scheduler_complete']=state_result.stdout
    target=root/'audits'/(site+'_final.json')
    temp=target.with_suffix('.json.tmp')
    temp.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    temp.replace(target)
    (root/'audits'/(site+'_scheduler.txt')).write_text(state_result.stdout,encoding='utf-8')
if errors:
    (root/'audits'/'refresh_errors.json').write_text(json.dumps(errors,indent=2)+'\n',encoding='utf-8')
    raise SystemExit('Refresh failed; previous dashboard retained: '+str(errors))
subprocess.run([sys.executable,str(root/'build_dashboard_data.py')],check=True)
if (app/'src/data.json').exists():
    p=app/'src/data.json'
    snapshot=json.loads(p.read_text(encoding='utf-8'))
    snapshot['buildStatus']='complete'
    p.write_text(json.dumps(snapshot,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    build=os.environ.get('DATA_APP_SCRIPT')
    if build:
        subprocess.run([os.environ.get('NODE', 'node'),build,'build','--project-dir',str(app),'--separate-data'],check=True)
print('Dashboard snapshot refreshed once. No uploads or job mutations. Set DATA_APP_SCRIPT to build the Data app.')
