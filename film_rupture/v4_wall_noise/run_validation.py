"""Four-core common-Brownian-path regression and dt comparison."""
from pathlib import Path
from array import array
import os,json,csv,subprocess,time,statistics,math,hashlib
root=Path(__file__).resolve().parent
cores={}
for c in sorted(os.sched_getaffinity(0)):
 p=Path(f'/sys/devices/system/cpu/cpu{c}/topology')
 cores.setdefault(((p/'physical_package_id').read_text().strip(),(p/'core_id').read_text().strip()),c)
cpus=list(cores.values())[:4];assert len(cpus)==4
env=dict(os.environ,OMP_NUM_THREADS='4',OMP_DYNAMIC='FALSE',OMP_PROC_BIND='TRUE',OMP_PLACES=','.join('{'+str(c)+'}' for c in cpus))
(root/'allocation.json').write_text(json.dumps(dict(job=os.getenv('SLURM_JOB_ID'),host=os.uname().nodename,physical_cpus=cpus),indent=2))
rows=[];fields={}
for repeat in range(3):
 for phi in (0,.0001,.001):
  configs=[('legacy',.001),('v4',.001),('v4',.005),('v4',.01)]
  if repeat%2:configs.reverse()
  for version,dt in configs:
   folder=root/'results'/f'{version}_phi{phi:g}_dt{dt:g}_r{repeat}';folder.mkdir(parents=True,exist_ok=False)
   exe=root/('legacy_control' if version=='legacy' else 'layer_v4')
   with (folder/'stdout.log').open('w') as out,(folder/'stderr.log').open('w') as err:
    run=subprocess.run([str(exe),str(phi),str(round(.3/dt)),str(round(.02/dt)),'result.json',str(dt)],cwd=folder,env=env,stdout=out,stderr=err,timeout=180)
   assert run.returncode==0,(folder,run.returncode)
   item=json.loads((folder/'result.json').read_text());history=list(csv.DictReader((folder/'history.csv').open()))
   assert len(history)==round(.32/dt) and 'equation_residual' in history[0]
   item.update(version=version,repeat=repeat,folder=folder.name,field_sha256=hashlib.sha256((folder/'final.bin').read_bytes()).hexdigest(),max_equation_residual=max(abs(float(h['equation_residual'])) for h in history),max_mass_drift=max(abs(float(h['relative_mass'])) for h in history))
   assert all(math.isfinite(float(v)) for h in history for v in h.values())
   assert abs(float(history[-1]['time'])-.32)<1e-10
   if version=='v4':assert item['max_equation_residual']<1e-8 and item['max_mass_drift']<1e-10,item
   a=array('d');a.frombytes((folder/'final.bin').read_bytes());fields[(version,phi,dt,repeat)]=a[2::3]
   assert len(a)==128*128*3
   rows.append(item);(root/'progress.json').write_text(json.dumps(rows,indent=2));print(json.dumps(item),flush=True)
summary=[]
for phi in (0,.0001,.001):
 ref=fields[('v4',phi,.001,0)];mean=sum(ref)/len(ref);rms=math.sqrt(sum((x-mean)**2 for x in ref)/len(ref))
 for version,dt in [('legacy',.001),('v4',.001),('v4',.005),('v4',.01)]:
  group=[r for r in rows if r['phi']==phi and r['version']==version and r['dt']==dt]
  seconds=statistics.median(r['seconds'] for r in group)
  a=fields[(version,phi,dt,0)];diff=[x-y for x,y in zip(a,ref)]
  summary.append(dict(version=version,phi=phi,dt=dt,threads=4,seconds=seconds,ms_per_step=seconds/round(.3/dt)*1000,t_star_per_hour=.3*3600/seconds,
    max_equation_residual=max(r['max_equation_residual'] for r in group),max_mass_drift=max(r['max_mass_drift'] for r in group),max_field_difference=max(map(abs,diff)),field_rms_error_over_fluctuation_rms=math.sqrt(sum(x*x for x in diff)/len(diff))/rms,
    noise_path_checksum=group[0]['noise_fingerprint']))
 for dt in (.005,.01):
  base=next(r for r in summary if r['version']=='v4' and r['phi']==phi and r['dt']==.001)
  cur=next(r for r in summary if r['version']=='v4' and r['phi']==phi and r['dt']==dt)
  assert abs(cur['noise_path_checksum']-base['noise_path_checksum'])<1e-12
with (root/'summary.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
(root/'summary.json').write_text(json.dumps(summary,indent=2))
(root/'DONE.json').write_text(json.dumps(dict(runs=len(rows),completed_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),source_sha256=hashlib.sha256((root/'layer_v4.c').read_bytes()).hexdigest(),v4_residual_gate=True,shared_Brownian_path_gate=True),indent=2))
