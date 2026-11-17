"""Unattended fixed-step run with atomic checkpoints and deterministic noise state."""
from pathlib import Path
from array import array
import os,json,subprocess,time,csv,math,hashlib,signal
root=Path(__file__).resolve().parent
cfg=json.loads((root/'long_config.json').read_text());stop=False
def request_stop(*args):
 global stop
 stop=True
signal.signal(signal.SIGUSR1,request_stop);signal.signal(signal.SIGTERM,request_stop)
started=time.monotonic();cores={}
for cpu in sorted(os.sched_getaffinity(0)):
 p=Path(f'/sys/devices/system/cpu/cpu{cpu}/topology');cores.setdefault(((p/'physical_package_id').read_text(),(p/'core_id').read_text()),cpu)
cpus=list(cores.values())[:4];assert len(cpus)==4
env=dict(os.environ,OMP_NUM_THREADS='4',OMP_DYNAMIC='FALSE',OMP_PROC_BIND='TRUE',OMP_PLACES=','.join('{'+str(c)+'}' for c in cpus))
out=root/'long_4core';out.mkdir(exist_ok=True)
assert hashlib.sha256((root/'layer_v4_long.c').read_bytes()).hexdigest()==cfg['source_sha256']
def atomic(path,obj):
 temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(obj,indent=2));os.replace(temp,path)
atomic(out/'allocation.json',dict(job=os.getenv('SLURM_JOB_ID'),host=os.uname().nodename,physical_cpus=cpus))
def field(path):
 a=array('d');a.frombytes(path.read_bytes());assert len(a)==128**2*3;return a[2::3]
def launch(folder,steps,resume=None,offset=0,original=False):
 folder.mkdir(exist_ok=False)
 args=[str(root/('layer_v4' if original else 'layer_v4_long')),str(cfg['phi']),str(steps),'0','result.json',str(cfg['dt'])]
 if resume is not None:args += [str(resume),str(offset)]
 with (folder/'stdout.log').open('w') as stdout,(folder/'stderr.log').open('w') as stderr:
  r=subprocess.run(args,cwd=folder,env=env,stdout=stdout,stderr=stderr,timeout=1200)
 assert r.returncode==0,(folder,r.returncode)
 rows=list(csv.DictReader((folder/'history.csv').open()));assert len(rows)==steps
 assert all(math.isfinite(float(v)) for row in rows for v in row.values())
 assert min(float(row['hmin']) for row in rows)>0
 assert max(abs(float(row['equation_residual'])) for row in rows)<cfg['max_equation_residual']
 return rows,json.loads((folder/'result.json').read_text())
try:
 seed=root/'trial_4core_t3/dt0.001/final.bin';seed_result=json.loads((seed.parent/'DONE.json').read_text())
 assert seed_result['source_sha256']==cfg['published_kernel_sha256'] and seed_result['steps']==3000
 assert hashlib.sha256(seed.read_bytes()).hexdigest()==seed_result['field_sha256']
 # Compare original and restart-capable drivers, then 2-step versus split 1+1.
 if not (out/'restart_gate.json').exists():
  launch(out/'gate_original',2,original=True);launch(out/'gate_long',2)
  assert max(abs(a-b) for a,b in zip(field(out/'gate_original/final.bin'),field(out/'gate_long/final.bin')))<1e-10
  launch(out/'gate_two',2,seed,3000);launch(out/'gate_one',1,seed,3000)
  launch(out/'gate_split',1,out/'gate_one/final.bin',3001)
  error=max(abs(a-b) for a,b in zip(field(out/'gate_two/final.bin'),field(out/'gate_split/final.bin')))
  assert error<1e-10,error
  atomic(out/'restart_gate.json',dict(passed=True,max_split_error=error,source_sha256=cfg['source_sha256']))
 if (out/'checkpoint.json').exists():
  cp=json.loads((out/'checkpoint.json').read_text());assert cp['source_sha256']==cfg['source_sha256']
  current=root/cp['field'];assert hashlib.sha256(current.read_bytes()).hexdigest()==cp['field_sha256']
  offset=cp['base_step'];initial_mean=cp['initial_mean']
 else:
  current=seed;offset=3000;h=field(seed);initial_mean=(sum(h)/len(h))/(1+seed_result['relative_mass'])
  cp=dict(status='READY',base_step=offset,t=offset*.001,field=str(seed.relative_to(root)),field_sha256=seed_result['field_sha256'],initial_mean=initial_mean,source_sha256=cfg['source_sha256'])
  atomic(out/'checkpoint.json',cp)
 while offset<round(cfg['target_t']/cfg['dt']) and not stop:
  if time.monotonic()-started>cfg['internal_budget_seconds']:break
  n=min(cfg['chunk_steps'],round(cfg['target_t']/cfg['dt'])-offset);folder=out/f'step_{offset:09d}'
  rows,result=launch(folder,n,current,offset)
  start_field=field(current);start_ratio=(sum(start_field)/len(start_field))/initial_mean
  max_mass=max(abs((1+float(row['relative_mass']))*start_ratio-1) for row in rows)
  assert max_mass<cfg['max_global_mass_drift'],max_mass
  h=field(folder/'final.bin');mass=(sum(h)/len(h))/initial_mean-1
  assert abs(mass)<cfg['max_global_mass_drift'],mass
  offset+=n;current=folder/'final.bin'
  cp=dict(status='RUNNING',base_step=offset,t=offset*.001,field=str(current.relative_to(root)),field_sha256=hashlib.sha256(current.read_bytes()).hexdigest(),initial_mean=initial_mean,source_sha256=cfg['source_sha256'],hmin=min(h),hmax=max(h),global_mass_drift=mass,max_equation_residual=max(abs(float(row['equation_residual'])) for row in rows),chunk_seconds=result['seconds'],t_star_per_hour=result['t_star_per_hour'],updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
  atomic(out/'checkpoint.json',cp);print(json.dumps(cp),flush=True)
 cp['status']='COMPLETED_T400' if offset>=400000 else 'CHECKPOINTED_BUDGET_OR_SIGNAL_STOP'
 atomic(out/'checkpoint.json',cp);atomic(out/'FINAL.json',cp)
except Exception as exc:
 atomic(out/'FAILED.json',dict(error=repr(exc),source_sha256=cfg['source_sha256'],time_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
 raise
