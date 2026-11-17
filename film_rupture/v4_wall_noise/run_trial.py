"""Bounded four-core v4 run after validation/publication; dt stays at .001."""
from pathlib import Path
from array import array
import os,json,subprocess,csv,math,time,hashlib
root=Path(__file__).resolve().parent
cores={}
for cpu in sorted(os.sched_getaffinity(0)):
 p=Path(f'/sys/devices/system/cpu/cpu{cpu}/topology')
 cores.setdefault(((p/'physical_package_id').read_text(),(p/'core_id').read_text()),cpu)
cpus=list(cores.values())[:4];assert len(cpus)==4
folder=root/'trial_4core_t3';folder.mkdir(exist_ok=False)
env=dict(os.environ,OMP_NUM_THREADS='4',OMP_DYNAMIC='FALSE',OMP_PROC_BIND='TRUE',OMP_PLACES=','.join('{'+str(c)+'}' for c in cpus))
(folder/'allocation.json').write_text(json.dumps(dict(job=os.getenv('SLURM_JOB_ID'),host=os.uname().nodename,physical_cpus=cpus),indent=2))
results=[];reference=None
for dt in (.001,.005,.01):
 case=folder/f'dt{dt:g}';case.mkdir()
 with (case/'stdout.log').open('w') as out,(case/'stderr.log').open('w') as err:
  r=subprocess.run([str(root/'layer_v4'),'.001',str(round(3/dt)),'0','result.json',str(dt)],cwd=case,env=env,stdout=out,stderr=err,timeout=600)
 assert r.returncode==0,r.returncode
 rows=list(csv.DictReader((case/'history.csv').open()));assert len(rows)==round(3/dt)
 assert all(math.isfinite(float(v)) for row in rows for v in row.values())
 mass=max(abs(float(row['relative_mass'])) for row in rows);res=max(abs(float(row['equation_residual'])) for row in rows)
 assert mass<1e-10 and res<1e-8 and min(float(row['hmin']) for row in rows)>0
 a=array('d');a.frombytes((case/'final.bin').read_bytes());field=a[2::3]
 if reference is None:reference=field
 mean=sum(reference)/len(reference);rms=math.sqrt(sum((x-mean)**2 for x in reference)/len(reference))
 error=math.sqrt(sum((x-y)**2 for x,y in zip(field,reference))/len(field))/rms
 result=json.loads((case/'result.json').read_text());result.update(completed_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),source_sha256=hashlib.sha256((root/'layer_v4.c').read_bytes()).hexdigest(),field_sha256=hashlib.sha256((case/'final.bin').read_bytes()).hexdigest(),max_mass_drift=mass,max_equation_residual=res,field_rms_error_over_fluctuation_rms=error)
 (case/'DONE.json').write_text(json.dumps(result,indent=2));results.append(result);print(json.dumps(result),flush=True)
assert max(r['noise_fingerprint'] for r in results)-min(r['noise_fingerprint'] for r in results)<1e-11
(folder/'DONE.json').write_text(json.dumps(results,indent=2))
