"""Create an isolated q17 observer/controller variant from the pinned q3 source.
The generated physical kernel must remain byte-identical. Python3.6/stdlib.
"""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys

EXPECTED_ENGINE='11653f58f068f265e3f927a2b509cfd168cf56dc65a1750026d2b909973d529d'
EXPECTED_PHYSICS='1f108bce15a106991904c0d4744ef142b4fd6fc6bbf50e01796635aaa1d46290'
MEASURE=r'''static double measure (double *minimum, double *maximum, double *mass) {
  double lo=HUGE_VAL, hi=-HUGE_VAL, vol=0.;
  foreach(reduction(min:lo) reduction(max:hi) reduction(+:vol)) {
    lo=min(lo,h[]); hi=max(hi,h[]); vol+=h[]*sq(Delta);
  }
  *minimum=lo; *maximum=hi; *mass=vol;
  if (lo<=event_h || lo==hi) return lo;
  /* Observer only: periodic global 17x17 mean, with no physical force changes.
     Global field avoids assuming a 9-cell MPI ghost layer exists. */
  double *field=full_field(), box=HUGE_VAL;
  if (pid()==0) {
    double tmp[NN*NN];
    for (int y=0;y<NN;y++) for (int x=0;x<NN;x++) {
      double sum=0.;
      for (int dx=-8;dx<=8;dx++) sum+=field[y*NN+(x+dx+NN)%NN];
      tmp[y*NN+x]=sum/17.;
    }
    for (int y=0;y<NN;y++) for (int x=0;x<NN;x++) {
      double sum=0.;
      for (int dy=-8;dy<=8;dy++) sum+=tmp[((y+dy+NN)%NN)*NN+x];
      box=min(box,sum/17.);
    }
  }
#if _MPI
  MPI_Bcast(&box,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
  free(field); box=max(box,lo);
  return lo+max(0.,min(1.,(lo-event_h)/(.2-event_h)))*(box-lo);
}
'''
def prepare(base,output):
    base=Path(base).resolve();out=Path(output).resolve();here=Path(__file__).resolve().parent
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    assert digest(base/'rare_events/engine.c')==EXPECTED_ENGINE
    assert digest(base/'rare_events/physics.h')==EXPECTED_PHYSICS
    out.mkdir(parents=True,exist_ok=False)
    files=list((base/'rare_events').glob('*.py'))
    files += [base/'rare_events'/n for n in ('engine.c','physics.h','build.sh')]
    files += [base/'film_rupture/source_yfix_h1.05.c']
    original={}
    for p in files:
        name=p.relative_to(base);dest=out/name;dest.parent.mkdir(parents=True,exist_ok=True)
        original[name.as_posix()]=digest(p);data=p.read_bytes()
        if p.suffix=='.py':data=data.decode('utf-8-sig').replace('q3','q17').replace('RESTFE26','RESTFW17').encode()
        dest.write_bytes(data)
    engine=(out/'rare_events/engine.c').read_text();first=engine.index('static double measure (');last=engine.index('\nint main (',first)
    engine=(engine[:first]+MEASURE+engine[last:]).replace('q3','q17').replace('RESTFE26','RESTFW17')
    engine=engine.replace("{'R','E','S','T','F','E','2','6'}","{'R','E','S','T','F','W','1','7'}")
    (out/'rare_events/engine.c').write_bytes(engine.encode())
    subprocess.check_call([sys.executable,str(out/'rare_events/prepare_sources.py')])
    assert digest(out/'rare_events/physics.h')==EXPECTED_PHYSICS
    shutil.copy2(here/'validate_q17.py',out/'validate_q17.py')
    identity=dict(status='prepared',score='q17',checkpoint_abi='RESTFW17',base_engine_sha256=EXPECTED_ENGINE,
       unchanged_physics_sha256=EXPECTED_PHYSICS,base_sources=original,
       generated_sources={p.relative_to(out).as_posix():digest(p) for p in out.rglob('*') if p.is_file()},
       scope='score-only isolated source; no production seeds/states copied and no dynamics changed')
    (out/'variant_identity.json').write_text(json.dumps(identity,indent=2)+'\n')
    print(json.dumps(dict(status='prepared',output=str(out),physics_unchanged=True)))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',required=True);p.add_argument('--output',required=True);a=p.parse_args();prepare(a.base,a.output)
