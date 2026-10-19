"""Independently recompute C q3 and raw volume from existing completed fields.
No solver is run, no field is normalized. Python 3.6, standard library only.
"""
from pathlib import Path
import argparse
import json
import math
import struct
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from rare_events import core


def q3(field):
    lo=min(field)
    box=min(math.fsum(field[((y+dy)%64)*64+(x+dx)%64] for dy in (-1,0,1) for dx in (-1,0,1))/9
            for y in range(64) for x in range(64))
    return lo+max(0.,min(1.,(lo-.01)/.19))*(box-lo)


def check(root):
    assert q3((1.,)*4096)==1.
    test=[1.]*4096;test[0]=.01;assert q3(test)==.01
    count=0;maxq=0.;maxmass=0.;records=0;maxdrift=0.
    for path in Path(root).rglob('summary.json'):
        final=path.parent/'final.chk'
        if not final.exists():continue
        summary=core.load(path);data=core.checkpoint(final)
        error=abs(q3(data['field'])-summary['q3']);maxq=max(maxq,error);assert error<1e-12
        header=core.HEADER.unpack_from(final.read_bytes());L=header[11]
        mass=math.fsum(data['field'])*(L/64)**2
        maxdrift=max(maxdrift,abs((mass-data['mass0'])/data['mass0']))
        error=abs((mass-data['mass0'])/data['mass0']-summary['mass_relative_drift'])
        maxmass=max(maxmass,error);assert error<1e-12
        count+=1
        table=path.parent/'records.tsv'
        if table.exists():
            for line in table.read_text().splitlines():
                step,t,score,name=line.split('\t');d=core.checkpoint(path.parent/name)
                assert d['step']==int(step) and abs(d['time']-float(t))<1e-8
                error=abs(q3(d['field'])-float(score));maxq=max(maxq,error);assert error<1e-12
                records+=1
    assert count>0,'No completed fields; do not report a fabricated zero error'
    return dict(status='PASS',completed_fields=count,saved_record_fields=records,
                max_q3_error=maxq,max_raw_mass_drift_recheck_error=maxmass,
                max_abs_raw_mass_relative_drift=maxdrift,
                scope='independent saved-field/score/volume audit only, no MPI or lifetime claim')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory');p.add_argument('--output')
    a=p.parse_args();result=check(a.directory)
    if a.output:core.save(a.output,result)
    print(json.dumps(result,indent=2))
