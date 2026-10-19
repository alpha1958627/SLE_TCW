from pathlib import Path
import csv
import json
from monitor_volume import read_log,field_volume
root=Path(__file__).resolve().parent
rows=list(csv.DictReader((root/'runs.csv').open(encoding='utf-8-sig')))
results=[]
for spec in rows:
    case=root/'smoke_v2'/('task_'+spec['task'])
    log,_=read_log(case/'log_restart')
    assert len(log)>=2, spec['task']
    length=float(spec['L0']); expected=length*length
    final=field_volume(case/'stop.hfield',length)
    raw_ratio=log[0]['V']/expected
    log_drift=max(abs(row['V']/log[0]['V']-1) for row in log)
    field_drift=final['V']/expected-1
    assert abs(raw_ratio-1)<1e-9, (spec['task'],'absolute volume',raw_ratio)
    assert log_drift<1e-9 and abs(field_drift)<1e-9, (spec['task'],log_drift,field_drift)
    results.append(dict(task=spec['task'],version=spec['version'],seed=spec['seed'],log_rows=len(log),raw_V0_over_L0_squared=raw_ratio,log_max_drift=log_drift,final_field_drift=field_drift))
assert len(results)==21
data=dict(status='SMOKE_PASS_21',steps_requested=1024,mpi_ranks=8,maximum_log_drift=max(x['log_max_drift'] for x in results),maximum_abs_final_field_drift=max(abs(x['final_field_drift']) for x in results),results=results)
(root/'smoke_verification.json').write_text(json.dumps(data,indent=2)+'\n')
(root/'SMOKE_OK').write_text('all 21 source/seed combinations passed\n')
print(json.dumps(data,indent=2))
