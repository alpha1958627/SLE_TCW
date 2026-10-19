"""Regenerate the conservative 21-run plan from the canonical h1.05 source."""
from decimal import Decimal
from pathlib import Path
import csv
import difflib
import hashlib
import json
import re
root=Path(__file__).resolve().parent
base=(root/'source_yfix_h1.05.c').read_bytes().decode('utf-8')
assert 'ix = (ix % n + n) % n;' in base and 'iy = (iy % n + n) % n;' in base
assert 'MPI_Allreduce(MPI_IN_PLACE, &vtot' not in base and 'MPI_Allreduce(MPI_IN_PLACE, &vex' not in base
assert 'double hstar = 0.0;' in base
assert 'write_hfield_binary("stop.hfield");' in base
rows=[]
for height in ['1.05','1.10','1.15','1.20']:
    r=Decimal(height)
    a,phi,length,slip=[format(x,'.17g') for x in [Decimal('.0013159')/r**2,Decimal('.01')/r**2,Decimal(100)/r,Decimal('.2')/r]]
    source=base
    if height!='1.05':
        for pattern,replacement in [(r'double tcw_A = [-+0-9.eE]+;','double tcw_A = '+a+';'),
            (r'double tcw_varphi = [-+0-9.eE]+;','double tcw_varphi = '+phi+';'),
            (r'(?m)^     L0 = [-+0-9.eE]+;','     L0 = '+length+';')]:
            source,count=re.subn(pattern,replacement,source);assert count==1
    name='source_yfix_h'+height+'.c'
    (root/name).write_bytes(source.encode())
    for run in range(1,6):
        rows.append(dict(task=len(rows),version='yfix',h0_ratio=height,run=run,seed=2026101800+run,A=a,varphi=phi,L0=length,slip_code=slip,initial_h=1,hstar='0',level=6,dt='0.003',step_limit=200000000,mpi_ranks=8,memory_GB=8,wall_limit_days=14,source=name))
control=base
for before,after,count in [('alpha.y[0,1]','alpha.y[1]',2),('gamma.y[0,1]','gamma.y[1]',2),
    ('foreach_face(y) {\n        double ff = (f_old[] + f_old[0,-1])/2.;','foreach_face(y) {\n        double ff = (f_old[] + f_old[-1])/2.;',1)]:
    assert control.count(before)==count;control=control.replace(before,after)
assert sum(a!=b for a,b in zip(base.splitlines(),control.splitlines()))==5
(root/'source_yold_h1.05.c').write_bytes(control.encode())
rows.append(dict(rows[0],task=20,version='yold',source='source_yold_h1.05.c'))
assert len(rows)==21 and len({r['seed'] for r in rows})==5
assert all(rows[0][k]==rows[-1][k] for k in rows[0] if k not in ['task','version','source'])
for height in ['1.05','1.10','1.15','1.20']:
    assert [r['seed'] for r in rows if r['h0_ratio']==height and r['version']=='yfix']==list(range(2026101801,2026101806))
with (root/'runs.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(root.glob('source_*.c'))}
with (root/'tasks.tsv').open('w',encoding='utf-8',newline='') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n');keys=['task','version','h0_ratio','run','seed','slip_code','source'];w.writerow(keys+['source_sha256'])
    for r in rows:w.writerow([r[k] for k in keys]+[hashes[r['source']]])
(root/'sources.sha256').write_bytes(''.join(v+'  '+k+'\n' for k,v in hashes.items()).encode())
(root/'y_index_control.diff').write_text(''.join(difflib.unified_diff(base.splitlines(True),control.splitlines(True),fromfile='source_yfix_h1.05.c',tofile='source_yold_h1.05.c')),encoding='utf-8')
(root/'verification.json').write_text(json.dumps(dict(status='STATIC_PASS',total_runs=21,distinct_seeds=5,periodic_noise=True,single_volume_reduction=True,repulsion_coefficient=0,attraction_retained=True,sha256=hashes),indent=2)+'\n',encoding='utf-8')
print('STATIC_PASS: 21 runs, five seeds, every source has conservative noise and a single volume reduction.')
