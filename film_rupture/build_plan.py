"""Build and statically verify the approval package; never runs a solver."""
from pathlib import Path
from decimal import Decimal
import csv
import difflib
import hashlib
import json
import re

root = Path(__file__).resolve().parent
old_bytes = (root / 'reference/sipf_old_h105.c').read_bytes()
assert hashlib.sha256(old_bytes).hexdigest() == '97a44c89e0b2137c47791bae7fd94532c0a6011349083de2d8090ad0a28e47d1'
old = old_bytes.decode('utf-8')
fixed = old
for before, after, expected in [
    ('alpha.y[1]', 'alpha.y[0,1]', 2),
    ('gamma.y[1]', 'gamma.y[0,1]', 2),
    ('foreach_face(y) {\n        double ff = (f_old[] + f_old[-1])/2.;',
     'foreach_face(y) {\n        double ff = (f_old[] + f_old[0,-1])/2.;', 1),
]:
    assert fixed.count(before) == expected, (before, fixed.count(before))
    fixed = fixed.replace(before, after)
assert len([1 for a, b in zip(old.splitlines(), fixed.splitlines()) if a != b]) == 5
assert len(old.splitlines()) == len(fixed.splitlines())

def replace_one(text, pattern, replacement):
    text, count = re.subn(pattern, replacement, text)
    assert count == 1, (pattern, count)
    return text

rows = []
for height in ['1.05', '1.10', '1.15', '1.20']:
    r = Decimal(height)
    values = [format(Decimal('0.0013159') / r**2, '.17g'),
              format(Decimal('0.01') / r**2, '.17g'),
              format(Decimal(100) / r, '.17g'),
              format(Decimal('0.2') / r, '.17g')]
    a, phi, length, slip = values
    source = fixed
    # Keep the matched h105 parameter literals byte-for-byte identical.
    if height != '1.05':
        source = replace_one(source, r'double tcw_A = [-+0-9.eE]+;', f'double tcw_A = {a};')
        source = replace_one(source, r'double tcw_varphi = [-+0-9.eE]+;', f'double tcw_varphi = {phi};')
        source = replace_one(source, r'(?m)^     L0 = [-+0-9.eE]+;', f'     L0 = {length};')
    name = f'source_yfix_h{height}.c'
    (root / name).write_bytes(source.encode('utf-8'))
    for run in range(1, 6):
        rows.append(dict(task=len(rows), version='yfix', h0_ratio=height, run=run,
                         seed=2026101800 + run, A=a, varphi=phi, L0=length,
                         slip_code=slip, initial_h=1, hstar='0.01', level=6,
                         dt='0.003', step_limit=200000000, mpi_ranks=8,
                         memory_GB=8, wall_limit_days=14, source=name))
control = dict(rows[0], task=20, version='old', source='reference/sipf_old_h105.c')
rows.append(control)
assert len(rows) == 21
assert len({row['seed'] for row in rows}) == 5
assert len({(row['version'], row['h0_ratio'], row['run']) for row in rows}) == 21
assert all(control[k] == rows[0][k] for k in control if k not in {'task', 'version', 'source'})
for height in ['1.05', '1.10', '1.15', '1.20']:
    selected = [row for row in rows if row['h0_ratio'] == height and row['version'] == 'yfix']
    assert [row['seed'] for row in selected] == list(range(2026101801, 2026101806))
with (root / 'runs.csv').open('w', encoding='utf-8-sig', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
(root / 'y_index_fix.diff').write_text(''.join(difflib.unified_diff(
    old.splitlines(True), fixed.splitlines(True), fromfile='sipf_old_h105.c',
    tofile='source_yfix_h1.05.c')), encoding='utf-8')
hashes = {str(p.relative_to(root)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest()
          for p in sorted(root.rglob('*.c'))}
result = dict(status='STATIC_CHECK_PASS', simulations_run=0, submitted_jobs=0,
              github_commit='b5898c3a3a4eb8a350264c21dd004dd67a6d346e',
              source_difference_lines=5, total_runs=21, distinct_seeds=5,
              compiler_checked=False, numerical_validation=False, sha256=hashes)
(root / 'verification.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print(json.dumps({k: v for k, v in result.items() if k != 'sha256'}, indent=2))
for row in rows[::5]:
    print(row['task'], row['version'], row['h0_ratio'], row['A'], row['varphi'], row['L0'], row['slip_code'])

with (root / 'tasks.tsv').open('w', encoding='utf-8', newline='') as f:
    writer = csv.writer(f, delimiter='\t', lineterminator='\n')
    writer.writerow(['task', 'version', 'h0_ratio', 'run', 'seed', 'slip_code', 'source', 'source_sha256'])
    for row in rows:
        writer.writerow([row[k] for k in ['task', 'version', 'h0_ratio', 'run', 'seed', 'slip_code', 'source']] + [hashes[row['source']]])
(root / 'sources.sha256').write_text(''.join(value + '  ' + name + '\n' for name, value in hashes.items()), encoding='utf-8')
