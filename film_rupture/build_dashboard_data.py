"""Create reviewed dashboard rows from the 21 retained simulation records."""
import csv
import datetime as dt
import json
from pathlib import Path
import math
import os

root = Path(__file__).resolve().parent
ledger = json.loads((root / 'deployment.json').read_text())
inputs = {s: json.loads((root / 'audits' / (s + '_final.json')).read_text()) for s in ledger['sites']}
status, history, fields = [], [], []
for spec in ledger['tasks']:
    site = spec['server']
    record = next(r for r in inputs[site]['runs'] if r['task'] == spec['task'])
    last = record.get('last', {})
    assert record['started'] and last.get('step', 0) > 0, spec['task']
    event = record.get('event')
    endpoint = event if event else last
    common = dict(task=int(spec['task']), ratio=spec['h0_ratio'], version=spec['version'],
                  run=int(spec['run']), seed=int(spec['seed']), server=site,
                  curve=('旧索引 r1' if spec['version'] == 'yold' else 'run ' + spec['run']))
    status.append(dict(common, job=spec['job_id'], state=record.get('scheduler_state','UNKNOWN'), step=endpoint['step'],
        t=endpoint['t'], hmin=endpoint['hmin'], stopReason='hmin≤0.01' if event else '未记录阈值事件',
        logDriftPct=100*record['log_drift_last'],
        maxLogDriftPct=100*record['log_drift_max_abs'], fieldDriftPct=100*record['latest_field_drift'],
        eta=None if event else record.get('estimated_step_limit_finish'), sourceTime=inputs[site]['time_Asia_Shanghai']))
    v0 = record['trace'][0]['V']
    for row in record['trace']:
        assert all(math.isfinite(row[k]) for k in ['t', 'hmin', 'V'])
        history.append(dict(common, t=row['t'], hmin=row['hmin'], step=row['step'],
                            volumeRaw=row['V'], volumeDriftPct=100*(row['V']/v0 - 1)))
    if event and event['t'] > record['trace'][-1]['t']:
        assert event['hmin'] <= .01 and event['step'] > record['trace'][-1]['step']
        history.append(dict(common, t=event['t'], hmin=event['hmin'], step=event['step'],
                            volumeRaw=None, volumeDriftPct=None, endpointSource='stderr RUPTURE'))
    for row in record['field_checks']:
        fields.append(dict(common, t=row['t'], volume=row['V'], driftPct=100*row['relative_to_expected_initial'], file=row['file']))
assert len(status) == 21 and len({r['task'] for r in status}) == 21
assert len([r for r in status if r['version'] == 'yold']) == 1
for task in range(21):
    selected = [r for r in history if r['task'] == task]
    assert selected and all(a['t'] < b['t'] for a, b in zip(selected, selected[1:]))

def source(files, definitions, caveats):
    return dict(provider='SSH / Slurm / Basilisk raw outputs', files=files,
        evidenceFlow=[dict(title='Retained trajectories', detail='deployment.json selects exactly 21 winning task handles; canceled cross-site duplicates are excluded.'),
                      dict(title='Source read', detail='Conservative-v2 cohort; SSH reads existing monitor/logs. Root: ' + ledger['base']),
                      dict(title='Current scope', detail='Fixed simulation cohort, Asia/Shanghai. Latest successful site snapshots are stored as audits/<site>_final.json for the sites in deployment.json.')],
        metricDefinitions=definitions, caveats=caveats)

snapshot = dict(surface='dashboard', title='ls = 0.2 | hmin - t', buildStatus='creating', status='observed', solverVersion='conservative-v2',
    generatedAt=min(x['time_Asia_Shanghai'] for x in inputs.values()), filters=[], queries={
    'film_history': dict(rows=history, source=source(['audits/'+s+'_final.json' for s in ledger['sites']]+['deployment.json'],
        [dict(label='hmin', definition='Global minimum of normalized h, logged every 1000 steps; exact RUPTURE hmin/t/step adds an event endpoint. No volume is imputed at that endpoint.'),
         dict(label='t', definition='Raw code time for each rescaled case. No cross-height physical-time conversion or time alignment is applied.'),
         dict(label='Volume drift (%)', definition='100 * (raw log V / first raw log V - 1); a constant MPI factor cancels, but this is not an independent absolute-volume measurement.')],
        ['Curves are observations, not fitted trajectories. Every curve ends at its own last recorded time.', 'Threshold events are numerical observations; volume drift prevents physical lifetime acceptance. Y-index control uses the same seed as h1.05 run1; both versions include periodic-noise and volume-reduction fixes.', 'Volume drift is observed; mass conservation, dt/grid convergence, and physical acceptance remain unverified.']),
        methods=[dict(language='python', code='build_dashboard_data.py: retain deployment task mapping; emit t and hmin unchanged from each trace; compute 100*(V/V_first-1).')]),
    'run_status': dict(rows=status, source=source(['deployment.json']+['audits/'+s+'_scheduler.txt' for s in ledger['sites']],
        [dict(label='Running', definition='Slurm RUNNING for the retained task handle at the recorded source snapshot.'),
         dict(label='ETA', definition='Extrapolation of measured steps/s to step limit 200000000, not a random rupture-time forecast.')],
        ['State and progress describe source snapshot time, not continuous live polling.', 'Five seeds are reused across heights for paired comparisons.'])),
    'field_volume': dict(rows=fields, source=source(['audits/'+s+'_final.json' for s in ledger['sites']],
        [dict(label='Volume', definition='Independent sum of all 64*64 saved cell heights multiplied by (L0/64)^2.'),
         dict(label='Field drift (%)', definition='100*(field volume/L0^2-1), using the specified initial flat film as reference.')],
        ['Conservative-v2: log V uses a single global reduction. Both yfix and yold controls have periodic-wrapped face-noise indices. Short tests passed; long-time/dt/grid convergence is not established.', 'Saved-field times can lag the latest log. No field normalization or pressure/time shifts are applied.']))})
(root / 'dashboard_snapshot.json').write_text(json.dumps(snapshot, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
app = Path(os.environ.get('DASHBOARD_APP', str(root/'dashboard')))
if (app / 'src/data.json').exists():
    prior = json.loads((app / 'src/data.json').read_text(encoding='utf-8'))
    for key in ['artifactId', 'id', 'theme', 'filters']:
        if key in prior:
            snapshot[key] = prior[key]
    snapshot['buildStatus'] = 'updating'
    (app / 'src/data.json').write_text(json.dumps(snapshot, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(json.dumps(dict(runs=len(status), points=len(history), field_checks=len(fields), snapshot=str(root / 'dashboard_snapshot.json')), ensure_ascii=False))
