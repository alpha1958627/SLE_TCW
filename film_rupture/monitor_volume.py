"""Read-only volume/startup audit, compatible with SIPF Python 3.6."""
import csv
import datetime as dt
import json
import math
from pathlib import Path
import re
import struct
import subprocess
import sys
import tempfile


def read_log(path):
    rows, marks = [], []
    if not path.exists():
        return rows, marks
    with path.open(errors='replace') as f:
        for line in f:
            match = re.match(r'# Time: Step (\d+): (\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)', line)
            if match:
                marks.append((int(match[1]), dt.datetime.strptime(match[2], '%Y-%m-%d %H:%M:%S')))
            parts = line.split()
            if len(parts) >= 8 and parts[0].isdigit():
                try:
                    rows.append(dict(step=int(parts[0]), hmin=float(parts[2]),
                                     t=float(parts[3]), dt=float(parts[4]), V=float(parts[6])))
                except ValueError:
                    continue
    return rows, marks


def field_volume(path, length):
    data = path.read_bytes()
    if len(data) < 4:
        raise ValueError('partial header')
    n, = struct.unpack('<i', data[:4])
    if n != 64 or len(data) != 4 + 8*n*n:
        raise ValueError('partial or unexpected field')
    h = struct.unpack('<{}d'.format(n*n), data[4:])
    if not all(math.isfinite(x) for x in h):
        raise ValueError('nonfinite field')
    return dict(file=path.name, n=n, V=math.fsum(h)*length*length/(n*n),
                hmin=min(h), hmax=max(h))


def field_time(path):
    return float(re.search(r'_t([0-9.]+)\.(?:bin|hfield)$', path.name)[1])


def audit_run(base, spec):
    case = base / spec['version'] / ('h0_' + spec['h0_ratio']) / ('run_0' + spec['run'])
    out = {k: spec[k] for k in ('task', 'version', 'h0_ratio', 'run', 'seed')}
    out.update(case=str(case), started=(case / '.started').exists(), volume_convergence='INSUFFICIENT_HISTORY')
    length = float(spec['L0'])
    expected = length*length
    out['expected_V_initial'] = expected
    rows, marks = read_log(case / 'log_restart')
    out['log_rows'] = len(rows)
    if rows:
        out['last'] = rows[-1]
        if all(math.isfinite(x[k]) for x in rows for k in ('t', 'V', 'hmin', 'dt')) and rows[0]['V'] > 0:
            v0 = rows[0]['V']
            out['raw_log_V0_over_expected_V0'] = v0 / expected
            out['log_drift_last'] = rows[-1]['V']/v0 - 1
            out['log_drift_max_abs'] = max(abs(x['V']/v0 - 1) for x in rows)
            window = rows[max(0, len(rows) - max(2, len(rows)//4)):]
            elapsed = window[-1]['t'] - window[0]['t']
            out['recent_drift_change'] = (window[-1]['V'] - window[0]['V']) / v0
            out['recent_relative_slope_per_t'] = out['recent_drift_change']/elapsed if elapsed > 0 else None
            out['volume_observation'] = 'DRIFT_ALERT' if out['log_drift_max_abs'] > 1e-3 else 'BELOW_0.1_PERCENT_ALERT'
            # A short flat trace cannot establish long-time or dt/grid convergence.
            stride = max(1, len(rows)//100)
            out['trace'] = rows[::stride]
            if out['trace'][-1] != rows[-1]:
                out['trace'].append(rows[-1])
        else:
            out['volume_observation'] = 'NONFINITE_OR_INVALID_LOG'
    if len(marks) >= 2:
        start, end = marks[max(0, len(marks)-101)], marks[-1]
        seconds = (end[1] - start[1]).total_seconds()
        if seconds > 0 and end[0] > start[0]:
            rate = (end[0] - start[0])/seconds
            remaining = max(0., (int(spec['step_limit']) - end[0])/rate)
            out.update(steps_per_second=rate, seconds_to_step_limit=remaining,
                       estimated_step_limit_finish=(end[1] + dt.timedelta(seconds=remaining)).isoformat(),
                       eta_scope='step_limit_only_not_random_event_time')
    fields = sorted(list(case.glob('hfield_t*.bin')) + list(case.glob('film_t*.hfield')), key=field_time)
    checked, errors = [], []
    for field in fields[:1] + fields[-3:]:
        if field.name in [x['file'] for x in checked]:
            continue
        try:
            result = field_volume(field, length)
            result.update(t=field_time(field), relative_to_expected_initial=result['V']/expected - 1)
            checked.append(result)
        except (ValueError, OSError, struct.error) as exc:
            errors.append(dict(file=field.name, error=str(exc)))
    out.update(field_checks=checked, field_read_errors=errors)
    if checked:
        out['latest_field_drift'] = checked[-1]['relative_to_expected_initial']
    errfile = case / 'stderr.log'
    stderr = ''
    if errfile.exists():
        with errfile.open('rb') as f:
            f.seek(max(0, errfile.stat().st_size - 65536))
            stderr = f.read().decode('utf-8', errors='replace')
    out['stderr_tail_warning_lines'] = [x for x in stderr.splitlines() if re.search(
        r'convergence|warning|nan\b|inf\b|error|floating point', x, re.I)][-8:]
    rupture = re.search(r'RUPTURE hmin=([\d.eE+-]+) t=([\d.eE+-]+) i=(\d+)', stderr)
    if rupture:
        out['event'] = dict(hmin=float(rupture[1]), t=float(rupture[2]), step=int(rupture[3]))
    done = case / 'CONTROL_DONE'
    if done.exists():
        out['control_done'] = done.read_text()
    return out


def self_check():
    with tempfile.TemporaryDirectory(dir=str(Path(__file__).resolve().parent), prefix='volume_check_') as d:
        assert Path(d).resolve().parent == Path(__file__).resolve().parent
        p = Path(d) / 'hfield_t0000000.0.bin'
        p.write_bytes(struct.pack('<i', 64) + struct.pack('<4096d', *([1.]*4096)))
        assert field_volume(p, 100.)['V'] == 10000.
        assert field_time(p) == 0.
        p.write_bytes(b'bad')
        try:
            field_volume(p, 100.)
        except ValueError:
            pass
        else:
            raise AssertionError('accepted partial field')
        p.write_text('i h1 h2 t dt lm V_total V_excess\n0 1 1 0 .003 6 80000 72000\n1000 1 .9 3 .003 6 80008 72008\n')
        rows, marks = read_log(p)
        assert len(rows) == 2 and rows[-1]['step'] == 1000 and not marks
        assert abs(rows[-1]['V']/rows[0]['V'] - 1 - .0001) < 1e-14
    print('volume monitor self-check PASS')


if __name__ == '__main__':
    if '--self-check' in sys.argv:
        self_check()
    else:
        base = Path(__file__).resolve().parent
        specs = list(csv.DictReader((base / 'runs.csv').open(encoding='utf-8-sig')))
        assert len(specs) == 21
        data = dict(time_Asia_Shanghai=dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(),
                    runs=[audit_run(base, spec) for spec in specs])
        jobfile = base / 'SUBMITTED_JOB_ID'
        if jobfile.exists():
            job = jobfile.read_text().strip()
            assert job.isdigit()
            for tool, command in [('squeue', ['squeue', '-r', '-j', job, '-o', '%i|%T|%M|%S|%e|%R']),
                                  ('sacct', ['sacct', '-X', '-j', job, '-P', '-n', '--format=JobIDRaw,State,ExitCode,Elapsed,Start,End,NodeList,AllocCPUS'])]:
                try:
                    data[tool] = subprocess.check_output(command, universal_newlines=True, stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError as exc:
                    data[tool + '_error'] = exc.output
        print(json.dumps(data, indent=2, allow_nan=False))
