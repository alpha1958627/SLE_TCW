"""Finite-window direct MC example; no scheduler, resampling, or implicit resume."""
import argparse
import csv
import json
import math
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--kernel', choices=['reference', 'native', 'basilisk', 'basilisk_fvsp'], default='reference')
    p.add_argument('--height', type=float, default=1.1)
    p.add_argument('--seed', type=int, default=2026102201)
    p.add_argument('--cap', type=float, default=.3)
    p.add_argument('--time', type=float, default=30.)
    p.add_argument('--library', type=Path)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    for name in ('height', 'cap', 'time'):
        if not math.isfinite(getattr(a, name)) or getattr(a, name) <= 0:
            p.error('--%s must be finite and positive' % name)
    if not 0 <= a.seed < 2**64:
        p.error('--seed must be an unsigned 64-bit integer')
    if a.kernel == 'native' and a.library is None:
        p.error('--library is required for native')
    if a.kernel == 'reference' and a.library is not None:
        p.error('--library does not apply to reference')
    # Same protocol as the frozen production core, with an explicit finite end.
    options = dict(dt_cap=a.cap, base_dt=.064, depth=18)
    if a.kernel in ('basilisk', 'basilisk_fvsp'):
        if a.kernel == 'basilisk_fvsp':
            from basilisk_fvsp_candidate import Trajectory
        else:
            from basilisk_midpoint_candidate import Trajectory
        trajectory = Trajectory(a.height, 64, a.seed, library=a.library, **options)
        rank = trajectory.film.rank
    else:
        from accelerated_runner import Trajectory
        trajectory = Trajectory(a.height, 64, a.seed, kernel=a.kernel,
                                library=a.library, **options)
        rank = 0
    end = trajectory.noise.ticks(a.time)
    if end < 1 or abs(end * trajectory.noise.quantum - a.time) > 1e-12:
        p.error('--time must be representable on the Brownian tick grid')
    # Every MPI rank checks the path before rank 0 creates it.
    if a.out.exists():
        p.error('--out must not already exist')
    if a.kernel in ('basilisk', 'basilisk_fvsp'):
        trajectory.film.library.bm_sync()
    if rank == 0:
        a.out.mkdir(parents=True, exist_ok=False)
    rows = [dict(time=0., hmin=float(trajectory.h.min()), dt=0., residual=0.)]
    while trajectory.status == 'RUNNING' and trajectory.ticks < end:
        record = trajectory.advance(stop_ticks=end)
        rows.append({key: record.get(key) for key in ('time', 'hmin', 'dt', 'residual')})
    state = trajectory.status if trajectory.status != 'RUNNING' else 'RIGHT_CENSORED_AT_FIXED_HORIZON'
    if rank == 0:
        trajectory.save(a.out / 'checkpoint.npz')
        summary = dict(outcome=state, time=trajectory.time, event=trajectory.event,
                       metadata=trajectory.metadata(), independent_mean=float(trajectory.h.mean()),
                       finite_positive=bool((trajectory.h > 0).all()))
        (a.out / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
        with (a.out / 'hmin.csv').open('w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=['time', 'hmin', 'dt', 'residual'])
            writer.writeheader()
            writer.writerows(rows)
        print(json.dumps(dict(outcome=state, time=trajectory.time, event=trajectory.event)))
    return 2 if state == 'NUMERICAL_FAILURE' else 0


if __name__ == '__main__':
    raise SystemExit(main())
