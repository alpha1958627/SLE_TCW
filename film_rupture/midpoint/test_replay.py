"""Local reference/native check. Does not substitute for the qcc/MPI gates."""
import argparse
import json
from pathlib import Path
import tempfile
import subprocess
import sys
import numpy as np
from adaptive_runner import BrownianTape, runtime
from accelerated_runner import Trajectory
from native_midpoint_candidate import build


def check(library):
    tape = BrownianTape(2026102201, (2, 8, 8), base_dt=.064, depth=18)
    a, b, c = 0, tape.ticks(.03), tape.ticks(.3)
    total = tape.increment(a, c)
    error = float(np.max(abs(total - tape.increment(a, b) - tape.increment(b, c))))
    assert error < 1e-14
    try:
        BrownianTape(-1, (2, 8, 8))
    except AssertionError:
        pass
    else:
        raise AssertionError('Negative seed was accepted')
    with tempfile.TemporaryDirectory() as folder:
        folder = Path(folder)
        options = dict(dt_cap=.3, base_dt=.064, depth=18)
        ref = Trajectory(1.1, 8, 2026102201, kernel='reference', **options)
        fast = Trajectory(1.1, 8, 2026102201, library=library, **options)
        maximum = 0.
        resumed = None
        for step in range(20):
            ref.advance()
            fast.advance()
            assert ref.status == fast.status == 'RUNNING' and ref.ticks == fast.ticks
            maximum = max(maximum, float(np.max(abs(ref.h - fast.h))))
            if step == 6:
                fast.save(folder / 'restart.npz')
                resumed = Trajectory.load(folder / 'restart.npz', library=library)
            elif step > 6:
                resumed.advance()
                assert np.array_equal(resumed.h, fast.h)
                assert resumed.metadata() == fast.metadata()
        assert maximum < 1e-10
        assert fast.ledger['max_residual'] <= 1e-10
        assert abs(float(fast.h.mean()) - 1.) < 1e-12
        return dict(PASS=True, runtime=runtime(), steps=20, n=8,
                    Brownian_additivity_error=error, invalid_seed_rejected=True,
                    reference_native_max_abs=maximum, restart_bitwise=True,
                    max_residual=fast.ledger['max_residual'],
                    volume_mean_error=abs(float(fast.h.mean()) - 1.))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--library', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.library:
        result = check(args.library)
    else:
        # The child exits before cleanup: Windows cannot delete a loaded DLL.
        with tempfile.TemporaryDirectory() as folder:
            library = build(Path(folder) / 'compiled')
            output = subprocess.check_output([sys.executable, str(Path(__file__).resolve()),
                                              '--library', str(library)], universal_newlines=True)
            result = json.loads(output)
    print(json.dumps(result, indent=2))
    if args.out:
        args.out.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
