"""Local CPU checks for the consolidated phi interface (no scheduler/GPU).
Run: python -m unittest discover -s film_rupture/midpoint/tests -v
"""
from pathlib import Path
import json
import sys
import tempfile
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'runtime'))
from midpoint_candidate import Film
from accelerated_runner import Trajectory


class RuntimeChecks(unittest.TestCase):
    def trajectory(self, phi_ref=.01):
        return Trajectory(1.2, 8, 2026102701, dt_cap=.001,
                          base_dt=.064, depth=18, kernel='reference', phi_ref=phi_ref)

    def test_default_and_noise_scaling(self):
        default, explicit = Film(1.2, 8), Film(1.2, 8, phi_ref=.01)
        h = 1 + np.random.default_rng(3).uniform(-.01, .01, (8, 8))
        dw = np.random.default_rng(4).normal(size=(2, 8, 8)) * .01
        np.testing.assert_array_equal(default.stochastic(h, dw), explicit.stochastic(h, dw))
        other = Film(1.2, 8, phi_ref=.013345)
        np.testing.assert_array_equal(default.force(h), other.force(h))
        np.testing.assert_allclose(other.stochastic(h, dw),
                                   default.stochastic(h, dw) * np.sqrt(.013345/.01), atol=1e-17)
        for invalid in (0, -1, float('nan'), float('inf')):
            with self.assertRaises(AssertionError):
                Film(1.2, 8, phi_ref=invalid)

    def test_restart_and_fork_both_phi_values(self):
        for phi in (.01, .013345):
            with self.subTest(phi=phi), tempfile.TemporaryDirectory() as folder:
                full, split = self.trajectory(phi), self.trajectory(phi)
                full.run_steps(8)
                split.run_steps(4)
                path = Path(folder) / 'state.npz'
                split.save(path)
                resumed = Trajectory.load(path)
                self.assertEqual(resumed.film.phi_ref, phi)
                resumed.run_steps(4)
                self.assertEqual(resumed.ticks, full.ticks)
                np.testing.assert_array_equal(resumed.h, full.h)
                child = resumed.fork(2026102702)
                self.assertEqual(child.film.phi_ref, phi)
                self.assertEqual(child.noise.seed, 2026102702)
                np.testing.assert_array_equal(child.h, resumed.h)
                self.assertLess(full.ledger['max_residual'], 1e-10)
                self.assertLess(abs(full.h.mean()-1), 1e-12)
                # A source identity mismatch must fail before any resume.
                with np.load(path, allow_pickle=False) as archive:
                    data = {key: archive[key] for key in archive.files}
                metadata = json.loads(str(data['metadata'].item()))
                metadata['fingerprint']['sources']['midpoint_candidate.py'] = '0'*64
                payload = json.dumps(metadata, sort_keys=True, separators=(',', ':'))
                import hashlib
                data['metadata'] = payload
                data['sha256'] = hashlib.sha256(payload.encode()+data['h'].tobytes()+data['noise_values'].tobytes()).hexdigest()
                with path.open('wb') as stream:
                    np.savez_compressed(stream, **data)
                with self.assertRaisesRegex(AssertionError, 'fingerprint mismatch'):
                    Trajectory.load(path)


if __name__ == '__main__':
    unittest.main()
