"""Recompute the published angles from the bundled, unchanged profiles."""
from pathlib import Path
import hashlib
import json
import math

from angle_methods import circular_arc_angle, parabola_angle, read_profile

DEFAULT_DATA = Path(__file__).resolve().parents[1] / 'results' / '2026-10-21'


def load_measurements(directory=DEFAULT_DATA):
    directory = Path(directory)
    manifest = json.loads((directory / 'snapshots.json').read_text(encoding='utf-8'))
    expected = json.loads((directory / 'k_theta_data.json').read_text(encoding='utf-8'))
    references = {row['case']: row for row in expected['rows']}
    measurements = []
    for snapshot in manifest['cases']:
        path = directory / snapshot['profile']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == snapshot['profile_sha256']
        metadata = dict(line.split('=', 1) for line in
                        (directory / snapshot['metadata']).read_text().splitlines() if '=' in line)
        assert float(metadata['hstar']) == .01 and float(metadata['noise']) == 0
        assert float(metadata['gamma']) == 1
        assert float(metadata['n_repulsive']) == 5 and float(metadata['m_attractive']) == 3
        assert abs(float(metadata['time']) - snapshot['time']) < 1e-8
        assert abs(float(metadata['k']) - snapshot['k']) < 1e-10
        s, hx, hy = read_profile(path)
        peak = max(max(hx), max(hy))
        fits = []
        for axis, heights in [('x', hx), ('y', hy)]:
            for sign in (-1, 1):
                for lo, hi in ((.2, .6), (.4, .8), (.6, .9)):
                    points = [(abs(x), h) for x, h in zip(s, heights)
                              if sign*x > 0 and .01 + lo*(peak-.01) < h < .01 + hi*(peak-.01)]
                    parabola = parabola_angle(points, .01)
                    circle = circular_arc_angle(points, .01)
                    assert parabola['status'] == circle['status'] == 'ok'
                    fits.append(dict(axis=axis, sign=sign, window=[lo, hi],
                                     parabola=parabola, circle=circle))
        angles = [f['parabola']['geometric_angle_deg'] for f in fits if f['window'] == [.4, .8]]
        theta = math.fsum(angles) / 4
        reference = references[snapshot['case']]
        assert abs(theta - reference['theta_deg']) < 1e-8
        theory = math.degrees(math.atan(math.sqrt(snapshot['k']*.01/2)))
        assert abs(theory - reference['theta_model_deg']) < 1e-10
        measurements.append(dict(snapshot, theta_deg=theta, theta_model_deg=theory,
                                 directional_min=min(angles), directional_max=max(angles),
                                 fits=fits, s=s, hx=hx, hy=hy))
    assert len(measurements) == len(references) == 6
    return measurements


def verify_assets(directory=DEFAULT_DATA):
    directory = Path(directory)
    checksums = json.loads((directory / 'sha256.json').read_text(encoding='utf-8'))
    for name, expected in checksums.items():
        path = (directory / name).resolve()
        assert path.is_relative_to(directory.resolve())
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, name
    provenance = json.loads((directory / 'provenance.json').read_text(encoding='utf-8'))
    source = Path(__file__).resolve().parents[1] / 'contact_angle.c'
    assert hashlib.sha256(source.read_bytes()).hexdigest() == provenance['solver_source_sha256']
    return len(checksums)
