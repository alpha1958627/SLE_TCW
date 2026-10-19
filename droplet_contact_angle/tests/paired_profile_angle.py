"""Paired estimator sensitivity on equal-time frames; NOT a steady-state gate."""
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[1]/'analysis'))
from angle_methods import read_profile, parabola_angle


def angles(case):
    s, hx, hy = read_profile(case/'profile_final.dat')
    peak = max(max(hx), max(hy))
    results = []
    for values in (hx, hy):
        for sign in (-1, 1):
            points = [(abs(x), h) for x, h in zip(s, values)
                      if sign*x > 0 and .01+.4*(peak-.01) < h < .01+.8*(peak-.01)]
            result = parabola_angle(points, .01)
            assert result['status'] == 'ok'
            results.append(result['geometric_angle_deg'])
    return results


root = Path(sys.argv[1])
fixed, control = (angles(root/p) for p in ('fixed', 'control'))
result = dict(stable_angle=False, purpose='paired estimator sensitivity, one equal-time frame',
              body_window=[.4, .8], fixed_four_directions_deg=fixed,
              control_four_directions_deg=control, fixed_mean_deg=math.fsum(fixed)/4,
              control_mean_deg=math.fsum(control)/4,
              max_paired_angle_difference_deg=max(abs(a-b) for a,b in zip(fixed,control)))
print(json.dumps(result, indent=2))
(root/'paired_angle_sensitivity.json').write_text(json.dumps(result, indent=2)+'\n')
