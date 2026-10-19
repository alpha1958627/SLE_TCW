import json
import math
from pathlib import Path
import sys
from audit import compare, field, rows

root = Path(sys.argv[1])
symmetry = compare(root/'xy', root/'yx', transpose=True)
omp = compare(root/'omp1', root/'omp8')
control = compare(root/'omp8', root/'control_drop')
variable = compare(root/'xy', root/'control_xy')
assert symmetry['max_height_difference'] < 1e-7, symmetry
assert omp['max_height_difference'] < 1e-7, omp
_, flat = field(root/'flat/hfield_final.bin')
assert max(abs(v-.01) for v in flat) < 1e-12
n, v0 = field(root/'sine/hfield_t0000.bin')
_, vf = field(root/'sine/hfield_final.bin')
def amplitude(v):
    return 2*math.fsum(math.fsum(v[x::n])*math.cos(16*math.pi*(x+.5)/n) for x in range(n))/(n*n)
rate = .01**3*((16*math.pi/32)**4+2*72.33918228/.01*(16*math.pi/32)**2)
expected = (1+.001*rate)**-20
actual = amplitude(vf)/amplitude(v0)
assert abs(actual-expected)/(1-expected) < .05
result = dict(status='PASS', stable_angle=False, transpose_symmetry=symmetry, omp_1_vs_8=omp,
              repaired_vs_control_droplet=control, repaired_vs_control_variable_mobility=variable,
              sine_amplitude_ratio=actual, sine_backward_euler_reference=expected,
              cases={p.name: json.loads((p/'verification.json').read_text())
                     for p in sorted(root.iterdir()) if (p/'verification.json').exists()})
(root/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
print(json.dumps(result, indent=2))
