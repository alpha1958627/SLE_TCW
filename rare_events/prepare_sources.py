"""Generate the shared kernel from the pinned, corrected film_rupture source.

No legacy AMS solver is imported. The no-repulsion delta is explicit and checked.
Python 3.6 compatible; generated files are deterministic UTF-8/LF.
"""
from pathlib import Path
import argparse
import difflib
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parent
REFERENCE = ROOT.parent / 'film_rupture' / 'source_yfix_h1.05.c'
BASE_SHA = 'e0f7f8103687df7459040b32564bc3f5c74bfed9a95fa0021b992963a8cd0c59'
BASE_COMMIT = 'fed7bb3bf963afe44466017e61199c61613b417d'
CURRENT_COMMIT = '77ce5cb'
CURRENT_SHA = 'dfadb17518a1aecba82efee97d5baaa9a11840c6c327e58afa50416b5c3ddb3d'


def no_repulsion(text):
    old = 'double hstar = pow(0.1, p.dd - 3.);'
    disabled = 'double hstar = 0.0; /* h^-5 repulsion disabled; h^-3 attraction retained. */'
    if disabled in text:
        assert text.count(disabled)==1 and old not in text
        text=text.replace(disabled,old)
    assert text.count(old) == 1
    text = text.replace(old, '/* User model: no h^-5 repulsive contribution. */')
    old = ('beta.x[] = (-p.n * COEFF_A * pow(ff, p.n - p.m - 1.) +\n'
           '                        p.dd * COEFF_A * hstar * pow(ff, p.n - p.dd - 1.)) * mob;')
    assert text.count(old) == 1
    text = text.replace(old, 'beta.x[] = -p.n * COEFF_A * pow(ff, p.n - p.m - 1.) * mob;')
    old = ('p.n * COEFF_A * (p.n - p.m - 1.) * pow(ff, p.n - p.m - 2.) * (f[] - f[-1]) / Delta +\n'
           '                         p.dd * COEFF_A * hstar * (p.n - p.dd - 1.) * pow(ff, p.n - p.dd - 2.) * (f[] - f[-1]) / Delta')
    assert text.count(old) == 1
    text = text.replace(old, 'p.n * COEFF_A * (p.n - p.m - 1.) * pow(ff, p.n - p.m - 2.) * (f[] - f[-1]) / Delta')
    assert 'hstar' not in text and 'p.dd * COEFF_A' not in text
    return text


def generate(check=False):
    digest=hashlib.sha256(REFERENCE.read_bytes()).hexdigest()
    assert digest in (BASE_SHA,CURRENT_SHA), 'Reference changed: review before updating'
    base = REFERENCE.read_text(encoding='utf-8')
    changed = no_repulsion(base)
    first = changed.index(' struct EquationParams {')
    last = changed.index(' /**\n  * @brief 自适应时间步长调整', first)
    solver = changed[first:last]
    first = changed.index(' static uint64_t tcw_mix64')
    last = changed.index('event stability(i++)', first)
    noise = changed[first:last]
    assert 'alpha.y[0,1]' in solver and 'gamma.y[0,1]' in solver
    assert 'f_old[0,-1]' in noise and 'ix = (ix % n + n) % n;' in noise
    header = ('/* Generated from film_rupture at ' + BASE_COMMIT + '.\n'
              ' * Only beta/gamma repulsion removal changes the physical force.\n'
              ' * Original source remains unchanged; see provenance.json. */\n' + solver + noise)
    diff = ''.join(difflib.unified_diff(base.splitlines(True), changed.splitlines(True),
                                      fromfile='film_rupture/reference', tofile='no_repulsion/reference'))
    provenance = dict(reference_commit=CURRENT_COMMIT if digest==CURRENT_SHA else BASE_COMMIT,
                      reference_file=REFERENCE.relative_to(ROOT.parent).as_posix(),reference_sha256=digest,
                      originally_gated_reference_commit=BASE_COMMIT,originally_gated_reference_sha256=BASE_SHA,
                      kernel_byte_identity='unchanged from the originally gated kernel; newer upstream only sets hstar=0',
                      kernel_sha256=hashlib.sha256(header.encode()).hexdigest(),
                      model='mu=-lap(h)+A/h^3; M=h^3+ell*h^2; no h^-5',
                      numerical_changes=['fixed multigrid backend', 'Jacobi relaxation',
                                         'explicit completed-step clock', 'fixed requested dt'],
                      physical_delta='remove both repulsive beta/gamma terms as explicitly requested')
    outputs = {'physics.h': header, 'no_repulsion.diff': diff,
               'provenance.json': json.dumps(provenance, indent=2, sort_keys=True) + '\n'}
    for name, value in outputs.items():
        if check:
            assert (ROOT/name).read_bytes() == value.encode(), name + ' is not reproducible'
        else:
            (ROOT/name).write_bytes(value.encode())
    # An independently event-driven full reference, not the segmented engine.
    native = changed.replace(' #include "grid/tree.h"', '')
    native = '#include "grid/multigrid.h"\n' + native
    old = '     tcw_noise_seed = seed;'
    assert native.count(old) == 1
    native = native.replace(old, old + '\n    double ref_ratio = argc > 6 ? strtod(argv[6], NULL) : 1.01;\n'
        '    tcw_A = 0.0013159/(ref_ratio*ref_ratio);\n    tcw_varphi = 0.01/(ref_ratio*ref_ratio);')
    native = native.replace('L0 = 95.238095238095238;', 'L0 = 100./ref_ratio;')
    native, count = re.subn(r'event stability\(i\+\+\) \{.*?\n\}',
                            'event stability(i++) { dt = dtnext(tcw_dt);\n}', native, count=1, flags=re.S)
    assert count == 1
    build = ROOT/'build'
    build.mkdir(exist_ok=True)
    (build/'reference.c').write_bytes(native.encode())
    return provenance


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--check', action='store_true')
    print(json.dumps(generate(p.parse_args().check), indent=2))
