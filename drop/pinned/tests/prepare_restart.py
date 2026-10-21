"""Prepare two read-only checkpoint controls in new empty directories."""
import hashlib
import json
import argparse
from pathlib import Path
import shutil
import sys
from audit import meta, rows
from native_dump import parse_dump

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('origin', type=Path)
parser.add_argument('root', type=Path)
parser.add_argument('--level', type=int, default=13)
parser.add_argument('--start', type=int, default=400)
parser.add_argument('--tag', default='final')
args = parser.parse_args()
origin, root = args.origin, args.root
assert (args.level, args.start, args.tag) in ((13, 400, 'final'), (14, 100, 't0100'))
checkpoint_name = 'native_final' if args.tag == 'final' else 'checkpoint_'+args.tag
checkpoint = origin/checkpoint_name
native = parse_dump(checkpoint.read_bytes())
assert native['step'] == args.start*100 and abs(native['time']-args.start) < 1e-8
info = meta(origin/('snapshot_'+args.tag+'.meta'))
assert abs(float(info['k'])-37.47688519) < 1e-7
assert float(info['hstar']) == .01 and float(info['noise']) == 0
assert int(info['level']) == args.level
files = [checkpoint_name, 'snapshot_'+args.tag+'.meta', 'volume_state.txt', 'contact_angle_source.c',
         'volume.csv', 'amr_native.csv']
hashes = {f: hashlib.sha256((origin/f).read_bytes()).hexdigest() for f in files}
assert hashes['contact_angle_source.c'] == 'd80bf828a51dd44870d867d45f617d8f50092cc37ba2dae435c718159d53afde'
source = (Path(__file__).parents[1]/'contact_angle.c').read_bytes()
assert source.count(b'M.y[0,1]') == 2
state = meta(origin/'volume_state.txt')
v0 = float(state['V0'])
# Inherit every available original segment's logged tail, including after this CP.
max_drift = max([float(state['max_rel_volume_drift'])] +
                [float(r['max_rel_drift']) for r in rows(origin/'volume.csv')] +
                [abs(float(r['V_after'])-v0)/v0 for r in rows(origin/'amr_native.csv')])
assert max_drift < 1e-8
for label in ('fixed', 'control'):
    dest = root/label
    dest.mkdir()
    shutil.copyfile(checkpoint, dest/'checkpoint_input')
    (dest/'volume_state.txt').write_text('V0=%.17g\nmax_rel_volume_drift=%.17g\n' % (v0, max_drift))
    (dest/'contact_angle.c').write_bytes(source if label == 'fixed' else source.replace(b'M.y[0,1]', b'M.y[1]'))
    assert hashlib.sha256((dest/'checkpoint_input').read_bytes()).hexdigest() == hashes[checkpoint_name]
(root/'origin_hashes.json').write_text(json.dumps(hashes, indent=2)+'\n')
