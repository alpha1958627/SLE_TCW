import hashlib
import json
from pathlib import Path
import sys
from audit import compare

root, origin = map(Path, sys.argv[1:])
expected = json.loads((root/'origin_hashes.json').read_text())
assert expected == {f: hashlib.sha256((origin/f).read_bytes()).hexdigest() for f in expected}
delta = compare(root/'fixed', root/'control')
assert delta['max_height_difference'] < 1e-7, delta
result = dict(status='PASS', stable_angle=False, original_files_unchanged=True,
              repaired_vs_control=delta,
              cases={p: json.loads((root/p/'verification.json').read_text()) for p in ('fixed', 'control')},
              origin_hashes=expected)
(root/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
print(json.dumps(result, indent=2))
