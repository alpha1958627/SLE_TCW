#!/bin/bash
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")/.."
: "${BASILISK:?Set BASILISK to a compatible Basilisk src directory}"
export CC99="${CC99:-mpicc}"
export PATH="$BASILISK:$PATH"
python3 rare_events/prepare_sources.py --check
SIGNATURE=$(python3 - <<'PY'
import hashlib
from pathlib import Path
p=Path('rare_events')
s=hashlib.sha256((p/'engine.c').read_bytes()+(p/'physics.h').read_bytes()).hexdigest()
(p/'build/signature.h').write_text('#define ENGINE_SIGNATURE "'+s+'"\n')
print(s)
PY
)
COMMON=(-DJACOBI=1 -disable-dimensions -D_GNU_SOURCE -O2 -Wall)
"$BASILISK/qcc" -D_MPI=0 "${COMMON[@]}" \
  rare_events/engine.c -o rare_events/build/serial -lm > rare_events/build/serial.log 2>&1
"$BASILISK/qcc" -D_MPI=1 "${COMMON[@]}" \
  rare_events/engine.c -o rare_events/build/mpi -lm > rare_events/build/mpi.log 2>&1
"$BASILISK/qcc" -D_MPI=0 "${COMMON[@]}" rare_events/build/reference.c \
  -o rare_events/build/reference_serial -L"$BASILISK/gl" -lglutils -lfb_tiny -lm > rare_events/build/reference.log 2>&1
python3 - "$SIGNATURE" <<'PY'
from pathlib import Path
import hashlib,json,sys
p=Path('rare_events'); names=['engine.c','physics.h','build/signature.h','build/serial','build/mpi','build/reference.c','build/reference_serial']
data=dict(engine_signature=sys.argv[1], reference=json.loads((p/'provenance.json').read_text()),
          files={n:hashlib.sha256((p/n).read_bytes()).hexdigest() for n in names})
(p/'build/build.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps(data,indent=2))
PY
