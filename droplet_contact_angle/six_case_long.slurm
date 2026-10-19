#!/bin/bash
#SBATCH -J hp001YfixLong
#SBATCH -p C64,C72
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=12G
#SBATCH -t 48:00:00
#SBATCH --array=0-5%6
#SBATCH -o long-%A_%a.out
#SBATCH -e long-%A_%a.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
BASE=$PWD
SOURCE_SHA=7ebff6509708e96af9afd0d9bf4eb9cfe636db827b7fc64cc9b2ff25fb8982e3
test "$(sha256sum contact_angle.c | cut -d' ' -f1)" = "$SOURCE_SHA"
SMOKE_ROOT=${SMOKE_ROOT:?completed six-case smoke root}
# Refuse to start any long member until ALL six identities passed.
for label in 15 25 35; do
 for lev in 13 14; do
  test -f "$SMOKE_ROOT/k${label}_L${lev}/SMOKE_PASS"
  python3 - "$SMOKE_ROOT/k${label}_L${lev}/verification.json" "$SOURCE_SHA" <<'PY'
import json, sys
from pathlib import Path
r = json.loads(Path(sys.argv[1]).read_text())
assert r['status'] == 'PASS' and r['source_sha256'] == sys.argv[2]
assert r['real_coupled_steps'] == 10 and r['final_step'] == 10
assert abs(r['final_time']-.1) < 1e-8
PY
 done
done
export BASILISK=${BASILISK:?set Basilisk src directory}
export PATH="$BASILISK:$PATH"
export OMP_NUM_THREADS=8 OMP_PROC_BIND=false OMP_PLACES=cores
index=${SLURM_ARRAY_TASK_ID:?}
angle=$((15+10*(index/2)))
level=$((13+index%2))
case "$angle" in
 15) A=1.362966948e-5;;
 25) A=3.747688519e-5;;
 35) A=7.233918228e-5;;
esac
dest="$BASE/runs/long_${SLURM_ARRAY_JOB_ID}/k${angle}_L${level}"
mkdir -p "$(dirname "$dest")"
mkdir "$dest"
cp contact_angle.c "$dest/"
cd "$dest"
python3 - "$A" "$angle" "$level" <<'PY'
import hashlib, sys
from pathlib import Path
params = dict(A=float(sys.argv[1]), k=float(sys.argv[1])/.01**3,
              nominal_angle_label=int(sys.argv[2]), level_max=int(sys.argv[3]),
              dt=.01, end_time=400, hstar=.01, noise=0, gamma=1, adapt_every=2,
              source_sha256=hashlib.sha256(Path('contact_angle.c').read_bytes()).hexdigest())
Path('input.meta').write_text(''.join('%s=%s\n' % pair for pair in params.items()))
PY
qcc -disable-dimensions -fopenmp -D_GNU_SOURCE -DCOEFF_A="$A" \
 -DHP001_ADAPT_EVERY=2 -DLOG_EVERY=100 -DRADIUS_EVERY=100 \
 -Wall -O2 -fno-tree-vectorize contact_angle.c -o solver \
 -L"$BASILISK/gl" -lglutils -lfb_tiny -lm > compile.log 2>&1
unset RESTART_FILE
# Fresh common initial field; no old run history or normalized volume is used.
set +e
srun -n 1 ./solver 2026101201 .01 400 0 1000000000 "$level" > run.log 2> stderr.log
rc=$?
set -e
printf 'exit_code=%s\n' "$rc" > solver_exit.txt
(( rc == 0 )) || exit "$rc"
python3 "$BASE/tests/audit.py" . --steps 40000 --end 400
echo PASS > NUMERICAL_ENDPOINT_PASS
# This marker does NOT assert a radius plateau or a stable contact angle.
