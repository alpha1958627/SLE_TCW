#!/bin/bash
#SBATCH -J hp001YfixSix
#SBATCH -p C64,C72
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=12G
#SBATCH -t 01:00:00
#SBATCH --array=0-5%6
#SBATCH -o six-%A_%a.out
#SBATCH -e six-%A_%a.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
BASE=$PWD
test -f "${PREFLIGHT_ROOT:?}/PREFLIGHT_PASS"
python3 - "$PREFLIGHT_ROOT/summary.json" <<'PY'
import hashlib, json, sys
from pathlib import Path
r = json.loads(Path(sys.argv[1]).read_text())
assert r['status'] == 'PASS'
assert r['cases']['omp8']['source_sha256'] == hashlib.sha256(Path('contact_angle.c').read_bytes()).hexdigest()
PY
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
dest="$BASE/runs/six_${SLURM_ARRAY_JOB_ID}/k${angle}_L${level}"
mkdir -p "$(dirname "$dest")"
mkdir "$dest"
cp contact_angle.c "$dest/"
cd "$dest"
python3 - "$A" "$angle" "$level" <<'PY'
import hashlib, sys
from pathlib import Path
params = dict(A=float(sys.argv[1]), k=float(sys.argv[1])/.01**3,
              nominal_angle_label=int(sys.argv[2]), level_max=int(sys.argv[3]),
              dt=.01, end_time=.1, hstar=.01, noise=0, gamma=1, adapt_every=2,
              source_sha256=hashlib.sha256(Path('contact_angle.c').read_bytes()).hexdigest())
Path('input.meta').write_text(''.join('%s=%s\n' % pair for pair in params.items()))
PY
qcc -disable-dimensions -fopenmp -D_GNU_SOURCE -DCOEFF_A="$A" \
 -DHP001_ADAPT_EVERY=2 -DLOG_EVERY=1 -DRADIUS_EVERY=1 \
 -Wall -O2 -fno-tree-vectorize contact_angle.c -o solver \
 -L"$BASILISK/gl" -lglutils -lfb_tiny -lm > compile.log 2>&1
unset RESTART_FILE
srun -n 1 ./solver 2026101201 .01 .1 0 1000000000 "$level" > run.log 2> stderr.log
printf 'exit_code=0\n' > solver_exit.txt
python3 "$BASE/tests/audit.py" . --steps 10 --end .1
echo PASS > SMOKE_PASS
