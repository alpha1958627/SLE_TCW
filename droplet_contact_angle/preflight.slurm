#!/bin/bash
#SBATCH -J hp001YfixCheck
#SBATCH -p C64,C72
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=12G
#SBATCH -t 01:00:00
#SBATCH -o preflight-%j.out
#SBATCH -e preflight-%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
BASE=$PWD
export BASILISK=${BASILISK:?set Basilisk src directory}
export PATH="$BASILISK:$PATH"
export OMP_PROC_BIND=false OMP_PLACES=cores
ROOT="$BASE/runs/preflight_${SLURM_JOB_ID}"
mkdir -p "$ROOT"
python3 tests/check_source.py
sha256sum contact_angle.c tests/*.py > "$ROOT/inputs.sha256"
build_run() {
  local name=$1 kind=$2 transpose=$3 threads=$4 amp=$5 source=$6 end=$7 dt=$8 level=$9
  local dest="$ROOT/$name"
  mkdir "$dest"
  cp "$source" "$dest/contact_angle.c"
  cd "$dest"
  qcc -disable-dimensions -fopenmp -D_GNU_SOURCE -DCOEFF_A=7.233918228e-5 \
    -DHP001_TEST_KIND="$kind" -DHP001_TEST_TRANSPOSE="$transpose" -DHP001_TEST_AMPLITUDE="$amp" \
    -DHP001_OUTPUT_LEVEL=10 -DHP001_ADAPT_EVERY=2 -DLOG_EVERY=1 -DRADIUS_EVERY=1 \
    -Wall -O2 -fno-tree-vectorize contact_angle.c -o solver \
    -L"$BASILISK/gl" -lglutils -lfb_tiny -lm > compile.log 2>&1
  export OMP_NUM_THREADS=$threads
  unset RESTART_FILE
  srun -n 1 ./solver 2026101201 "$dt" "$end" 0 1000000000 "$level" > run.log 2> stderr.log
  printf 'exit_code=0\n' > solver_exit.txt
  python3 "$BASE/tests/audit.py" . --steps "$(python3 -c "print(round($end/$dt))")" --end "$end"
  cd "$BASE"
}
# A control executable differs only in the two active multigrid y offsets.
# It is NEVER a production solver and is generated only inside this job.
python3 - "$ROOT/control.c" <<'PY'
import sys
from pathlib import Path
s = Path('contact_angle.c').read_bytes()
assert s.count(b'M.y[0,1]') == 2
Path(sys.argv[1]).write_bytes(s.replace(b'M.y[0,1]', b'M.y[1]'))
PY
build_run flat 1 0 8 .00001 "$BASE/contact_angle.c" .02 .001 14
build_run sine 2 0 8 .00001 "$BASE/contact_angle.c" .02 .001 14
build_run xy 3 0 1 .003 "$BASE/contact_angle.c" .2 .01 14
build_run yx 3 1 1 .003 "$BASE/contact_angle.c" .2 .01 14
build_run control_xy 3 0 1 .003 "$ROOT/control.c" .2 .01 14
build_run omp1 0 0 1 .00001 "$BASE/contact_angle.c" .12 .001 13
build_run omp8 0 0 8 .00001 "$BASE/contact_angle.c" .12 .001 13
build_run control_drop 0 0 8 .00001 "$ROOT/control.c" .12 .001 13
python3 tests/compare_preflight.py "$ROOT"
echo PASS > "$ROOT/PREFLIGHT_PASS"
