#!/bin/bash
#SBATCH -J hp001YfixRestart
#SBATCH -p C64,C72
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=12G
#SBATCH -t 01:00:00
#SBATCH -o restart-%j.out
#SBATCH -e restart-%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
BASE=$PWD
export BASILISK=${BASILISK:?set Basilisk src directory}
export PATH="$BASILISK:$PATH"
export OMP_NUM_THREADS=8 OMP_PROC_BIND=false OMP_PLACES=cores
ORIGIN=${ORIGIN:?existing k25 checkpoint directory, read only}
RESTART_LEVEL=${RESTART_LEVEL:-13}
RESTART_START=${RESTART_START:-400}
RESTART_TAG=${RESTART_TAG:-final}
END=$(python3 -c "print(int('$RESTART_START')+.1)")
ROOT="$BASE/runs/restart_${SLURM_JOB_ID}"
mkdir -p "$ROOT"
python3 tests/prepare_restart.py "$ORIGIN" "$ROOT" --level "$RESTART_LEVEL" --start "$RESTART_START" --tag "$RESTART_TAG"
for label in fixed control; do
 cd "$ROOT/$label"
 qcc -disable-dimensions -fopenmp -D_GNU_SOURCE -DCOEFF_A=3.747688519e-5 \
  -DHP001_ADAPT_EVERY=2 -DLOG_EVERY=1 -DRADIUS_EVERY=1 \
  -Wall -O2 -fno-tree-vectorize contact_angle.c -o solver \
  -L"$BASILISK/gl" -lglutils -lfb_tiny -lm > compile.log 2>&1
 RESTART_FILE="$PWD/checkpoint_input" srun -n 1 ./solver 2026101201 .01 "$END" 0 1000000000 "$RESTART_LEVEL" > run.log 2> stderr.log
 printf 'exit_code=0\n' > solver_exit.txt
 python3 "$BASE/tests/audit.py" . --steps 10 --end "$END"
done
cd "$BASE"
python3 tests/compare_restart.py "$ROOT" "$ORIGIN"
echo PASS > "$ROOT/RESTART_COMPARE_PASS"
