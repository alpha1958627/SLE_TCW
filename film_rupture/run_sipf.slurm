#!/bin/bash
#SBATCH -J ls02C21
#SBATCH -p C64
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -c 1
#SBATCH --array=0-20%21
#SBATCH --mem=8G
#SBATCH -t 14-00:00:00
#SBATCH --no-requeue
#SBATCH -o slurm-%A_%a.out
#SBATCH -e slurm-%A_%a.err
if [[ -f /etc/profile.d/dawning.sh ]]; then source /etc/profile.d/dawning.sh; fi
set -euo pipefail
BASE=${SLURM_SUBMIT_DIR:?Submit from film_rupture}
export OMP_NUM_THREADS=1
export TZ=Asia/Shanghai
test -f "$BASE/BUILD_OK"
test -f "$BASE/SMOKE_OK"
TASK=${SLURM_ARRAY_TASK_ID:?}
[[ "$TASK" =~ ^[0-9]+$ ]] && (( TASK >= 0 && TASK <= 20 ))
IFS=$'\t' read -r ID VERSION H0 RUN SEED SLIP SOURCE SOURCE_SHA < <(sed -n "$((TASK + 2))p" "$BASE/tasks.tsv")
test "$TASK" = "$ID"
NAME=$(basename "$SOURCE" .c)
CASE="$BASE/$VERSION/h0_$H0/run_0$RUN"
mkdir -p "$CASE"
cd "$CASE"
# Refuse reuse of a started trajectory, including accidental resubmission.
mkdir .started
cp "$BASE/$SOURCE" source.c
printf '%s  source.c\n' "$SOURCE_SHA" | sha256sum -c -
sha256sum source.c "$BASE/bin/$NAME" > source_binary.sha256
cp "$BASE/tasks.tsv" input_tasks.tsv
printf 'task=%s\nversion=%s\nh0_ratio=%s\nrun=%s\nseed=%s\nslip_code=%s\ninitial_h=1\nlevel=6\ndt=0.003\nstep_limit=200000000\nstop_hmin=0.01\nh_floor=0.005\nhstar=0.01\nmpi_ranks=8\nvolume_reference=L0_squared\njob_id=%s\narray_id=%s\nnode=%s\nstarted=%s\n' \
 "$TASK" "$VERSION" "$H0" "$RUN" "$SEED" "$SLIP" "$SLURM_JOB_ID" "$SLURM_ARRAY_JOB_ID" "$SLURMD_NODENAME" "$(date -Is)" > meta.txt
set +e
/usr/bin/time -p -o walltime.txt srun --input=none --mpi=pmix_v3 --cpu-bind=cores -n 8 "$BASE/bin/$NAME" \
    "$SEED" 0.003 200000000 "$SLIP" 1.0 > run.log 2> stderr.log
RC=$?
set -e
printf 'exit_code=%s\nfinished=%s\n' "$RC" "$(date -Is)" > CONTROL_DONE
exit "$RC"
