#!/bin/bash
#SBATCH -J ls02NR21check
#SBATCH -p C64
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=8G
#SBATCH -t 00:10:00
#SBATCH --no-requeue
#SBATCH -o smoke-%j.out
#SBATCH -e smoke-%j.err
if [[ -f /etc/profile.d/dawning.sh ]]; then source /etc/profile.d/dawning.sh; fi
set -euo pipefail
BASE=${SLURM_SUBMIT_DIR:?Submit from film_rupture}
export OMP_NUM_THREADS=1
cd "$BASE"
mkdir smoke_v3
while IFS=$'\t' read -r ID VERSION H0 RUN SEED SLIP SOURCE SOURCE_SHA; do
    mkdir "$BASE/smoke_v3/task_$ID"
    cd "$BASE/smoke_v3/task_$ID"
    NAME=$(basename "$SOURCE" .c)
    srun --input=none --mpi=pmix_v3 --cpu-bind=cores -n 8 "$BASE/bin/$NAME" "$SEED" .003 1024 "$SLIP" 1 >run.log 2>stderr.log
done < <(tail -n +2 "$BASE/tasks.tsv")
cd "$BASE"
python3 check_smoke.py
