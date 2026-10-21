#!/bin/bash
#SBATCH -J q17-score-gate
#SBATCH -p C64
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH --no-requeue
#SBATCH -o gate-%j.out
#SBATCH -e gate-%j.err
set -eo pipefail
source /etc/profile.d/dawning.sh
source /public/software/profile.d/mpi_openmpi-gnu-4.0.3.sh
set -u
cd "${SLURM_SUBMIT_DIR:?}"
export OMP_NUM_THREADS=1
export BASILISK="${BASILISK:?Verified Basilisk source directory required}"
bash rare_events/build.sh
python3 validate_q17.py --workdir "validation/gate_${SLURM_JOB_ID}" --baseline "${BASELINE_ROOT:?}"
