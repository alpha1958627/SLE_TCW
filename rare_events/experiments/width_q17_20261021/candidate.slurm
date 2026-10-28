#!/bin/bash
#SBATCH -p C64
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=6G
#SBATCH -t 4-00:00:00
#SBATCH --no-requeue
set -eo pipefail
source /etc/profile.d/dawning.sh
source /public/software/profile.d/mpi_openmpi-gnu-4.0.3.sh
set -u
cd "${SLURM_SUBMIT_DIR:?}"
export OMP_NUM_THREADS=1
python3 admin/run_candidate.py --root "$PWD" --mode "${MODE:?}"
