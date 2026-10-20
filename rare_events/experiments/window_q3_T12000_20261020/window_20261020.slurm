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
cd "${SLURM_SUBMIT_DIR:?Submit from frozen repository root}"
export OMP_NUM_THREADS=1
folder=admin/window_q3_T12000_20261020
python3 "$folder/run_window_20261020.py" --repo "$PWD" --plan "$folder/plan.json" --mode "${MODE:?}"
