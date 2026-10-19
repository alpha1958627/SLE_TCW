#!/bin/bash
#SBATCH -J ruptureREcheck
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
cd "${SLURM_SUBMIT_DIR:?Submit from repository root}"
export OMP_NUM_THREADS=1
python3 -m rare_events.gates --workdir "validation/gate_${SLURM_JOB_ID}"
