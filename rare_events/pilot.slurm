#!/bin/bash
#SBATCH -J ruptureREpilot
#SBATCH -p C64
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --mem=24G
#SBATCH -t 3-00:00:00
#SBATCH --no-requeue
#SBATCH -o pilot-%j.out
#SBATCH -e pilot-%j.err
set -eo pipefail
source /etc/profile.d/dawning.sh
source /public/software/profile.d/mpi_openmpi-gnu-4.0.3.sh
set -u
cd "${SLURM_SUBMIT_DIR:?Submit from repository root}"
: "${VALIDATED_GATE:?Set VALIDATED_GATE to the successful numerical gate directory}"
export OMP_NUM_THREADS=1
python3 -m rare_events.campaign --workdir "runs/job_${SLURM_JOB_ID}" --gate "$VALIDATED_GATE"
