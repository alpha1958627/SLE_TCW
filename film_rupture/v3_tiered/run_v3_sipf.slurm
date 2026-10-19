#!/bin/bash
# 2026-10-20 v3 批: ls0.2 (s1=0.2) × h0=0.5-1.2 (步长0.1) × 5 run, 换算口径, 分档 dt
# 分档 cap (dt-limit test 2026-10-19): hmin>=0.3 -> 1e-1; <0.3 -> 1e-2; <0.1 -> 3e-3
# 输出: hfield 每100t + 截面 hmin每0.01一档(~99帧) + 每1000t(~100帧) + 破碎快照
#SBATCH -J s02v3
#SBATCH -p C64
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -c 1
#SBATCH --array=0-39%40
#SBATCH --mem=8G
#SBATCH -t 14-00:00:00
#SBATCH -o /public/home/ljh/20261020_ls02_v3_h05_h12/slurm-%A_%a.out
#SBATCH -e /public/home/ljh/20261020_ls02_v3_h05_h12/slurm-%A_%a.err

set -euo pipefail

BASE=/public/home/ljh/20261020_ls02_v3_h05_h12
BASILISK=/public/home/ljh/pnas_sharp_interface_20260812/basilisk/src
export BASILISK OMP_NUM_THREADS=1
export CC99=/public/software/mpi/openmpi/gnu/4.0.3/bin/mpicc
export PATH="$BASILISK:/usr/bin:/bin:$PATH"

H0S=(0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2)
SLIPC=(0.4 0.33333333333333333 0.28571428571428571 0.25 0.22222222222222222 0.2 0.18181818181818182 0.16666666666666667)
DT=1e-1
STEPS=200000000

TASK=${SLURM_ARRAY_TASK_ID:?}
HIDX=$((TASK / 5))
RUN=$((TASK % 5 + 1))
H0=${H0S[$HIDX]}
SLIP=${SLIPC[$HIDX]}
SEED=$((2026112100 + HIDX * 10 + RUN))

CASE="$BASE/h0_${H0}/run_0${RUN}"
mkdir -p "$CASE"
cd "$CASE"
cp "$BASE/source_v3_h${H0}.c" source.c
sha256sum source.c > source.sha256

"$BASILISK/qcc" -D_GNU_SOURCE -D_MPI=1 -Wall -O2 source.c -o thinfilm \
    -L"$BASILISK/gl" -lglutils -lfb_tiny -lm

printf 'date=2026-10-20\nversion=v3\npurpose=ls0.2 h0=0.5-1.2 step 0.1 x 5 run; tiered dt cap (1e-1/1e-2/3e-3 by hmin 0.3/0.1); all fixes (yfix2+y-index+conservation)\nconvention=code_h_init=1.0 A=0.0013159/h0^2 varphi=0.01/h0^2 L0=100/h0 hstar=0.01\nslip_s1=0.2\nslip_code=%s\nh0=%s\nlevel=6\ndt tiered: hmin>=0.3->1e-1, <0.3->1e-2, <0.1->3e-3\nstep_limit=%s\nnrank=8\nseed=%s\nrun=%s\nlog_every=1000\nstop_hmin=0.01\nhfield_periodic=every 100 t\nsections=hmin every 0.01 (sect_hmin) + every 1000 t (sect_t)\nrupture_snapshot=per-rank\ncase=%s\n' \
  "$SLIP" "$H0" "$STEPS" "$SEED" "$RUN" "$CASE" > meta.txt

set +e
/usr/bin/time -p -o walltime.txt srun --mpi=pmix_v3 --cpu-bind=cores -n 8 ./thinfilm \
    "$SEED" "$DT" "$STEPS" "$SLIP" 1.0 > run.log 2> stderr.log
RC=$?
set -e
printf 'exit_code=%s\n' "$RC" > CONTROL_DONE
exit "$RC"
