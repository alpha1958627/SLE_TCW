#!/bin/sh
# Explicit spatial variant; default remains the original FV curvature.
set -eu
cd "$(dirname "$0")"
root=$(pwd)
: "${BASILISK:?Set BASILISK to the compatible Basilisk src directory}"
case "${CURVATURE:-fv}" in
  fv) source_dir="$root"; source=sle_smid_fv_basilisk.c; library=basilisk_midpoint ;;
  fvsp) source_dir="$root/experiments/fvsp"; source=sle_smid_fvsp_basilisk.c; library=basilisk_fvsp ;;
  *) echo 'CURVATURE must be fv or fvsp' >&2; exit 2 ;;
esac
mkdir -p compiled
test ! -e "compiled/$library.so"
case "${MPI:-0}" in
  0) set -- ;;
  1) : "${CC99:?Set CC99 to the MPI C compiler}"; set -- -D_MPI=1 ;;
  *) echo 'MPI must be 0 or 1' >&2; exit 2 ;;
esac
(cd "$source_dir"; "$BASILISK/qcc" -I"$BASILISK" -D_GNU_SOURCE "$@" -O3 -fno-fast-math \
  -ffp-contract=off -shared -fPIC "$source" -o "$root/compiled/$library.so" -lm)
sha256sum "$source_dir/$source" "compiled/$library.so" > "compiled/$library.sha256"
