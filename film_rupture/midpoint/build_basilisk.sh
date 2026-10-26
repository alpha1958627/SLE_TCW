#!/bin/sh
# Run from any directory. Core source stays byte-for-byte frozen.
set -eu
cd "$(dirname "$0")"
: "${BASILISK:?Set BASILISK to the compatible Basilisk src directory}"
mkdir -p compiled
test ! -e compiled/basilisk_midpoint.so
case "${MPI:-0}" in
  0) set -- ;;
  1) : "${CC99:?Set CC99 to the MPI C compiler}"; set -- -D_MPI=1 ;;
  *) echo 'MPI must be 0 or 1' >&2; exit 2 ;;
esac
"$BASILISK/qcc" -I"$BASILISK" -D_GNU_SOURCE "$@" -O3 -fno-fast-math \
  -ffp-contract=off -shared -fPIC basilisk_midpoint.c \
  -o compiled/basilisk_midpoint.so -lm
sha256sum basilisk_midpoint.c compiled/basilisk_midpoint.so > compiled/basilisk.sha256
