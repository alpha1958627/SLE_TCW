#!/bin/bash
set -euo pipefail
if [[ -f /etc/profile.d/dawning.sh ]]; then source /etc/profile.d/dawning.sh; fi
cd "$(dirname "$(readlink -f "$0")")"
: "${BASILISK:?Set BASILISK to the compatible Basilisk src directory}"
export CC99="${CC99:-mpicc}"
export PATH="$BASILISK:$PATH"
export OMP_NUM_THREADS=1
sha256sum -c sources.sha256
mkdir -p bin
for src in source_yfix_h1.05.c source_yfix_h1.10.c source_yfix_h1.15.c source_yfix_h1.20.c source_yold_h1.05.c; do
    name=$(basename "$src" .c)
    "$BASILISK/qcc" -D_GNU_SOURCE -D_MPI=1 -Wall -O2 "$src" -o "bin/$name" \
        -L"$BASILISK/gl" -lglutils -lfb_tiny -lm > "bin/$name.build.log" 2>&1
    ldd "bin/$name" > "bin/$name.ldd.txt"
    if grep -q 'not found' "bin/$name.ldd.txt"; then exit 3; fi
done
sha256sum bin/source_yfix_h* bin/source_yold_h1.05 > binaries.sha256
date -Is > BUILD_OK
