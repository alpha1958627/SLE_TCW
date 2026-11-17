#!/bin/bash
set -euo pipefail
: "${BASILISK:?Set BASILISK to the Basilisk src directory}"
export PATH="$BASILISK:$PATH" CC99=gcc
mkdir -p build_tmp
export TMPDIR="$PWD/build_tmp"
qcc -disable-dimensions -D_GNU_SOURCE -Wall -O2 -g -fopenmp layer_v4.c -o layer_v4 -lm >build.out 2>build.err
qcc -disable-dimensions -D_GNU_SOURCE -Wall -O2 -g -fopenmp legacy_control.c -o legacy_control -lm >legacy_build.out 2>legacy_build.err
sha256sum layer_v4.c layer_v4 legacy_control.c legacy_control run_validation.py >inputs.sha256
