#!/bin/sh
set -eu
cd "$(dirname "$0")"
: "${BASILISK:?Set BASILISK to the Basilisk source directory}"
: "${CC99:=gcc}"
export BASILISK CC99
gcc -std=c99 -Wall -Wextra -O2 src/test_activation.c -lm -o test_activation
./test_activation
(cd src && "$BASILISK/qcc" -disable-dimensions -D_GNU_SOURCE -O2 droplet.c -o ../droplet -lm)
sha256sum src/droplet.c src/persistent_activation.h droplet > build.sha256
