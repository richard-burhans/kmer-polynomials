#!/usr/bin/env bash
# Build km-peeler (Shibuya, Belazzougui & Kucherov, WABI 2022; github.com/yhhshb/km-peeler, MIT licence) at a pinned
# commit, for a given payload length.
#   recon/build_kmpeeler.sh <dir> <algorithm syncmers|minimizers> <k> <z>
# km-peeler's configure.py writes the payload header and then runs plain `make`, which fails on gcc >= 10: its tentative
# definitions collide under the -fno-common default. We run configure for the header, then make with -fcommon added.
# Their source is not modified.
set -euo pipefail
DIR=${1:?dir}; ALG=${2:?algorithm}; K=${3:?k}; Z=${4:?z}
PIN=b8ec39b789366f5f5f6d22121dc409c65900a408
[ -d "$DIR/.git" ] || git clone -q https://github.com/yhhshb/km-peeler.git "$DIR"
git -C "$DIR" checkout -q "$PIN"
cd "$DIR"
python3 configure.py sequences -a "$ALG" -k "$K" -z "$Z" > configure.log 2>&1 || true   # its make step is expected to fail
make clean > /dev/null 2>&1 || true
make ibltseq CFLAGS="-pedantic -Wall -Werror -Wno-format -O3 -DNDEBUG=1 -fcommon" > make.log 2>&1
echo "km-peeler $PIN built: $ALG k=$K z=$Z -> $DIR/ibltseq"
