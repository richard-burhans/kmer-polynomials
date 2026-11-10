#!/usr/bin/env bash
# Build the ntHash baseline for bench.cpp, pinned, with the SAME flags as kp7 (Makefile OPT).
#   hash/build_nthash.sh [workdir]     # default: bench/data/kmer_poly_work
# Produces <workdir>/ntHash/build_o3/libnthash.a; the Makefile's NTHASH points at <workdir>/ntHash.
# A vendored <workdir>/ntHash with a .git_commit file (cannabis-genome scripts/laila_perf.py --git) is built in place.
set -euo pipefail
W=${1:-/home/rcb112/project/claude/bench/data/kmer_poly_work}
COMMIT=c26bd4572a19de81e30d55042dbd33c1fd21d4b6   # bcgsc/ntHash main, 2026-04-01 (after v2.3.0)
OPT="-O3 -march=native"                              # must match hash/Makefile OPT
mkdir -p "$W"
if [ -f "$W/ntHash/.git_commit" ]; then
  # a VENDORED tree (git archive, e.g. in a laila_perf bundle: no git, no network in the job); trust only the pin
  [ "$(cat "$W/ntHash/.git_commit")" = "$COMMIT" ] || { echo "vendored ntHash is not at $COMMIT" >&2; exit 1; }
  cd "$W/ntHash"
else
  [ -d "$W/ntHash/.git" ] || git clone -q https://github.com/bcgsc/ntHash.git "$W/ntHash"
  cd "$W/ntHash"
  git fetch -q origin "$COMMIT" 2>/dev/null || true
  git checkout -q "$COMMIT"
  [ "$(git rev-parse HEAD)" = "$COMMIT" ] || { echo "ntHash is not at $COMMIT" >&2; exit 1; }
fi
mkdir -p build_o3
g++ $OPT -std=c++17 -Iinclude -c src/kmer.cpp -o build_o3/kmer.o
g++ $OPT -std=c++17 -Iinclude -c src/seed.cpp -o build_o3/seed.o
rm -f build_o3/libnthash.a
ar rcs build_o3/libnthash.a build_o3/kmer.o build_o3/seed.o
echo "ntHash $COMMIT built with: g++ $OPT -std=c++17 -> $W/ntHash/build_o3/libnthash.a"
