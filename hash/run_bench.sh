#!/usr/bin/env bash
# Run the KP-6/KP-7 vs ntHash benchmark on cs10, end to end.
#   hash/run_bench.sh [k=31] [rounds=5] [throughput]   (throughput: skip the stats; output file gets a .throughput suffix)
# Builds ntHash (pinned, build_nthash.sh), builds and tests kp7, extracts cs10 if needed, and writes the raw
# output to results/kp7_rolling_vs_nthash_<date>.raw.txt. Single-threaded.
set -euo pipefail
cd "$(dirname "$0")"
K=${1:-31}; ROUNDS=${2:-5}; MODE=${3:-}
B=/home/rcb112/project/claude/bench/data
W=$B/kmer_poly_work
FA=$W/cs10.fa
./build_nthash.sh "$W"
make -s test_kp7 bench
./test_kp7 | tail -1
[ -s "$FA" ] || "$B/kent-bin/twoBitToFa" "$B/hub-build/2bit/cs10_NCBI_RefSeq_softmasked.2bit" "$FA"
OUT=../results/kp7_rolling_vs_nthash_$(date +%F)${MODE:+.$MODE}.raw.txt
{
  echo "# $(date -Is) host $(hostname) cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
  echo "# gcc: $(gcc --version | head -1); flags: $(make -s -p -n 2>/dev/null | awk -F' = ' '/^OPT = /{print $2; exit}')"
  echo "# kp7 commit: $(git rev-parse --short HEAD)$(git diff --quiet -- . || echo '+dirty'); input: $FA ($(stat -c %s "$FA") bytes)"
  ./bench --nthash-classes
  ./bench "$FA" "$K" "$ROUNDS" $MODE
} | tee "$OUT"
