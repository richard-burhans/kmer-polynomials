#!/usr/bin/env bash
# Run the hybrid v2 (kp7pack) benchmark on cs10, end to end.
#   hash/run_bench_pack.sh [k=31] [rounds=5] [throughput|""] [label]
# Builds ntHash (pinned, build_nthash.sh), builds and tests kp7pack, extracts cs10 if needed, and writes the raw output
# to results/kp7_pack_<date>[.<mode>][.<label>].raw.txt. REFUSES to overwrite an existing file (pass a label).
set -euo pipefail
cd "$(dirname "$0")"
K=${1:-31}; ROUNDS=${2:-5}; MODE=${3:-}; LABEL=${4:-}
B=/home/rcb112/project/claude/bench/data
W=$B/kmer_poly_work
FA=$W/cs10.fa
OUT=../results/kp7_pack_$(date +%F)${MODE:+.$MODE}${LABEL:+.$LABEL}.raw.txt
[ -e "$OUT" ] && { echo "refusing to overwrite $OUT: pass a label (4th argument)" >&2; exit 1; }
./build_nthash.sh "$W"
make -s test_kp7pack bench_pack
./test_kp7pack | tail -1
[ -s "$FA" ] || "$B/kent-bin/twoBitToFa" "$B/hub-build/2bit/cs10_NCBI_RefSeq_softmasked.2bit" "$FA"
{
  echo "# $(date -Is) host $(hostname) cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)"
  echo "# gcc: $(gcc --version | head -1); flags: $(make -s -p -n 2>/dev/null | awk -F' = ' '/^OPT = /{print $2; exit}')"
  echo "# kp7pack commit: $(git rev-parse --short HEAD)$(git diff --quiet -- . || echo '+dirty'); input: $FA ($(stat -c %s "$FA") bytes)"
  ./bench_pack "$FA" "$K" "$ROUNDS" $MODE
} | tee "$OUT"
