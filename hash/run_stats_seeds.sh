#!/usr/bin/env bash
# Multi-seed sampling statistics on cs10 (KP-7 next steps, item 3).
#   hash/run_stats_seeds.sh [k=31] [seeds=5] [label]
# Writes results/kp7_stats_seeds_<date>[.<label>].tsv; REFUSES to overwrite. Untimed, so a shared host is fine.
set -euo pipefail
cd "$(dirname "$0")"
K=${1:-31}; SEEDS=${2:-5}; LABEL=${3:-}
B=/home/rcb112/project/claude/bench/data
W=$B/kmer_poly_work
FA=$W/cs10.fa
OUT=../results/kp7_stats_seeds_$(date +%F)${LABEL:+.$LABEL}.tsv
[ -e "$OUT" ] && { echo "refusing to overwrite $OUT: pass a label (3rd argument)" >&2; exit 1; }
./build_nthash.sh "$W"
make -s test_kp7pack stats_seeds
./test_kp7pack | tail -1
[ -s "$FA" ] || "$B/kent-bin/twoBitToFa" "$B/hub-build/2bit/cs10_NCBI_RefSeq_softmasked.2bit" "$FA"
{
  echo "# $(date -Is) kp7 commit $(git rev-parse --short HEAD)$(git diff --quiet -- . || echo '+dirty'); input $FA"
  ./stats_seeds "$FA" "$K" "$SEEDS"
} | tee "$OUT"
