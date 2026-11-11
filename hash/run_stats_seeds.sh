#!/usr/bin/env bash
# Multi-seed sampling statistics on cs10 (KP-7 next steps, item 3).
#   hash/run_stats_seeds.sh [k=31] [seeds=5] [label]
# Writes results/kp7_stats_seeds_<date>[.<label>].tsv; REFUSES to overwrite. Untimed, so a shared host is fine.
# Also logs every long gap (<out>.gaps.tsv) and merges them into distinct loci across runs (hash/merge_gaps.py:
# <out>.loci.tsv and <out>.loci.json), since one repeat gaps in many runs (round-1 adjudication L2 F1).
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
  ./stats_seeds "$FA" "$K" "$SEEDS" "gaps=${OUT%.tsv}.gaps.tsv"
} | tee "$OUT"
python3 merge_gaps.py "${OUT%.tsv}.gaps.tsv" --out "${OUT%.tsv}.loci.tsv" --json "${OUT%.tsv}.loci.json"
