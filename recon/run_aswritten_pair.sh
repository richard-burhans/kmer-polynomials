#!/usr/bin/env bash
# KP-10 on as-written k-mers for one genome pair (recon/aswritten.py): build the value sets, reconcile, verify.
#
#   bash recon/run_aswritten_pair.sh <genomeA> <genomeB> <config>...
#
# Each config is k:mode:bits:max_rounds, e.g. 31:sym:8:12 (strand-symmetric S, 2^8 buckets), 31:fwd:8:1 (the
# strand-specific control, same sketch, one round), 63:sym:auto:12 (k > 32, bucket count from the exact value D).
# Builds run in parallel per config, then the reconciliations run in parallel. Progress: one line per step in
# $LOGS/progress.txt (no process listing needed); each step's output in $LOGS/<step>.log; each run's JSON is copied to
# $RESULTS. Environment: PYTHON (default python3), SEED_X / SEED_LAM (default 0x5eed / 0x1a3d), LOGS, RESULTS.
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
A=$1 B=$2; shift 2
PY=${PYTHON:-python3}
SX=${SEED_X:-0x5eed} SL=${SEED_LAM:-0x1a3d}
LOGS=${LOGS:-/home/rcb112/project/claude/bench/data/kmer_poly_work/round1_owed/logs}
RESULTS=${RESULTS:-$ROOT/results/round1_owed_2026-11-11}
mkdir -p "$LOGS" "$RESULTS"
cd "$ROOT"
step() { # name, command...
  local name=$1; shift
  echo "$(date -u +%FT%TZ) start $name" >> "$LOGS/progress.txt"
  "$@" > "$LOGS/$name.log" 2>&1
  local rc=$?
  echo "$(date -u +%FT%TZ) done  $name rc=$rc" >> "$LOGS/progress.txt"
  return $rc
}
make -s -C hash kp7e_stream
for cfg in "$@"; do
  IFS=: read -r k mode bits rounds <<< "$cfg"
  for g in "$A" "$B"; do
    step "build_${g}_k${k}_${mode}" "$PY" -m recon.aswritten build "$g" --k "$k" --mode "$mode" --seed-x "$SX" \
      --seed-lam "$SL" &
  done
  wait
done
for cfg in "$@"; do
  IFS=: read -r k mode bits rounds <<< "$cfg"
  tag="${A}__${B}_k${k}_${mode}_b${bits}_r${rounds}"
  step "run_$tag" "$PY" -m recon.aswritten run "$A" "$B" --k "$k" --mode "$mode" --seed-x "$SX" --seed-lam "$SL" \
    --bits "$bits" --max-rounds "$rounds" --out "$RESULTS/run_$tag.json" &
done
wait
echo "$(date -u +%FT%TZ) all runs finished" >> "$LOGS/progress.txt"
