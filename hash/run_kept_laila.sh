#!/usr/bin/env bash
# Kept-k-mer evaluation on laila: ONE laila_perf job times every method of a plan, interleaved, with hardware counters.
#   hash/run_kept_laila.sh <plan> [commit=HEAD] [repeats=3]
#   then: python3 hash/kept_table.py <results dir>          (results land in bench/data/laila_perf/<stamp>_<label>)
# Plans:
#   sweep    every evaluator (batched) at 1/1000, 1/100, 1/16, 1/4, the current samplers, the bases (rate 0), batched
#            ntHash and the round-1 baselines (UH on every k-mer; gate + ntHash / gate + UH on kept k-mers)
#   confirm  the best evaluators against batched ntHash at the sketch rates, branchy and batched, and batch sizes
#            (5 repeats recommended)
# The bundle is `git archive <commit>` of this repository (laila_perf --git), so the timed code is exactly <commit>.
# Every binary is built in the job by hash/Makefile (OPT = -O3 -march=native for kp7, ntbatch and the baselines alike).
# Each repeat runs the plan's methods in order; odd repeats run it reversed, so drift does not favour one end.
set -euo pipefail
PLAN=${1:?plan: sweep or confirm}
REPO=$(cd "$(dirname "$0")/.." && pwd)
COMMIT=$(git -C "$REPO" rev-parse "${2:-HEAD}")
REPEATS=${3:-3}
B=/home/rcb112/project/claude/bench
H='cd kp/hash'

M=()
case "$PLAN" in
  sweep)
    M+=("ntHash batched" "uh all")
    for r in 0 1/1000 1/100 1/16 1/4; do M+=("lib-e 8 bat $r"); done
    for r in 0 1/1000 1/100; do M+=("lib-e 8 br $r"); done
    for r in 1/16 1/4; do M+=("lib-m 8 bat $r"); done
    M+=("ec 8 bat 0")
    for ev in "ec 8" "ec 11" "ec 13" "ec 16" "ep 11" "ep 16" "gp 11" "gp 16" "mp 11" "mp 16" "nt 8" "nt 16" "uh 0"; do
      for r in 1/1000 1/100 1/16 1/4; do M+=("$ev bat $r"); done
    done
    ;;
  confirm)
    # chosen from the sweep (job 3dddf48cdf81c0ca; results/kp7_kept_eval_2026-11-11.md): the current samplers, the
    # best chunk widths, batch sizes, CM31's tables, and the gate + ntHash / gate + UH baselines, branchy at sparse rates
    M+=("ntHash batched" "uh all")
    for r in 0 1/1000 1/100; do M+=("lib-e 8 br $r"); done
    for r in 0 1/1000 1/100 1/16 1/4; do M+=("lib-e 8 bat $r"); done
    for r in 1/1000 1/100; do M+=("ec 13 br $r"); done
    for r in 1/1000 1/100 1/16 1/4; do M+=("ec 13 bat $r"); done
    for r in 1/16 1/4; do M+=("ec 11 bat $r" "ec 8 bat $r b64" "ec 8 bat $r b4096" "lib-m 8 bat $r"); done
    for ev in "nt 8" "uh 0"; do
      for r in 1/1000 1/100; do M+=("$ev br $r"); done
      for r in 1/16 1/4; do M+=("$ev bat $r"); done
    done
    ;;
  *) echo "unknown plan $PLAN" >&2; exit 2 ;;
esac

runs=()
add() { runs+=(--run "$H && ./bench_kept \"\$DATA\" 31 \"only=$1\""); }
for ((rep = 0; rep < REPEATS; rep++)); do
  add __load_only__
  if ((rep % 2 == 0)); then for x in "${M[@]}"; do add "$x"; done
  else for ((j = ${#M[@]} - 1; j >= 0; j--)); do add "${M[j]}"; done; fi
done
add __load_only__
runs+=(--run "$H && ./bench_kept \"\$DATA\" 31 check")
[ -n "${DRY:-}" ] && { printf '%s\n' "${M[@]}"; exit 0; }   # DRY=1: list the plan, launch nothing
echo "plan $PLAN: ${#M[@]} methods x $REPEATS repeats + $((REPEATS + 1)) load runs + check; kmer-polynomials at $COMMIT" >&2

cd "$B/repos/cannabis-genome"
python3 -u scripts/laila_perf.py run --label "kp_kept_$PLAN" --yes \
  --git "$REPO@$COMMIT:kp" \
  --build 'cd kp/hash && make bench_kept test_kp7kept && ./test_kp7kept | tail -2 && make mutation_kept' \
  "${runs[@]}" \
  --data "$B/data/kmer_poly_work/cs10.fa"
