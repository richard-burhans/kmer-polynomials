#!/usr/bin/env bash
# The paper's frozen benchmark: ONE laila job times every number the manuscript quotes, at ONE commit, on cs10 at
# k = 31, so no timing in the paper is compared across jobs (between-job variance was 3-6%).
#   manuscript/run_paper_bench.sh [commit=HEAD] [repeats=5]
#   then: python3 manuscript/paper_bench_table.py <results dir>     (results land in bench/data/laila_perf/<stamp>_<label>)
# The bundle is `git archive <commit>` (laila_perf --git), so the timed code is exactly <commit>. Every binary is built
# in the job with the same OPT = -O3 -march=native (hash/Makefile, ideas/eisenstein/Makefile), ntHash at c26bd45.
# Each repeat runs a __load_only__ baseline per binary and every method once, one counted process each, in order; odd
# repeats run the whole list reversed, so drift does not favour one end. One more load run per binary closes the job.
# The table takes medians over repeats. bench_spaced_e runs once with 5 interleaved rounds and prints its own medians.
# Runs are sequential and single-threaded: one of the job's 16 slots is used.
#   bench_e     every-k-mer hashes and the KP-7 samplers as the paper first timed them (job 08653d9e6a7453a2)
#   bench       the CM31 SIMD rollers
#   bench_kept  the samplers in the kept-k-mer harness, the `ec 13` variant, and the round-1 baselines: a universal
#               hash on every k-mer, and the same packed-code gate followed by ntHash or the universal hash on kept k-mers
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
COMMIT=$(git -C "$REPO" rev-parse "${1:-HEAD}")
REPEATS=${2:-5}
B=/home/rcb112/project/claude/bench
H='cd kp/hash'
E_METHODS=(
  "ntHash batched (ntbatch)" "ntHash library roll"
  "e sym_seq L1" "e sym_seq L2" "gauss61 sym_seq L2"
  "e branchy 1/1000" "e branchy 1/100" "e batched 1/16" "e batched 1/4"
)
M_METHODS=("m31 avx2 8-lane" "m31 avx512 16-lane")
K_METHODS=(
  "ntHash batched" "uh all"
  "lib-e 8 br 1/1000" "lib-e 8 br 1/100" "lib-e 8 bat 1/1000" "lib-e 8 bat 1/100" "lib-e 8 bat 1/16" "lib-e 8 bat 1/4"
  "ec 13 br 1/1000" "ec 13 br 1/100"
  "nt 8 br 1/1000" "nt 8 br 1/100" "nt 8 bat 1/16" "nt 8 bat 1/4"
  "uh 0 br 1/1000" "uh 0 br 1/100" "uh 0 bat 1/16" "uh 0 bat 1/4"
)
L=()
L+=("e|__load_only__"); for x in "${E_METHODS[@]}"; do L+=("e|$x"); done
L+=("m|__load_only__"); for x in "${M_METHODS[@]}"; do L+=("m|$x"); done
L+=("k|__load_only__"); for x in "${K_METHODS[@]}"; do L+=("k|$x"); done
cmd() {
  case "${1%%|*}" in
    e) echo "$H && ./bench_e \"\$DATA\" 31 1 \"only=${1#*|}\"" ;;
    m) echo "$H && ./bench \"\$DATA\" 31 1 throughput \"only=${1#*|}\"" ;;
    k) echo "$H && ./bench_kept \"\$DATA\" 31 \"only=${1#*|}\"" ;;
  esac
}
runs=()
for ((rep = 0; rep < REPEATS; rep++)); do
  if ((rep % 2 == 0)); then for x in "${L[@]}"; do runs+=(--run "$(cmd "$x")"); done
  else for ((j = ${#L[@]} - 1; j >= 0; j--)); do runs+=(--run "$(cmd "${L[j]}")"); done; fi
done
for b in e m k; do runs+=(--run "$(cmd "$b|__load_only__")"); done
runs+=(--run "$H && ./bench_kept \"\$DATA\" 31 check")
runs+=(--run "cd kp/ideas/eisenstein && ./bench_spaced_e \"\$DATA\" 60000000 5")
[ -n "${DRY:-}" ] && { printf '%s\n' "${runs[@]}" | grep -v '^--run$'; exit 0; }   # DRY=1: list the runs, launch nothing
echo "paper bench: ${#L[@]} entries x $REPEATS repeats + 3 load runs + check + spaced; kmer-polynomials at $COMMIT" >&2

cd "$B/repos/cannabis-genome"
python3 -u scripts/laila_perf.py run --label kp_paper_bench --yes \
  --git "$REPO@$COMMIT:kp" \
  --git "$B/data/kmer_poly_work/ntHash@c26bd4572a19de81e30d55042dbd33c1fd21d4b6:ntHash" \
  --build 'cd kp/hash && ./build_nthash.sh "$PWD/../.." && make bench bench_e bench_kept test_kp7kept NTHASH="$PWD/../../ntHash" && ./test_kp7kept | tail -2 && cd ../ideas/eisenstein && make bench_spaced_e NTHASH="$PWD/../../../ntHash"' \
  "${runs[@]}" \
  --data "$B/data/kmer_poly_work/cs10.fa"
