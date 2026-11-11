#!/usr/bin/env bash
# The paper's consolidated benchmark: ONE laila job times every method the manuscript reports, on cs10 at k = 31, so no
# timing in the paper is compared across jobs (between-job variance was 3-6%). Hardware counters per command.
#   manuscript/run_paper_bench.sh            then: python3 manuscript/paper_bench_table.py <results dir>
# Order: a __load_only__ baseline per binary, then every method, twice, interleaved, with a load baseline between
# repeats and at the end. The spaced-seed bench runs once (3 interleaved rounds inside; it prints its own rates).
set -euo pipefail
B=/home/rcb112/project/claude/bench
H='cd kp/hash'
E_METHODS=(
  "ntHash batched (ntbatch)" "ntHash library roll"
  "e scalar roller (reference)" "e sym_seq L1" "e sym_seq L2"
  "gauss61 sym_seq L1" "gauss61 sym_seq L2" "m31 simd (avx512 else avx2)"
  "e branchy 1/1000" "e branchy 1/100" "e branchy 1/16" "e branchy 1/4"
  "e batched 1/1000" "e batched 1/100" "e batched 1/16" "e batched 1/4"
  "m31 branchy 1/1000" "m31 branchy 1/100" "m31 batched 1/16" "m31 batched 1/4"
)
M_METHODS=("m31 avx2 8-lane" "m31 avx512 16-lane")
e() { runs+=(--run "$H && ./bench_e \"\$DATA\" 31 1 \"only=$1\""); }
m() { runs+=(--run "$H && ./bench \"\$DATA\" 31 1 throughput \"only=$1\""); }
runs=()
for rep in 1 2; do
  e __load_only__
  for x in "${E_METHODS[@]}"; do e "$x"; done
  m __load_only__
  for x in "${M_METHODS[@]}"; do m "$x"; done
done
e __load_only__
m __load_only__
runs+=(--run "cd kp/ideas/eisenstein && ./bench_spaced_e \"\$DATA\" 60000000 3")
cd "$B/repos/cannabis-genome"
python3 -u scripts/laila_perf.py run --label kp_paper_bench --yes \
  --bundle ../kmer-polynomials:kp \
  --git "$B/data/kmer_poly_work/ntHash@c26bd4572a19de81e30d55042dbd33c1fd21d4b6:ntHash" \
  --build 'cd kp/hash && ./build_nthash.sh "$PWD/../.." && make bench bench_e NTHASH="$PWD/../../ntHash" && cd ../ideas/eisenstein && make bench_spaced_e NTHASH="$PWD/../../../ntHash"' \
  "${runs[@]}" \
  --data "$B/data/kmer_poly_work/cs10.fa"
