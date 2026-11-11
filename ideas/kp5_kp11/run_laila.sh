#!/usr/bin/env bash
# Laila timing for KP-5 (inverted repeats: Eisenstein hash vs packed codes) and KP-11 (mismatch counts: (Z_4)^k NTT vs
# direct ball enumeration) on cs10 NC_044378.1. One job, counters per command; each program also prints its own
# phase timings (load excluded), which are the numbers to quote.
#   ideas/kp5_kp11/run_laila.sh
set -euo pipefail
B=/home/rcb112/project/claude/bench
F=$B/data/kmer_poly_work/cs10_NC_044378.1.fa
[ -s "$F" ] || awk '/^>/{p=($1==">NC_044378.1")} p' "$B/data/kmer_poly_work/cs10.fa" > "$F"
cd "$B/repos/cannabis-genome"
E='cd kp/ideas/kp5_kp11'
runs=()
for spec in "palin 20" "inv 20 100" "inv 40 100"; do for m in hash pack; do runs+=(--run "$E && ./ir \"\$DATA\" - $spec $m"); done; done
for kd in "10 3" "12 1" "12 2" "12 3"; do runs+=(--run "$E && ./mismatch \"\$DATA\" - $kd both"); done
runs+=(--run "$E && ./mismatch \"\$DATA\" - 13 3 fft")
python3 -u scripts/laila_perf.py run --label kp5_kp11 --yes \
  --bundle ../kmer-polynomials:kp \
  --build 'cd kp/ideas/kp5_kp11 && make ir mismatch && ./test.sh' \
  "${runs[@]}" \
  --data "$F"
