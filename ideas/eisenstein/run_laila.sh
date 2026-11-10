#!/usr/bin/env bash
# Laila timing for the Eisenstein revisit: KP-3 (de Bruijn navigation at k = 31, 47, 63; first 20 Mb) and KP-14 (spaced
# seeds: Eisenstein vs Gaussian vs ntHash2 SeedNtHash, 58.8 Mb, 3 interleaved rounds). One job, counters per command.
#   ideas/eisenstein/run_laila.sh
set -euo pipefail
B=/home/rcb112/project/claude/bench
cd "$B/repos/cannabis-genome"
E='cd kp/ideas/eisenstein'
python3 -u scripts/laila_perf.py run --label kp_eisenstein_ideas --yes \
  --bundle ../kmer-polynomials:kp \
  --git "$B/data/kmer_poly_work/ntHash@c26bd4572a19de81e30d55042dbd33c1fd21d4b6:ntHash" \
  --build 'cd kp/hash && ./build_nthash.sh "$PWD/../.." && cd ../ideas/eisenstein && make dbg_e31 dbg_e47 dbg_e63 bench_spaced_e NTHASH="$PWD/../../../ntHash"' \
  --run "$E && ./dbg_e31 \"\$DATA\" 20000000" \
  --run "$E && ./dbg_e47 \"\$DATA\" 20000000" \
  --run "$E && ./dbg_e63 \"\$DATA\" 20000000" \
  --run "$E && ./bench_spaced_e \"\$DATA\" 60000000 3" \
  --data "$B/data/kmer_poly_work/cs10_head60M.fa"
