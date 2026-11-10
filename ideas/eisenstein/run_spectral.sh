#!/usr/bin/env bash
# Eisenstein revisit of KP-4 and KP-13 (IDEAS.md): the same scripts and inputs as results/ideas_spectral_2026-11-10.md,
# with --encoding eisenstein. Same --rng, so KP-13 samples the same windows as the quarter-encoding run.
#   ideas/eisenstein/run_spectral.sh [kp4|kp13|all]
set -euo pipefail
cd "$(dirname "$0")/../spectral"
W=/home/rcb112/project/claude/bench/data/kmer_poly_work
PY=${PY:-$W/.venv/bin/python}
OUT=../../results/ideas_eisenstein
mkdir -p "$OUT"
what=${1:-all}
if [[ $what == kp4 || $what == all ]]; then
  "$PY" kp4_dft_bound.py "$W/spectral/cs10.NC_044378.1.fa" --encoding eisenstein --out "$OUT/kp4_eisenstein_cs10_NC_044378.1.json" > /dev/null
  "$PY" kp4_dft_bound.py "$W/spectral/cs10.NC_044378.1.fa" --encoding quarter --out "$OUT/kp4_quarter_cs10_NC_044378.1.json" > /dev/null
fi
if [[ $what == kp13 || $what == all ]]; then
  "$PY" kp13_fft_diag.py "$W/spectral/cs10.NC_044378.1.fa" "$W/KOMPa.chr7.fa" \
    /home/rcb112/project/claude/bench/data/multiz-lowmem/cache/KOMPa.over.chain \
    --encoding eisenstein --lastz-windows 0 --genome-search 0 --out "$OUT/kp13_eisenstein_cs10_vs_KOMPa_chr7.json" \
    > "$OUT/kp13_eisenstein_run.log"
fi
