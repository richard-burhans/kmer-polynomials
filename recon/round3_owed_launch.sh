#!/usr/bin/env bash
# Launch round 3's owed runs (recon/round3_owed_job.sh) on laila through cannabis-genome's scripts/galaxy_run.py, as
# two jobs that may run side by side (each within its own 16-slot grant), or dry-run every phase locally on tiny
# stand-in genomes.
#
#   bash recon/round3_owed_launch.sh kp|iblt                 one job on laila; results in bench/data/galaxy_run/<stamp>_<label>
#   bash recon/round3_owed_launch.sh dry kp|iblt <fa_A> <fa_giga> <fa_draft> <dir>
#
# The inputs tarball is the 2026-11-11 run's (recon/kp10_vs_iblt_launch.sh: 2bit genomes, twoBitToFa, KMC 3.2.4).
# The bundle holds hash/ and recon/ (git's view) and km-peeler at the pinned commit (git archive of the local clone).
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BENCH=/home/rcb112/project/claude/bench/data
KMP_REPO=$BENCH/kmer_poly_work/km-peeler
KMP_PIN=b8ec39b789366f5f5f6d22121dc409c65900a408
INPUTS=$BENCH/kmer_poly_work/kp10_vs_iblt/inputs.tar.gz
GALAXY_RUN=/home/rcb112/project/claude/bench/repos/cannabis-genome/scripts/galaxy_run.py
declare -A PHASES=([kp]="setup kp_seeds kp_sweep k63" [iblt]="setup iblt_validate kmp_trials iblt_proto ext tally")

case ${1:-} in
  kp|iblt)
    job=$1
    [ -s "$INPUTS" ] || { echo "missing $INPUTS (bash recon/kp10_vs_iblt_launch.sh builds it)" >&2; exit 1; }
    runs=()
    for p in ${PHASES[$job]}; do runs+=(--run "bash recon/round3_owed_job.sh $p"); done
    python3 "$GALAXY_RUN" run --host laila --label "round3_owed_$job" --bundle "$ROOT/hash:hash" --bundle "$ROOT/recon:recon" \
      --git "$KMP_REPO@$KMP_PIN:kmpeeler" --build "bash recon/round3_owed_job.sh build" "${runs[@]}" \
      --data "$INPUTS:tar.gz" --collect out --max-hours 10 --yes
    ;;
  dry)
    job=$2 FA=$3 FG=$4 FD=$5 D=$6
    mkdir -p "$D/tb" "$D/st/2bit" "$D/st/kent" "$D/st/kmc/bin" "$D/st/kmc/lib"
    "$BENCH/kent-bin/faToTwoBit" "$FA" "$D/st/2bit/Cannbio-2.2bit"
    "$BENCH/kent-bin/faToTwoBit" "$FG" "$D/st/2bit/Cannbio2_GigaDB.2bit"
    "$BENCH/kent-bin/faToTwoBit" "$FD" "$D/st/2bit/Cannbio-2_draft.2bit"
    tar -xzf "$INPUTS" -C "$D/st" kent kmc
    tar -czf "$D/inputs.tar.gz" -C "$D/st" 2bit kent kmc
    rm -rf "$D/w" && mkdir -p "$D/w/kmpeeler"
    cp -r "$ROOT/hash" "$ROOT/recon" "$D/w/"
    git -C "$KMP_REPO" archive "$KMP_PIN" | tar -x -C "$D/w/kmpeeler"
    echo "$KMP_PIN" > "$D/w/kmpeeler/.git_commit"
    cd "$D/w"
    export DATA_1=$D/inputs.tar.gz GALAXY_SLOTS=2 NSEEDS=2 NSEEDS_POW2=1 PLAN_DIV=50 REPS_FIXED=4 REPS_PROTO=2 EXT_TRIALS=2 K63_MEM_GB=2
    bash recon/round3_owed_job.sh build > "$D/run.log" 2>&1
    for p in ${PHASES[$job]}; do echo "=== $p" >> "$D/run.log"; bash recon/round3_owed_job.sh $p >> "$D/run.log" 2>&1 || echo "phase $p exit $?" >> "$D/run.log"; done
    ;;
  *) sed -n 2,10p "$0"; exit 2 ;;
esac
