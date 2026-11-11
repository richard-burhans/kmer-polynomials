#!/usr/bin/env bash
# Launch the KP-10 vs km-peeler comparison (recon/kp10_vs_iblt_job.sh) on laila through cannabis-genome's
# scripts/laila_perf.py, or dry-run every phase locally on tiny stand-in genomes.
#
#   bash recon/kp10_vs_iblt_launch.sh laila             full size, on laila (results in bench/data/laila_perf/<stamp>_kp10_vs_iblt)
#   bash recon/kp10_vs_iblt_launch.sh dry <fa_A> <fa_giga> <fa_draft> <dir>
#                                                        same phases, same bundle and inputs layout, local, small
#
# The inputs tarball (gzip -1) holds 2bit/<genome>.2bit for the three assemblies, kent/twoBitToFa, and KMC 3.2.4 from
# the round-1 bioconda env (bin/kmc, bin/kmc_tools and the lib/ it loads through $ORIGIN/../lib).
# The bundle holds hash/ and recon/ of this repo (git's view: tracked + untracked-not-ignored) and km-peeler at the
# pinned commit, exported by `git archive` from the local clone of github.com/yhhshb/km-peeler.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BENCH=/home/rcb112/project/claude/bench/data
KMP_REPO=$BENCH/kmer_poly_work/km-peeler
KMP_PIN=b8ec39b789366f5f5f6d22121dc409c65900a408
KMC_ENV=$BENCH/kmer_poly_work/round1_owed/env_kmc
LAILA_PERF=/home/rcb112/project/claude/bench/repos/cannabis-genome/scripts/laila_perf.py
GENOMES="Cannbio-2 Cannbio2_GigaDB Cannbio-2_draft"
PHASES=("setup" "kp10_timed giga" "iblt_timed giga" "kp10_timed draft" "iblt_timed draft" "kp10_sweep" "iblt_sweep"
        "iblt_protocol" "cleanup")

inputs() { # out.tar.gz twobit_dir
  local out=$1 tb=$2 st
  st=$(mktemp -d "$(dirname "$out")/inputs.XXXX")
  mkdir -p "$st/2bit" "$st/kent" "$st/kmc/bin" "$st/kmc/lib"
  for g in $GENOMES; do ln -s "$tb/$g.2bit" "$st/2bit/$g.2bit"; done
  ln -s "$BENCH/kent-bin/twoBitToFa" "$st/kent/twoBitToFa"
  ln -s "$KMC_ENV/bin/kmc" "$KMC_ENV/bin/kmc_tools" "$st/kmc/bin/"
  ln -s "$KMC_ENV/lib/libstdc++.so.6" "$KMC_ENV/lib/libgcc_s.so.1" "$st/kmc/lib/"
  tar -chf - -C "$st" 2bit kent kmc | gzip -1 > "$out"   # -h: the files, not the links
  rm -rf "$st"
}

case ${1:-} in
  laila)
    W=$BENCH/kmer_poly_work/kp10_vs_iblt
    mkdir -p "$W"
    [ -s "$W/inputs.tar.gz" ] || inputs "$W/inputs.tar.gz" "$BENCH/hub-build/2bit"
    runs=()
    for p in "${PHASES[@]}"; do runs+=(--run "bash recon/kp10_vs_iblt_job.sh $p"); done
    python3 "$LAILA_PERF" run --bundle "$ROOT/hash:hash" --bundle "$ROOT/recon:recon" --git "$KMP_REPO@$KMP_PIN:kmpeeler" \
      --build "bash recon/kp10_vs_iblt_job.sh build" "${runs[@]}" --data "$W/inputs.tar.gz" --data-ext tar.gz \
      --label kp10_vs_iblt --max-hours 12 --yes
    ;;
  dry)
    FA=$2 FG=$3 FD=$4 D=$5
    mkdir -p "$D/tb" "$D/w"
    "$BENCH/kent-bin/faToTwoBit" "$FA" "$D/tb/Cannbio-2.2bit"
    "$BENCH/kent-bin/faToTwoBit" "$FG" "$D/tb/Cannbio2_GigaDB.2bit"
    "$BENCH/kent-bin/faToTwoBit" "$FD" "$D/tb/Cannbio-2_draft.2bit"
    inputs "$D/inputs.tar.gz" "$D/tb"
    rm -rf "$D/w" && mkdir -p "$D/w"
    cp -r "$ROOT/hash" "$ROOT/recon" "$D/w/"
    mkdir -p "$D/w/kmpeeler"
    git -C "$KMP_REPO" archive "$KMP_PIN" | tar -x -C "$D/w/kmpeeler"
    echo "$KMP_PIN" > "$D/w/kmpeeler/.git_commit"
    cd "$D/w"
    export DATA=$D/inputs.tar.gz NSEEDS=${NSEEDS:-2} TRIALS=${TRIALS:-2} JOBS=${JOBS:-1} PROTO_REPS=${PROTO_REPS:-1} \
      KMC_THREADS=2 KMC_MEM_GB=2
    bash recon/kp10_vs_iblt_job.sh build
    for p in "${PHASES[@]}"; do echo "=== $p"; bash recon/kp10_vs_iblt_job.sh $p; done
    ;;
  *) sed -n 2,13p "$0"; exit 2 ;;
esac
