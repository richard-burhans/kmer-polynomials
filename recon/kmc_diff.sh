#!/usr/bin/env bash
# Exact canonical k-mer set difference of two hub genomes by KMC 3 (any k up to 256): the independent truth for
# recon/aswritten.py at k > 32.
#
#   bash recon/kmc_diff.sh <genomeA> <genomeB> <k> <workdir>
#
# KMC is run with -ci1 (keep singletons: we want the SET of k-mers, not the solid ones) and -cs255 (counts are irrelevant, but kmc_tools cannot subtract counter-less sets made with -cs1),
# canonical (KMC's default), from FASTA made by twoBitToFa. KMC is case-blind and skips k-mers with N, as
# recon/kmerset.py and hash/kp7e_stream.c do. Outputs in <workdir>/k<k>/: <g>.kmc_{pre,suf}, amb.txt / bma.txt (dumped
# k-mers of A\B and B\A), and summary.json with the counts. KMC_BIN defaults to the bioconda env made for round 1;
# TWOBIT (the <g>.2bit directory) and KENT (twoBitToFa's directory) default to bench and are overridden on laila.
set -euo pipefail
A=$1 B=$2 K=$3 W=$4
KMC_BIN=${KMC_BIN:-/home/rcb112/project/claude/bench/data/kmer_poly_work/round1_owed/env_kmc/bin}
TWOBIT=${TWOBIT:-/home/rcb112/project/claude/bench/data/hub-build/2bit}
KENT=${KENT:-/home/rcb112/project/claude/bench/data/kent-bin}
THREADS=${THREADS:-16}
MEM_GB=${MEM_GB:-6}
D=$W/k$K
mkdir -p "$D/tmp" "$W/fa"
for g in "$A" "$B"; do
  [ -s "$W/fa/$g.fa" ] || "$KENT/twoBitToFa" "$TWOBIT/$g.2bit" "$W/fa/$g.fa"
  if [ ! -s "$D/$g.kmc_suf" ]; then
    "$KMC_BIN/kmc" -k"$K" -ci1 -cs255 -fm -m"$MEM_GB" -t"$THREADS" "$W/fa/$g.fa" "$D/$g" "$D/tmp" > "$D/$g.kmc.log" 2>&1
  fi
done
"$KMC_BIN/kmc_tools" -t"$THREADS" simple "$D/$A" "$D/$B" kmers_subtract "$D/amb" reverse_kmers_subtract "$D/bma" \
  > "$D/subtract.log" 2>&1
"$KMC_BIN/kmc_tools" transform "$D/amb" dump "$D/amb.txt" > /dev/null 2>&1
"$KMC_BIN/kmc_tools" transform "$D/bma" dump "$D/bma.txt" > /dev/null 2>&1
uniq_of() { grep -m1 "No. of unique counted k-mers" "$D/$1.kmc.log" | awk '{print $NF}'; }
NA=$(uniq_of "$A") NB=$(uniq_of "$B") AMB=$(wc -l < "$D/amb.txt") BMA=$(wc -l < "$D/bma.txt")
VER=$("$KMC_BIN/kmc" 2>&1 | grep -m1 -o "ver\. [0-9.]*" || true)
cat > "$D/summary.json" <<EOF
{"tool": "KMC $VER", "A": "$A", "B": "$B", "k": $K, "distinct_A": $NA, "distinct_B": $NB, "A_minus_B": $AMB,
 "B_minus_A": $BMA, "D": $((AMB + BMA)), "amb": "$D/amb.txt", "bma": "$D/bma.txt"}
EOF
cat "$D/summary.json"
