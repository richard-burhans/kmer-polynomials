#!/usr/bin/env bash
# KP-10 (recon/kp10c.c) against Shibuya, Belazzougui & Kucherov's km-peeler IBLT (github.com/yhhshb/km-peeler, pinned
# b8ec39b), whole genome, on laila. Launched by recon/kp10_vs_iblt_launch.sh through cannabis-genome's
# scripts/laila_perf.py: each phase below is one laila_perf run, so pcount reports its wall and CPU time separately.
# Runs from the bundle root (hash/, recon/, kmpeeler/); $DATA is the inputs tarball (2bit genomes, twoBitToFa, KMC).
#
#   bash recon/kp10_vs_iblt_job.sh build                      compile kp10c, kmtools and km-peeler
#   bash recon/kp10_vs_iblt_job.sh setup                      FASTA, canonical sets, exact differences, KMC truth
#   bash recon/kp10_vs_iblt_job.sh kp10_timed PAIR            one KP-10 run at the true D, fixed seeds, ONE thread
#   bash recon/kp10_vs_iblt_job.sh iblt_timed PAIR            one km-peeler trial at n = 1.05 D, every step timed
#   bash recon/kp10_vs_iblt_job.sh kp10_sweep                 D-estimate sweep + NSEEDS OS-random seeds, both pairs
#   bash recon/kp10_vs_iblt_job.sh iblt_sweep                 TRIALS seeds per IBLT ratio, both pairs, in parallel
#   bash recon/kp10_vs_iblt_job.sh iblt_protocol              resend-at-double-size from n = 1.05 f D
#   bash recon/kp10_vs_iblt_job.sh cleanup
#
# PAIR is giga (Cannbio-2 vs Cannbio2_GigaDB) or draft (Cannbio-2 vs Cannbio-2_draft). Every result is one line
# "@@ <kind> <json>" in run.log; recon/kp10_vs_iblt_summarize.py reads them. Decisions inside a protocol never read
# the truth: KP-10 splits on its own failed checks, the IBLT resends when `ibltseq list` exits non-zero; the truth
# (kmtools / KMC 3 canonical set difference) only scores the outcome.
set -uo pipefail
ROOT=$PWD
W=$ROOT/work
BIN=$ROOT/bin
K=31
SEED_X=0x5eed
SEED_LAM=0x1a3d
NSEEDS=${NSEEDS:-30}
TRIALS=${TRIALS:-30}
JOBS=${JOBS:-48}
GENOMES="Cannbio-2 Cannbio2_GigaDB Cannbio-2_draft"
RATIOS_giga="0.95 1.0 1.05 1.1 1.25"
RATIOS_draft="1.0 1.05 1.25 1.5 2.0"
DEST_FACTORS="0.125 0.25 0.5 1 2 4"
PROTO_FACTORS="0.125 0.25 0.5 1"
PROTO_REPS=${PROTO_REPS:-10}
export TIMEFORMAT='%R %U %S'

pair_a() { echo Cannbio-2; }
pair_b() { case $1 in giga) echo Cannbio2_GigaDB ;; draft) echo Cannbio-2_draft ;; *) echo "bad pair $1" >&2; exit 2 ;; esac; }
true_d() { cat "$W/truth/$1/D"; }
scale() { awk -v f="$1" -v d="$2" 'BEGIN { printf "%d", f * d + 0.5 }'; }
emit() { printf '@@ %s %s\n' "$1" "$2"; }
now() { date +%s.%N; }
since() { awk -v a="$1" -v b="$(now)" 'BEGIN { printf "%.3f", b - a }'; }

build() {
  mkdir -p "$BIN"
  gcc --version | head -1
  gcc -O3 -march=native -Wall -Wextra -fopenmp -Ihash recon/kp10c.c hash/kp7e.c -o "$BIN/kp10c" || return 1
  gcc -O3 -march=native -Wall -Wextra recon/kmtools.c -o "$BIN/kmtools" || return 1
  # km-peeler: the header configure.py writes for `sequences -a kmers|syncmers -k 31` (configure_for_fixed_length_
  # sequences), byte for byte, then their make with -fcommon (gcc >= 10; recon/build_kmpeeler.sh). Source unmodified.
  printf '#ifndef COMPILE_OPTIONS_H\n#define COMPILE_OPTIONS_H\n\n#define STORE_SEQUENCES %d\n\n/*#define DEBUG*/\n\n#endif\n' $K \
    > kmpeeler/compile_options.h
  make -C kmpeeler ibltseq CFLAGS="-pedantic -Wall -Werror -Wno-format -O3 -DNDEBUG=1 -fcommon" > kmpeeler/make.log 2>&1 \
    || { tail -20 kmpeeler/make.log; return 1; }
  cp kmpeeler/ibltseq "$BIN/ibltseq"
  echo "km-peeler commit $(cat kmpeeler/.git_commit)"
  sha256sum recon/kp10c.c recon/kmtools.c hash/kp7e.c hash/kp7e.h recon/kp10_vs_iblt_job.sh
}

setup() {
  df -h . | tail -1
  echo "nproc $(nproc); $(free -g | sed -n 2p)"
  mkdir -p in "$W/fa" "$W/sets" "$W/truth" "$W/kmc/fa"
  tar -xf "$DATA" -C in || return 1
  md5sum in/2bit/*.2bit
  local t0
  for g in $GENOMES; do in/kent/twoBitToFa "in/2bit/$g.2bit" "$W/fa/$g.fa" & done
  wait
  # canonical k-mer sets: the IBLT's input and the first truth
  for g in $GENOMES; do
    ( { time "$BIN/kmtools" kset $K < "$W/fa/$g.fa" > "$W/sets/$g.u64" 2> "$W/sets/$g.kset.json"; } 2> "$W/sets/$g.kset.time" ) &
  done
  wait
  for g in $GENOMES; do
    read -r re us sy < "$W/sets/$g.kset.time"
    emit kset "{\"genome\": \"$g\", \"wall\": $re, \"user\": $us, \"sys\": $sy, \"out\": $(cat "$W/sets/$g.kset.json")}"
  done
  for p in giga draft; do
    local a b
    a=$(pair_a); b=$(pair_b $p)
    mkdir -p "$W/truth/$p"
    "$BIN/kmtools" diff "$W/sets/$a.u64" "$W/sets/$b.u64" "$W/truth/$p/amb.u64" "$W/truth/$p/bma.u64" 2> "$W/truth/$p/diff.json"
    emit diff "{\"pair\": \"$p\", \"out\": $(cat "$W/truth/$p/diff.json")}"
    echo $(( $(stat -c %s "$W/truth/$p/amb.u64") / 8 + $(stat -c %s "$W/truth/$p/bma.u64") / 8 )) > "$W/truth/$p/D"
  done
  # KMC 3 (recon/kmc_diff.sh, unchanged logic): the independent truth. One workdir, so Cannbio-2 is counted once.
  for g in $GENOMES; do ln -sf "$W/fa/$g.fa" "$W/kmc/fa/$g.fa"; done
  for p in giga draft; do
    local a b d
    a=$(pair_a); b=$(pair_b $p); d=$W/kmc/k$K
    t0=$(now)
    TWOBIT=$ROOT/in/2bit KENT=$ROOT/in/kent KMC_BIN=$ROOT/in/kmc/bin THREADS=${KMC_THREADS:-16} MEM_GB=${KMC_MEM_GB:-32} \
      bash recon/kmc_diff.sh "$a" "$b" $K "$W/kmc" > "$W/truth/$p/kmc_summary.json" || { echo "kmc_diff failed for $p"; tail -5 "$d"/*.log; }
    "$BIN/kmtools" kmcdump $K < "$d/amb.txt" > "$W/truth/$p/kmc_amb.u64" 2> "$W/truth/$p/kmcdump_amb.json"
    "$BIN/kmtools" kmcdump $K < "$d/bma.txt" > "$W/truth/$p/kmc_bma.u64" 2> "$W/truth/$p/kmcdump_bma.json"
    emit kmc "{\"pair\": \"$p\", \"wall\": $(since "$t0"), \"summary\": $(tr '\n' ' ' < "$W/truth/$p/kmc_summary.json"), \
\"same_amb\": $("$BIN/kmtools" same "$W/truth/$p/amb.u64" "$W/truth/$p/kmc_amb.u64"), \
\"same_bma\": $("$BIN/kmtools" same "$W/truth/$p/bma.u64" "$W/truth/$p/kmc_bma.u64")}"
  done
  rm -rf "$W/kmc/k$K"
  df -h . | tail -1
}

kp10_timed() {
  local p=$1 a b D
  a=$(pair_a); b=$(pair_b "$p"); D=$(true_d "$p")
  "$BIN/kp10c" --a "$W/fa/$a.fa" --b "$W/fa/$b.fa" --k $K --m 64 --c 2 --dest "$D" --seed-x $SEED_X --seed-lam $SEED_LAM \
    --threads 1 --truth-amb "$W/truth/$p/amb.u64" --truth-bma "$W/truth/$p/bma.u64" --label "timed_$p" \
    | while read -r l; do emit kp10_timed "$l"; done
}

# one km-peeler trial; prints one JSON object. Args: pair n seed tag
iblt_once() {
  local p=$1 n=$2 seed=$3 tag=$4 a b d t0 tb ta tdf tl rc bytes score
  a=$(pair_a); b=$(pair_b "$p")
  d=$W/tr/$tag; mkdir -p "$d"
  t0=$(now); "$BIN/kmtools" text $K < "$W/sets/$a.u64" 2>/dev/null | "$BIN/ibltseq" build -n "$n" -r 3 -s "$seed" -l $K -o "$d/a.iblt"; ta=$(since "$t0")
  t0=$(now); "$BIN/kmtools" text $K < "$W/sets/$b.u64" 2>/dev/null | "$BIN/ibltseq" build -n "$n" -r 3 -s "$seed" -l $K -o "$d/b.iblt"; tb=$(since "$t0")
  t0=$(now); "$BIN/ibltseq" diff -i "$d/a.iblt" -j "$d/b.iblt" -o "$d/d.iblt"; tdf=$(since "$t0")
  t0=$(now); "$BIN/ibltseq" list -i "$d/d.iblt" > "$d/list.txt" 2> "$d/list.err"; rc=$?; tl=$(since "$t0")
  bytes=$(stat -c %s "$d/a.iblt")
  score=$("$BIN/kmtools" listscore $K "$W/truth/$p/amb.u64" "$W/truth/$p/bma.u64" < "$d/list.txt")
  rm -rf "$d"
  echo "{\"pair\": \"$p\", \"n\": $n, \"seed\": $seed, \"iblt_bytes\": $bytes, \"list_exit\": $rc, \"wall_build_A\": $ta, \
\"wall_build_B\": $tb, \"wall_diff\": $tdf, \"wall_list\": $tl, \"score\": $score}"
}

iblt_timed() {
  local p=$1 a b D n d
  a=$(pair_a); b=$(pair_b "$p"); D=$(true_d "$p"); n=$(scale 1.05 "$D")
  d=$W/tr/timed_$p; mkdir -p "$d"
  { time "$BIN/kmtools" text $K < "$W/sets/$a.u64" > /dev/null 2>&1; } 2> "$d/t_text"
  { time bash -c "'$BIN/kmtools' text $K < '$W/sets/$a.u64' 2>/dev/null | '$BIN/ibltseq' build -n $n -r 3 -s 1 -l $K -o '$d/a.iblt'"; } 2> "$d/t_a"
  { time bash -c "'$BIN/kmtools' text $K < '$W/sets/$b.u64' 2>/dev/null | '$BIN/ibltseq' build -n $n -r 3 -s 1 -l $K -o '$d/b.iblt'"; } 2> "$d/t_b"
  { time "$BIN/ibltseq" diff -i "$d/a.iblt" -j "$d/b.iblt" -o "$d/d.iblt"; } 2> "$d/t_diff"
  { time "$BIN/ibltseq" list -i "$d/d.iblt" > "$d/list.txt" 2> "$d/list.err"; } 2> "$d/t_list"
  local rc=$?
  local js="{\"pair\": \"$p\", \"D\": $D, \"n\": $n, \"seed\": 1, \"list_exit\": $rc, \"iblt_bytes\": $(stat -c %s "$d/a.iblt")"
  for s in text a b diff list; do
    read -r re us sy < "$d/t_$s"
    js="$js, \"$s\": {\"real\": $re, \"user\": $us, \"sys\": $sy}"
  done
  js="$js, \"score\": $("$BIN/kmtools" listscore $K "$W/truth/$p/amb.u64" "$W/truth/$p/bma.u64" < "$d/list.txt")}"
  emit iblt_timed "$js"
  rm -rf "$d"
}

kp10_one() { # pair label threads dests seed_x seed_lam
  local p=$1 a b
  a=$(pair_a); b=$(pair_b "$p")
  "$BIN/kp10c" --a "$W/fa/$a.fa" --b "$W/fa/$b.fa" --k $K --m 64 --c 2 --dest "$4" --seed-x "$5" --seed-lam "$6" \
    --threads "$3" --truth-amb "$W/truth/$p/amb.u64" --truth-bma "$W/truth/$p/bma.u64" --label "$2" 2> /dev/null
}

kp10_sweep() {
  for p in giga draft; do
    (
      D=$(true_d "$p")
      dests=$(for f in $DEST_FACTORS; do scale "$f" "$D"; printf ','; done | sed 's/,$//')
      kp10_one "$p" "dest_sweep_$p" 8 "$dests" $SEED_X $SEED_LAM | while read -r l; do emit kp10_dest "$l"; done
      seq 1 "$NSEEDS" | xargs -P 4 -I{} bash "$0" kp10_one "$p" "os_seed_{}_$p" 4 "$D" os os \
        | while read -r l; do emit kp10_seed "$l"; done
    ) &
  done
  wait
}

iblt_sweep() {
  mkdir -p "$W/tr"
  for p in giga draft; do
    D=$(true_d "$p"); eval "rs=\$RATIOS_$p"
    for r in $rs; do
      n=$(scale "$r" "$D")
      for t in $(seq 1 "$TRIALS"); do echo "$p $n $((1000 + t)) $r"; done
    done
  done | xargs -P "$JOBS" -n 4 bash "$0" iblt_trial | while read -r l; do emit iblt_trial "$l"; done
}
iblt_trial() { # pair n seed ratio
  local js
  js=$(iblt_once "$1" "$2" "$3" "sw_$1_$2_$3")
  echo "{\"ratio\": $4, ${js#\{}"
}

# resend at double size until `ibltseq list` exits 0 (at most 8 rounds); the truth only scores the final listing
iblt_proto() { # pair f rep
  local p=$1 f=$2 rep=$3 D n round=0 total=0 js rc ns=""
  D=$(true_d "$p"); n=$(scale "$(awk -v f="$f" 'BEGIN { print 1.05 * f }')" "$D")
  while :; do
    round=$((round + 1))
    js=$(iblt_once "$p" "$n" "$(awk -v f="$f" -v r="$rep" -v k="$round" 'BEGIN { printf "%d", 100000 * k + 1000 * r + f * 1000 }')" "pr_${p}_${f}_${rep}_$round")
    total=$((total + $(echo "$js" | sed 's/.*"iblt_bytes": \([0-9]*\).*/\1/')))
    rc=$(echo "$js" | sed 's/.*"list_exit": \([0-9]*\).*/\1/')
    ns="$ns${ns:+, }$n"
    if [ "$rc" = 0 ] || [ $round -ge 8 ]; then break; fi
    n=$((2 * n))
  done
  echo "{\"pair\": \"$p\", \"f\": $f, \"rep\": $rep, \"D\": $D, \"rounds\": $round, \"n_per_round\": [$ns], \"total_bytes\": $total, \"last\": $js}"
}
iblt_protocol() {
  mkdir -p "$W/tr"
  for p in giga draft; do for f in $PROTO_FACTORS; do for r in $(seq 1 "$PROTO_REPS"); do echo "$p $f $r"; done; done; done \
    | xargs -P "$JOBS" -n 3 bash "$0" iblt_proto | while read -r l; do emit iblt_protocol "$l"; done
}

cleanup() { df -h . | tail -1; rm -rf "$W" in; df -h . | tail -1; }

cmd=$1; shift
case $cmd in
  build|setup|kp10_timed|iblt_timed|kp10_sweep|iblt_sweep|iblt_protocol|cleanup|kp10_one|iblt_trial|iblt_proto) "$cmd" "$@" ;;
  *) echo "unknown phase $cmd" >&2; exit 2 ;;
esac
