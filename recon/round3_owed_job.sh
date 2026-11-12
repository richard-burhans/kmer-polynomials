#!/usr/bin/env bash
# Round 3's owed runs (reports/rounds/round-2/round.json "open"; results/round3_owed_2026-11-12.md), on laila through
# cannabis-genome's scripts/galaxy_run.py (recon/round3_owed_launch.sh). Runs from the bundle root (hash/, recon/,
# kmpeeler/); $DATA_1 is the inputs tarball of the 2026-11-11 run (2bit genomes, twoBitToFa, KMC 3.2.4).
#
#   job kp:   build setup kp_seeds kp_sweep k63          the sketch (runs 1, 3, 4, 6)
#   job iblt: build setup iblt_validate kmp_trials iblt_proto ext tally      the table (runs 2, 3, 5, 7, 8)
#
#   build          compile kp10c, kmtools, iblt_grow, km-peeler for 31-mers and for 51-base extended syncmers
#   setup          FASTA, canonical 31-mer sets and the exact differences of giga, draft and dvg
#   kp_seeds       run 1 + 4: kp10c, buckets sized exactly from D (--rule exact, no power of two), one-way protocol,
#                  30 OS-random seeds per pair; the power-of-two rule, 10 seeds; one fixed seed one-way vs two-way
#   kp_sweep       run 1: fixed seeds, first-round buckets from loads 0.5..1.0 per m in {48, 64, 80, 88, 96, 128},
#                  m at 2^16 (giga) / 2^8 (draft), and the D-estimate sweep 1/8..4 under the exact rule
#   k63            run 6: KMC 3 truth at k = 63, then kp10c k = 63 at the stated 0.8 m rule (2^9), 30 OS seeds
#   iblt_validate  whole-genome checks that make the cheap runs exact: km-peeler's d.iblt from the full sets equals
#                  its d.iblt from the difference sets (cmp), recon/iblt_grow.c writes km-peeler's a.iblt byte for
#                  byte, and its peeler lists what `ibltseq list` lists
#   kmp_trials     run 2 + 3: km-peeler itself with independent seeds, many trials per ratio, on giga, draft, the
#                  held-out pair dvg and held-out random subsets of giga's difference; max |counter| of every d.iblt
#   iblt_proto     run 5: recon/iblt_grow.c fixed / resend / grow / subtract on giga and draft
#   ext            run 7: km-peeler extended syncmers (their Section 3.5) on the whole genome, 10 independent seeds,
#                  with a trace of every missed k-mer and of every differing syncmer the listing did not return
#   tally          run 8: per-record tally of each difference
#
# Pairs: giga = Cannbio-2 vs Cannbio2_GigaDB, draft = Cannbio-2 vs Cannbio-2_draft, dvg = Cannbio-2_draft vs
# Cannbio2_GigaDB (held out: used by no earlier run). Every result is one line "@@ <kind> <json>" in run.log; files
# for the write-up go to out/ (collected). Parallelism stays within $GALAXY_SLOTS (16 on laila): at most 4 kp10c
# processes of 4 threads, 14 km-peeler trials, 16 iblt_grow threads.
set -uo pipefail
ROOT=$PWD
W=$ROOT/work
BIN=$ROOT/bin
OUT=$ROOT/out
K=31
SLOTS=${GALAXY_SLOTS:-16}
PKP=$(( SLOTS / 4 > 0 ? SLOTS / 4 : 1 ))     # kp10c processes of 4 threads; full-set km-peeler builds (validate)
PKMP=$(( SLOTS > 2 ? SLOTS - 2 : 1 ))       # km-peeler trials (build, build, diff, list: one process at a time each)
PEXT=$(( SLOTS / 3 > 0 ? SLOTS / 3 : 1 ))   # extended-syncmer trials
DATA=${DATA_1:-${DATA:-}}
SEED_X=0x5eed
SEED_LAM=0x1a3d
NSEEDS=${NSEEDS:-30}
NSEEDS_POW2=${NSEEDS_POW2:-10}
PLAN_DIV=${PLAN_DIV:-1}         # dry runs: divide every trial count (at least 1 each)
REPS_FIXED=${REPS_FIXED:-100}
REPS_PROTO=${REPS_PROTO:-30}
EXT_TRIALS=${EXT_TRIALS:-10}
GENOMES="Cannbio-2 Cannbio2_GigaDB Cannbio-2_draft"
export TIMEFORMAT='%R %U %S'
mkdir -p "$OUT"

pair_a() { case $1 in giga|draft) echo Cannbio-2 ;; dvg) echo Cannbio-2_draft ;; *) echo "bad pair $1" >&2; exit 2 ;; esac; }
pair_b() { case $1 in giga|dvg) echo Cannbio2_GigaDB ;; draft) echo Cannbio-2_draft ;; *) echo "bad pair $1" >&2; exit 2 ;; esac; }
true_d() { cat "$W/truth/$1/D"; }
emit() { printf '@@ %s %s\n' "$1" "$2"; }
ceil_div() { awk -v d="$1" -v q="$2" 'BEGIN { x = d / q; n = int(x); if (n < x) n++; if (n < 1) n = 1; printf "%d", n }'; }
scale() { awk -v f="$1" -v d="$2" 'BEGIN { printf "%d", f * d + 0.5 }'; }

build() {
  mkdir -p "$BIN"
  gcc --version | head -1
  gcc -O3 -march=native -Wall -Wextra -fopenmp -Ihash recon/kp10c.c hash/kp7e.c -o "$BIN/kp10c" -lm || return 1
  gcc -O3 -march=native -Wall -Wextra recon/kmtools.c -o "$BIN/kmtools" || return 1
  gcc -O3 -march=native -fopenmp -Ikmpeeler recon/iblt_grow.c kmpeeler/murmur3.c -o "$BIN/iblt_grow" -lm || return 1
  # km-peeler, source unmodified: the header its configure.py writes for a 31-base payload (k-mers) and for a 51-base
  # payload (extended syncmers, 2k - z with k = 31, z = 11), then its make with -fcommon (gcc >= 10)
  cp -r kmpeeler kmpeeler_ext
  local d p
  for d in kmpeeler:$K kmpeeler_ext:51; do
    p=${d#*:}; d=${d%%:*}
    printf '#ifndef COMPILE_OPTIONS_H\n#define COMPILE_OPTIONS_H\n\n#define STORE_SEQUENCES %d\n\n/*#define DEBUG*/\n\n#endif\n' "$p" \
      > "$d/compile_options.h"
    make -C "$d" ibltseq CFLAGS="-pedantic -Wall -Werror -Wno-format -O3 -DNDEBUG=1 -fcommon" > "$d/make.log" 2>&1 \
      || { tail -20 "$d/make.log"; return 1; }
  done
  cp kmpeeler/ibltseq "$BIN/ibltseq"
  cp kmpeeler_ext/ibltseq "$BIN/ibltseq_ext"
  echo "km-peeler commit $(cat kmpeeler/.git_commit 2>/dev/null)"
  sha256sum recon/kp10c.c recon/kmtools.c recon/iblt_grow.c recon/round3_plan.py recon/round3_owed_job.sh hash/kp7e.c
}

setup() {
  df -h . | tail -1
  mkdir -p in "$W/fa" "$W/sets" "$W/truth"
  tar -xf "$DATA" -C in || return 1
  md5sum in/2bit/*.2bit
  for g in $GENOMES; do in/kent/twoBitToFa "in/2bit/$g.2bit" "$W/fa/$g.fa" & done
  wait
  for g in $GENOMES; do "$BIN/kmtools" kset $K < "$W/fa/$g.fa" > "$W/sets/$g.u64" 2> "$W/sets/$g.kset.json" & done
  wait
  for g in $GENOMES; do emit kset "{\"genome\": \"$g\", \"out\": $(cat "$W/sets/$g.kset.json")}"; done
  local p a b
  for p in giga draft dvg; do
    a=$(pair_a $p); b=$(pair_b $p)
    mkdir -p "$W/truth/$p"
    "$BIN/kmtools" diff "$W/sets/$a.u64" "$W/sets/$b.u64" "$W/truth/$p/amb.u64" "$W/truth/$p/bma.u64" 2> "$W/truth/$p/diff.json"
    emit diff "{\"pair\": \"$p\", \"out\": $(cat "$W/truth/$p/diff.json")}"
    echo $(( $(stat -c %s "$W/truth/$p/amb.u64") / 8 + $(stat -c %s "$W/truth/$p/bma.u64") / 8 )) > "$W/truth/$p/D"
  done
  df -h . | tail -1
}

# ---------------- the sketch ----------------
kp() { # pair label threads dests seed_x seed_lam extra...
  local p=$1 l=$2 t=$3 d=$4 sx=$5 sl=$6 a b
  shift 6
  a=$(pair_a "$p"); b=$(pair_b "$p")
  "$BIN/kp10c" --a "$W/fa/$a.fa" --b "$W/fa/$b.fa" --k $K --dest "$d" --seed-x "$sx" --seed-lam "$sl" --threads "$t" \
    --truth-amb "$W/truth/$p/amb.u64" --truth-bma "$W/truth/$p/bma.u64" --label "$l" "$@" 2> /dev/null
}
kp_one() { kp "$@" | while read -r l; do emit kp "$l"; done; }

kp_seeds() {
  local p D
  {
    for p in giga draft; do
      D=$(true_d $p)
      for s in $(seq 1 "$NSEEDS"); do echo "$p exact_os_${s}_$p 4 $D os os --rule exact --oneway 1"; done
      for s in $(seq 1 "$NSEEDS_POW2"); do echo "$p pow2_os_${s}_$p 4 $D os os --rule pow2 --oneway 1"; done
      echo "$p fixed_exact_oneway_$p 4 $D $SEED_X $SEED_LAM --rule exact --oneway 1"
      echo "$p fixed_exact_twoway_$p 4 $D $SEED_X $SEED_LAM --rule exact --oneway 0"
    done
  } | xargs -P "$PKP" -L 1 bash "$0" kp_one
}

kp_sweep() {
  local p D m L dests pw
  {
    for p in giga draft; do
      D=$(true_d $p)
      [ $p = giga ] && pw=65536 || pw=256
      for m in 48 64 80 88 96 128; do
        if [ $m = 64 ]; then
          dests=""
          for L in 0.5 0.6 0.7 0.8 0.9 1.0; do dests="$dests,$D@$(ceil_div "$D" "$(awk -v l=$L -v m=$m 'BEGIN { print l * m }')")"; done
          dests="${dests#,},$D@$pw"
          echo "$p sweep_m${m}_$p 4 $dests $SEED_X $SEED_LAM --m $m --rule exact --oneway 1"
          dests=""
          for f in 0.125 0.25 0.5 2 4; do dests="$dests,$(scale $f "$D")"; done
          echo "$p dest_exact_m${m}_$p 4 ${dests#,} $SEED_X $SEED_LAM --m $m --rule exact --oneway 1"
        else
          dests=""
          for L in 0.6 0.8 1.0; do dests="$dests,$D@$(ceil_div "$D" "$(awk -v l=$L -v m=$m 'BEGIN { print l * m }')")"; done
          echo "$p sweep_m${m}_$p 4 ${dests#,},$D@$pw $SEED_X $SEED_LAM --m $m --rule exact --oneway 1"
        fi
      done
    done
  } | xargs -P "$PKP" -L 1 bash "$0" kp_one
}

k63() {
  local a=Cannbio-2 b=Cannbio-2_draft d=$W/kmc k=63
  mkdir -p "$d/fa"
  ln -sf "$W/fa/$a.fa" "$d/fa/$a.fa"; ln -sf "$W/fa/$b.fa" "$d/fa/$b.fa"
  TWOBIT=$ROOT/in/2bit KENT=$ROOT/in/kent KMC_BIN=$ROOT/in/kmc/bin THREADS=$SLOTS MEM_GB=${K63_MEM_GB:-48} \
    bash recon/kmc_diff.sh "$a" "$b" $k "$d" > "$W/kmc63_summary.json" || { echo "kmc_diff failed"; return 1; }
  emit kmc63 "$(tr '\n' ' ' < "$W/kmc63_summary.json")"
  local D
  D=$(( $(wc -l < "$d/k$k/amb.txt") + $(wc -l < "$d/k$k/bma.txt") ))
  echo "k = 63 truth D = $D"
  {
    for s in $(seq 1 "$NSEEDS"); do echo "k63_os_$s 4 $D os os"; done
    echo "k63_fixed 4 $D,$D@1024,$D@$(ceil_div "$D" 51.2) $SEED_X $SEED_LAM"
  } | xargs -P "$PKP" -L 1 bash "$0" k63_one "$d/k$k"
  rm -rf "$d/k$k"/*.kmc_* "$d/k$k/tmp"
}
k63_one() { # kmcdir label threads dests seed_x seed_lam
  "$BIN/kp10c" --a "$W/fa/Cannbio-2.fa" --b "$W/fa/Cannbio-2_draft.fa" --k 63 --dest "$4" --seed-x "$5" --seed-lam "$6" \
    --threads "$3" --truth-kmc-amb "$1/amb.txt" --truth-kmc-bma "$1/bma.txt" --label "$2" --rule pow2 --oneway 1 2> /dev/null \
    | while read -r l; do emit kp63 "$l"; done
}

# ---------------- the table ----------------
textof() { [ -s "$2" ] || "$BIN/kmtools" text $K < "$1" > "$2" 2> /dev/null; }

iblt_validate() {
  local p=giga a b D
  a=$(pair_a $p); b=$(pair_b $p); D=$(true_d $p)
  mkdir -p "$W/val"
  textof "$W/truth/$p/amb.u64" "$W/truth/$p/amb.txt"
  textof "$W/truth/$p/bma.u64" "$W/truth/$p/bma.txt"
  local i=0 r s
  python3 recon/round3_plan.py seeds 4 > "$W/val/seeds"
  for r in 1.0 1.0 1.05 1.05; do
    i=$((i + 1)); s=$(sed -n ${i}p "$W/val/seeds")
    echo "$i $(scale $r "$D") $s $r"
  done | xargs -P "$PKP" -L 1 bash "$0" iblt_validate_one
}
iblt_validate_one() { # i n seed ratio
  local i=$1 n=$2 s=$3 r=$4 p=giga a b v=$W/val/$1 t0
  a=$(pair_a $p); b=$(pair_b $p)
  mkdir -p "$v"
  t0=$(date +%s)
  "$BIN/kmtools" text $K < "$W/sets/$a.u64" 2> /dev/null | "$BIN/ibltseq" build -n "$n" -r 3 -s "$s" -l $K -o "$v/a.iblt"
  "$BIN/kmtools" text $K < "$W/sets/$b.u64" 2> /dev/null | "$BIN/ibltseq" build -n "$n" -r 3 -s "$s" -l $K -o "$v/b.iblt"
  "$BIN/ibltseq" diff -i "$v/a.iblt" -j "$v/b.iblt" -o "$v/dfull.iblt"
  "$BIN/ibltseq" build -i "$W/truth/$p/amb.txt" -n "$n" -r 3 -s "$s" -l $K -o "$v/da.iblt"
  "$BIN/ibltseq" build -i "$W/truth/$p/bma.txt" -n "$n" -r 3 -s "$s" -l $K -o "$v/db.iblt"
  "$BIN/ibltseq" diff -i "$v/da.iblt" -j "$v/db.iblt" -o "$v/ddiff.iblt"
  "$BIN/iblt_grow" write --in "$W/sets/$a.u64" --n "$n" --seed "$s" --out "$v/mine_a.iblt" 2> /dev/null
  local same_d same_a rc
  cmp -s "$v/dfull.iblt" "$v/ddiff.iblt" && same_d=true || same_d=false
  cmp -s "$v/a.iblt" "$v/mine_a.iblt" && same_a=true || same_a=false
  "$BIN/ibltseq" list -i "$v/dfull.iblt" > "$v/list.txt" 2> /dev/null; rc=$?
  emit iblt_validate "{\"i\": $i, \"pair\": \"$p\", \"ratio\": $r, \"n\": $n, \"seed\": $s, \"a_iblt_bytes\": $(stat -c %s "$v/a.iblt"), \
\"d_full_equals_d_difference_only\": $same_d, \"iblt_grow_write_equals_kmpeeler_a\": $same_a, \"kmp_list_exit\": $rc, \
\"kmp_score\": $("$BIN/kmtools" listscore $K "$W/truth/$p/amb.u64" "$W/truth/$p/bma.u64" < "$v/list.txt"), \
\"iblt_grow_peel\": $("$BIN/iblt_grow" peel --in "$v/dfull.iblt"), \"stat\": $("$BIN/iblt_grow" stat --in "$v/dfull.iblt"), \
\"wall_s\": $(( $(date +%s) - t0 ))}"
  rm -rf "$v"
}

# one km-peeler trial on a difference (ambtxt, bmatxt); its d.iblt is the whole-set one (iblt_validate)
kmp_once() { # label set ratio n seed
  local l=$1 set=$2 r=$3 n=$4 s=$5 d=$W/tr/$1_$5 sd=$W/diffsets/$2 rc t0 tl
  mkdir -p "$d"
  "$BIN/ibltseq" build -i "$sd/amb.txt" -n "$n" -r 3 -s "$s" -l $K -o "$d/a.iblt"
  "$BIN/ibltseq" build -i "$sd/bma.txt" -n "$n" -r 3 -s "$s" -l $K -o "$d/b.iblt"
  "$BIN/ibltseq" diff -i "$d/a.iblt" -j "$d/b.iblt" -o "$d/d.iblt"
  t0=$(date +%s.%N)
  timeout 1800 "$BIN/ibltseq" list -i "$d/d.iblt" > "$d/list.txt" 2> /dev/null; rc=$?
  tl=$(awk -v a="$t0" -v b="$(date +%s.%N)" 'BEGIN { printf "%.3f", b - a }')
  emit kmp_trial "{\"label\": \"$l\", \"set\": \"$set\", \"ratio\": $r, \"n\": $n, \"seed\": $s, \"D\": $(cat "$sd/D"), \
\"iblt_bytes\": $(stat -c %s "$d/a.iblt"), \"list_exit\": $rc, \"wall_list\": $tl, \
\"score\": $("$BIN/kmtools" listscore $K "$sd/amb.u64" "$sd/bma.u64" < "$d/list.txt"), \"stat\": $("$BIN/iblt_grow" stat --in "$d/d.iblt")}"
  rm -rf "$d"
}

kmp_trials() {
  mkdir -p "$W/tr" "$W/diffsets"
  local p sz sd
  for p in giga draft dvg; do
    sd=$W/diffsets/$p; mkdir -p "$sd"
    cp "$W/truth/$p/amb.u64" "$W/truth/$p/bma.u64" "$sd/"; true_d $p > "$sd/D"
  done
  for sz in 1000 10000 100000 1000000; do
    [ "$sz" -lt "$(true_d giga)" ] || continue
    sd=$W/diffsets/giga_sub$sz; mkdir -p "$sd"
    emit subset "$(python3 recon/round3_plan.py subset "$W/truth/giga/amb.u64" "$W/truth/giga/bma.u64" $sz $sz "$sd/amb.u64" "$sd/bma.u64")"
    echo $sz > "$sd/D"
  done
  for sd in "$W"/diffsets/*; do textof "$sd/amb.u64" "$sd/amb.txt"; textof "$sd/bma.u64" "$sd/bma.txt"; done
  # (set, ratio, trials)
  local plan="giga 0.95 30
giga 1.0 100
giga 1.05 300
draft 1.0 100
draft 1.05 300
dvg 1.0 100
dvg 1.05 100
giga_sub1000 1.0 100
giga_sub1000 1.05 100
giga_sub10000 1.0 100
giga_sub10000 1.05 100
giga_sub100000 1.0 100
giga_sub100000 1.05 100
giga_sub1000000 1.0 100
giga_sub1000000 1.05 100"
  plan=$(echo "$plan" | while read -r set r t; do
    [ -d "$W/diffsets/$set" ] && echo "$set $r $(( t / PLAN_DIV > 0 ? t / PLAN_DIV : 1 ))"; done)
  local total
  total=$(echo "$plan" | awk '{ s += $3 } END { print s }')
  python3 recon/round3_plan.py seeds "$total" > "$W/tr/seeds"
  echo "$plan" | while read -r set r t; do
    D=$(cat "$W/diffsets/$set/D")
    for i in $(seq 1 "$t"); do echo "kmp_${set}_r$r $set $r $(scale "$r" "$D")"; done
  done | paste -d ' ' - "$W/tr/seeds" | xargs -P "$PKMP" -L 1 bash "$0" kmp_once
}

iblt_proto() {
  local p
  for p in giga draft; do
    "$BIN/iblt_grow" proto --amb "$W/truth/$p/amb.u64" --bma "$W/truth/$p/bma.u64" --proto fixed --f 0.95,1.0,1.05 \
      --reps $REPS_FIXED --threads "$SLOTS" --label "${p}_fixed" | while read -r l; do emit iblt_proto "$l"; done
    for pr in resend subtract; do
      "$BIN/iblt_grow" proto --amb "$W/truth/$p/amb.u64" --bma "$W/truth/$p/bma.u64" --proto $pr \
        --f 0.125,0.25,0.5,0.95,1 --reps $REPS_PROTO --threads "$SLOTS" --max-rounds 8 --label "${p}_$pr" \
        | while read -r l; do emit iblt_proto "$l"; done
    done
    "$BIN/iblt_grow" proto --amb "$W/truth/$p/amb.u64" --bma "$W/truth/$p/bma.u64" --proto grow \
      --f 0.125,0.25,0.5,0.95,1 --reps $REPS_PROTO --threads "$SLOTS" --max-rounds 6 --label "${p}_grow" \
      | while read -r l; do emit iblt_proto "$l"; done
  done
}

# extended syncmers, k = 31, z = 11 (2k - z = 51), as recon/kmpeeler_run.py ext ran them on 2026-11-10
ext_set() { # genome
  local g=$1 e=$W/ext
  "$BIN/ibltseq_ext" syncmers -k $K -m 11 -g 20 -i "$W/fa/$g.fa" > "$e/$g.raw"
  LC_ALL=C sort -u -S 8G -T "$e" "$e/$g.raw" | tr a-z A-Z | LC_ALL=C sort -u -S 8G -T "$e" > "$e/$g.ext.txt"
  rm -f "$e/$g.raw"
}
ext_trial() { # pair i n seed
  local p=$1 i=$2 n=$3 s=$4 a b e=$W/ext d=$W/ext/t_$1_$2 rc t0 tl
  a=$(pair_a "$p"); b=$(pair_b "$p")
  mkdir -p "$d"
  t0=$(date +%s.%N)
  "$BIN/ibltseq_ext" build -i "$e/$a.ext.txt" -n "$n" -r 3 -s "$s" -l 51 -o "$d/a.iblt"
  "$BIN/ibltseq_ext" build -i "$e/$b.ext.txt" -n "$n" -r 3 -s "$s" -l 51 -o "$d/b.iblt"
  "$BIN/ibltseq_ext" diff -i "$d/a.iblt" -j "$d/b.iblt" -o "$d/d.iblt"
  timeout 7200 "$BIN/ibltseq_ext" list -i "$d/d.iblt" > "$d/list.txt" 2> "$d/list.err"; rc=$?
  tl=$(awk -v a="$t0" -v b="$(date +%s.%N)" 'BEGIN { printf "%.1f", b - a }')
  emit ext_trial "{\"pair\": \"$p\", \"trial\": $i, \"n\": $n, \"seed\": $s, \"D_ext\": $(cat "$e/D_$p"), \
\"iblt_bytes\": $(stat -c %s "$d/a.iblt"), \"list_exit\": $rc, \"wall_s\": $tl, \"listed\": $(wc -l < "$d/list.txt"), \
\"score\": $("$BIN/kmtools" extscore $K "$W/truth/$p/amb.u64" "$W/truth/$p/bma.u64" "$d/miss.txt" < "$d/list.txt")}"
  md5sum < "$d/miss.txt" | cut -c1-32 > "$d/miss.md5"
  mkdir -p "$OUT/ext_${p}_trial$i"
  cp "$d/miss.txt" "$d/miss.md5" "$OUT/ext_${p}_trial$i/"
  if [ "$i" = 1 ]; then cp "$d/list.txt" "$e/list_$p.txt"; fi
  rm -rf "$d"
}
ext() {
  mkdir -p "$W/ext"
  local g p a b D n i
  for g in $GENOMES; do ext_set "$g" & done
  wait
  for p in giga draft; do
    a=$(pair_a $p); b=$(pair_b $p)
    LC_ALL=C comm -3 "$W/ext/$a.ext.txt" "$W/ext/$b.ext.txt" | tr -d '\t' | LC_ALL=C sort > "$W/ext/diff_$p.txt"
    wc -l < "$W/ext/diff_$p.txt" > "$W/ext/D_$p"
    emit ext_set "{\"pair\": \"$p\", \"ext_A\": $(wc -l < "$W/ext/$a.ext.txt"), \"ext_B\": $(wc -l < "$W/ext/$b.ext.txt"), \
\"D_ext\": $(cat "$W/ext/D_$p"), \"ext_diff_with_non_ACGT\": $(grep -c '[^ACGT]' "$W/ext/diff_$p.txt")}"
  done
  # giga: 10 trials with independent seeds, 5 at a time
  D=$(cat "$W/ext/D_giga"); n=$(scale 1.05 "$D")
  python3 recon/round3_plan.py seeds $((EXT_TRIALS + 3)) > "$W/ext/seeds"
  for i in $(seq 1 "$EXT_TRIALS"); do echo "giga $i $n $(sed -n ${i}p "$W/ext/seeds")"; done | xargs -P "$PEXT" -L 1 bash "$0" ext_trial
  ext_trace giga
  # draft vs final: the run that died out of memory on 2026-11-10; 3 trials, one at a time, 2 h listing cap each
  D=$(cat "$W/ext/D_draft"); n=$(scale 1.05 "$D")
  for i in 1 2 3; do ext_trial draft $i "$n" "$(sed -n $((EXT_TRIALS + i))p "$W/ext/seeds")"; done
}
ext_trace() { # pair: where the missed k-mers are, and which differing syncmers the listing left out
  local p=$1 a b e=$W/ext t=$OUT/ext_trace_$1
  a=$(pair_a "$p"); b=$(pair_b "$p")
  mkdir -p "$t"
  emit ext_miss_md5 "{\"pair\": \"$p\", \"distinct_miss_sets\": $(cat "$OUT"/ext_${p}_trial*/miss.md5 | sort -u | wc -l)}"
  cut -f1 "$OUT/ext_${p}_trial1/miss.txt" | head -200 > "$t/miss_kmers.txt"   # 7 on 2026-11-10; capped for the greps
  "$BIN/kmtools" locate $K "$t/miss_kmers.txt" < "$W/fa/$a.fa" > "$t/locate_$a.tsv" 2> "$t/locate_$a.json"
  "$BIN/kmtools" locate $K "$t/miss_kmers.txt" < "$W/fa/$b.fa" > "$t/locate_$b.tsv" 2> "$t/locate_$b.json"
  # each missed k-mer, both strands, in the extended-syncmer sets of each genome
  local km rc
  while read -r km; do
    rc=$(echo "$km" | rev | tr ACGT TGCA)
    printf '%s\t%s\t%s\t%s\n' "$km" "$(grep -c -F -e "$km" -e "$rc" "$e/$a.ext.txt")" \
      "$(grep -c -F -e "$km" -e "$rc" "$e/$b.ext.txt")" "$(grep -c -F -e "$km" -e "$rc" "$e/diff_$p.txt")"
    grep -F -e "$km" -e "$rc" "$e/diff_$p.txt" | sed "s/^/$km\t/" >> "$t/miss_in_diff_syncmers.tsv"
  done < "$t/miss_kmers.txt" > "$t/miss_in_ext_sets.tsv"
  # the differing extended syncmers that trial 1's listing did not return
  cut -c3- "$e/list_$p.txt" | LC_ALL=C sort > "$e/listed_$p.txt"
  LC_ALL=C comm -23 "$e/diff_$p.txt" "$e/listed_$p.txt" > "$t/unlisted_syncmers.txt"
  LC_ALL=C comm -13 "$e/diff_$p.txt" "$e/listed_$p.txt" > "$t/listed_not_in_diff.txt"
  emit ext_trace "{\"pair\": \"$p\", \"missed\": $(wc -l < "$t/miss_kmers.txt"), \"unlisted\": $(wc -l < "$t/unlisted_syncmers.txt"), \
\"unlisted_with_non_ACGT\": $(grep -c '[^ACGT]' "$t/unlisted_syncmers.txt"), \"listed_not_in_diff\": $(wc -l < "$t/listed_not_in_diff.txt"), \
\"locate_A\": $(cat "$t/locate_$a.json"), \"locate_B\": $(cat "$t/locate_$b.json")}"
}

tally() {
  "$BIN/kmtools" tally $K "$W/truth/giga/amb.u64" < "$W/fa/Cannbio-2.fa" > "$OUT/tally_giga_A_Cannbio-2.tsv" 2> "$W/t1.json" &
  "$BIN/kmtools" tally $K "$W/truth/giga/bma.u64" < "$W/fa/Cannbio2_GigaDB.fa" > "$OUT/tally_giga_B_Cannbio2_GigaDB.tsv" 2> "$W/t2.json" &
  "$BIN/kmtools" tally $K "$W/truth/draft/bma.u64" < "$W/fa/Cannbio-2_draft.fa" > "$OUT/tally_draft_B_Cannbio-2_draft.tsv" 2> "$W/t3.json" &
  wait
  emit tally "{\"giga_A\": $(cat "$W/t1.json"), \"giga_B\": $(cat "$W/t2.json"), \"draft_B\": $(cat "$W/t3.json")}"
}

cmd=$1; shift
case $cmd in
  build|setup|kp_seeds|kp_sweep|k63|iblt_validate|kmp_trials|iblt_proto|ext|tally|kp_one|k63_one|iblt_validate_one|kmp_once|ext_trial)
    "$cmd" "$@" ;;
  *) echo "unknown phase $cmd" >&2; exit 2 ;;
esac
