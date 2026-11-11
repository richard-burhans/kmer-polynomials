#!/usr/bin/env bash
# Tests for ideas/kp5_kp11: hash vs packed-code reference vs Python brute force (KP-5), and fft vs direct vs brute force
# (KP-11). Exit 0 iff everything matches.
set -euo pipefail
cd "$(dirname "$0")"
T=$(mktemp -d "${TMPDIR:-/tmp}/kp511.XXXXXX"); trap 'rm -rf "$T"' EXIT
fail=0
python3 gen_test.py fasta "$T/t.fa" 11
field() { sed -n "s/.* $1=\([^ ]*\).*/\1/p"; }
# KP-5: palindromes (even k) and inverted repeats with arms on both sides of 32 and at the 64 limit
for spec in "palin 12" "palin 8" "inv 12 40" "inv 20 40" "inv 33 40" "inv 40 40" "inv 64 40" "inv 6 10"; do
  set -- $spec
  h=$(./ir "$T/t.fa" - $spec hash); p=$(./ir "$T/t.fa" - $spec pack)
  hc=$(echo "$h" | field count); pc=$(echo "$p" | field count)
  hs=$(echo "$h" | field checksum); ps=$(echo "$p" | field checksum)
  b=$(python3 gen_test.py "$1" "$T/t.fa" "${@:2}")
  if [ "$hc" = "$pc" ] && [ "$hs" = "$ps" ] && [ "$hc" = "$b" ]; then echo "ok   $spec: $hc"; else echo "FAIL $spec: hash $hc/$hs pack $pc/$ps brute $b"; fail=1; fi
done
# KP-11: fft == direct over all 4^k targets, and both == brute force on sampled targets
for kd in "6 0" "6 1" "6 2" "6 3" "7 2" "9 3"; do
  set -- $kd
  out=$(./mismatch "$T/t.fa" - $1 $2 both)
  [ "$(echo "$out" | field identical)" = yes ] && echo "ok   mismatch k=$1 d=$2 fft == direct" || { echo "FAIL mismatch k=$1 d=$2: $out"; fail=1; }
done
# brute force on sampled targets, against direct (which fft equals)
python3 gen_test.py mismatch "$T/t.fa" 7 2 25 > "$T/bf.txt"
if ./mismatch "$T/t.fa" - 7 2 check < "$T/bf.txt" > "$T/cmp.txt"; then echo "ok   mismatch k=7 d=2 brute force on 25 targets"; else echo "FAIL brute force"; fail=1; fi
exit $fail
