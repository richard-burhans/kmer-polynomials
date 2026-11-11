#!/usr/bin/env bash
# Run hash/nteq on a genome and write its output with provenance (ntHash commit, compiler, flags, host, time).
#   bash hash/run_nteq.sh <fasta> <out.txt> [k=31]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
FA=$1 OUT=$2 K=${3:-31}
NTHASH=${NTHASH:-/home/rcb112/project/claude/bench/data/kmer_poly_work/ntHash}
make -s -C "$HERE" nteq NTHASH="$NTHASH"
{
  echo "# hash/nteq $FA k=$K"
  echo "# repo $(git -C "$HERE" rev-parse --short HEAD)$(git -C "$HERE" diff --quiet HEAD -- . || echo '+uncommitted'); ntHash $(git -C "$NTHASH" rev-parse HEAD 2>/dev/null || cat "$NTHASH/.git_commit")"
  echo "# $(g++ --version | head -1); flags $(grep -m1 '^OPT' "$HERE/Makefile"); host $(uname -n); $(date -u +%FT%TZ)"
  echo "# fasta sha256 $(sha256sum "$FA" | cut -c1-16)"
  /usr/bin/time -f "# wall %e s, max RSS %M KB, exit %x" "$HERE/nteq" "$FA" "$K" 2>&1
} > "$OUT"
tail -2 "$OUT"
