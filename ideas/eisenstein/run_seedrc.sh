#!/usr/bin/env bash
# Build and run ideas/eisenstein/seedrc, writing its output (stdout and ntHash's stderr warnings) with provenance.
#   bash ideas/eisenstein/run_seedrc.sh <out.tsv> [n=200000] [seed=7]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=$1
NTHASH=${NTHASH:-/home/rcb112/project/claude/bench/data/kmer_poly_work/ntHash}
make -s -C "$HERE" seedrc NTHASH="$NTHASH"
{
  echo "# ideas/eisenstein/seedrc; repo $(git -C "$HERE" rev-parse --short HEAD); ntHash $(git -C "$NTHASH" rev-parse --short HEAD); host $(uname -n); $(date -u +%FT%TZ)"
  "$HERE/seedrc" "${2:-200000}" "${3:-7}" 2> >(sed 's/^/# stderr: /' >&2)
} > "$OUT" 2>&1
grep -v '^#' "$OUT"
