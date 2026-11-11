#!/usr/bin/env bash
# Every test suite and mutant check in the repository, in one run, with its full output.
#
#   bash tests/run_all.sh [OUT]          # OUT defaults to stdout
#   PYTHON=/path/to/python bash tests/run_all.sh results/x.txt
#
# C: hash/ (make test test_pack test_e; make mutation mutation_pack mutation_e) and ideas/{algebra,eisenstein,kp5_kp11}
# (make test; make mutation). Python: pytest (tests/), the two ctypes ports' self-tests, ideas/algebra/roots.py.
# Every mutant must be killed; the script exits nonzero if any suite fails or any mutant survives.
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
PY=${PYTHON:-python3}
OUT=${1:-/dev/stdout}
fail=0
run() { # label, dir, command...
  local label=$1 dir=$2; shift 2
  echo "===== $label ($dir: $*)"
  (cd "$ROOT/$dir" && "$@") 2>&1 | grep -v '^gcc \|^g++ '
  local rc=${PIPESTATUS[0]}
  echo "----- $label: exit $rc"
  [ "$rc" -eq 0 ] || fail=1
}
{
  echo "run_all.sh at $(git -C "$ROOT" rev-parse --short HEAD)$(git -C "$ROOT" diff --quiet HEAD -- hash ideas recon prototype tests || echo '+uncommitted'), $(date -u +%FT%TZ), $(uname -n)"
  run "hash build" hash make test_kp7 test_kp7m31 test_kp7pack test_kp7e
  run "hash tests" hash make test test_pack test_e
  run "hash mutants" hash make mutation mutation_pack mutation_e
  for d in algebra eisenstein kp5_kp11; do
    run "ideas/$d tests" "ideas/$d" make test
    run "ideas/$d mutants" "ideas/$d" make mutation
  done
  run "pytest" . "$PY" -m pytest -q
  run "kp7e_hash self-test" . "$PY" recon/kp7e_hash.py
  run "kp7hash self-test" . "$PY" recon/kp7hash.py
  run "roots.py" ideas/algebra "$PY" roots.py
  echo "===== overall: $([ $fail -eq 0 ] && echo PASS || echo FAIL)"
} > "$OUT"
exit $fail
