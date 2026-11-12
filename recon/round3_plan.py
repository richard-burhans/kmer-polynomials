"""Trial plans for the round-3 owed km-peeler runs (recon/round3_owed_job.sh, phase kmp_trials). Standard library only
(the laila job image has python3 but no pip).

    python3 recon/round3_plan.py seeds <n>                      n km-peeler root seeds, one per line
    python3 recon/round3_plan.py subset <amb.u64> <bma.u64> <size> <seed> <out_amb> <out_bma>

seeds: 32-bit roots drawn from os.urandom. km-peeler derives its r = 3 hash functions as root, root + 1, root + 2
(ibflib.c:27-29), so two roots closer than 3 share a hash function; the draw rejects any root whose triple overlaps
an earlier one, so no two trials share a hash function (round-2 adjudications L1 F4, L4 F3, L7 F6).

subset: a uniformly random subset of `size` elements of the difference A\\B U B\\A, each kept on its own side, for the
held-out sizing check at other difference sizes. A table built from a subset is the difference table of two sets
differing by exactly that subset (linearity, recon/iblt_grow.c), so this is a smaller real difference, not a model.
"""

from __future__ import annotations

import array
import os
import random
import sys


def seeds(n: int) -> list[int]:
    used: set[int] = set()
    out: list[int] = []
    while len(out) < n:
        r = int.from_bytes(os.urandom(4), "little")
        if r > 0xFFFFFFFF - 2 or any((r + j) in used for j in range(3)):
            continue
        used.update(r + j for j in range(3))
        out.append(r)
    return out


def read(path: str) -> array.array:
    a = array.array("Q")
    with open(path, "rb") as fh:
        a.frombytes(fh.read())
    return a


def subset(amb: str, bma: str, size: int, seed: int, out_a: str, out_b: str) -> None:
    a, b = read(amb), read(bma)
    pick = sorted(random.Random(seed).sample(range(len(a) + len(b)), size))
    sa = array.array("Q", (a[i] for i in pick if i < len(a)))
    sb = array.array("Q", (b[i - len(a)] for i in pick if i >= len(a)))
    with open(out_a, "wb") as fh:
        sa.tofile(fh)
    with open(out_b, "wb") as fh:
        sb.tofile(fh)
    print(f'{{"subset": {size}, "seed": {seed}, "A_minus_B": {len(sa)}, "B_minus_A": {len(sb)}}}')


def main() -> int:
    if len(sys.argv) == 3 and sys.argv[1] == "seeds":
        print("\n".join(map(str, seeds(int(sys.argv[2])))))
        return 0
    if len(sys.argv) == 8 and sys.argv[1] == "subset":
        subset(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), sys.argv[6], sys.argv[7])
        return 0
    print(__doc__, file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
