"""Exact enumeration of the strand-symmetric candidates F*R, F+R and (F+lam)(R+lam) over the RING, not a field.

    python3 ideas/algebra/fr_exact.py [k ...]        # default: 7 8; exits nonzero if a structural check fails

Each k-mer s gives F_s(x) = sum_t v(s_t) x^t and R_s(x) = F_{rc(s)}(x), polynomials with coefficients in
  fourth-root  Z[i]:  A = 1, C = i, G = -i, T = -1   (prototype/kp7_symmetric_check.py, hash/kp7.c, hash/kp7m31.c)
  Eisenstein   Z[w]:  A = 1, C = w, G = -w, T = -1, w^2 = -1 - w   (hash/kp7e.c)
Two k-mers have the same F*R (or F+R) for EVERY x exactly when the product (or sum) polynomials are equal, so counting
distinct coefficient tuples gives the structural count, with no random point and no chance collision (round-1
adjudication L4 F2 and L1 F9: the paper's caption ran F*R over Z[i]/(2^31 - 1) at one random x). (F+lam)(R+lam) =
F*R + lam (F+R) + lam^2, so as a polynomial in (x, lam) it is determined by the pair (F*R, F+R).

Also counted: double-stranded classes {s, rc(s)}, and orbits under the group {identity, reverse, complement, reverse
complement}, which F*R cannot separate: R_s = -F_{rev s} and R_{rev s} = -F_s, so reversal leaves F*R unchanged, and
complement negates both factors. The excess of orbits over distinct F*R values is the number of extra merges by k-mers with the same autocorrelation
that are not related by reversal (homometric pairs). One example pair is printed.
"""

from __future__ import annotations

import itertools
import sys

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
# ring elements as integer pairs (a, b) = a + b*u, u = i or w
ENC = {
    "fourth-root Z[i]": ({"A": (1, 0), "C": (0, 1), "G": (0, -1), "T": (-1, 0)}, "i"),
    "Eisenstein Z[w]": ({"A": (1, 0), "C": (0, 1), "G": (0, -1), "T": (-1, 0)}, "w"),
}


def rmul(a, b, u):
    """(a0 + a1 u)(b0 + b1 u) with u^2 = -1 (i) or u^2 = -1 - u (w)."""
    p0, p1, p2 = a[0] * b[0], a[0] * b[1] + a[1] * b[0], a[1] * b[1]
    return (p0 - p2, p1) if u == "i" else (p0 - p2, p1 - p2)


def poly(s, val):
    return tuple(val[c] for c in s)


def pmul(f, g, u):
    out = [(0, 0)] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        for j, b in enumerate(g):
            m = rmul(a, b, u)
            out[i + j] = (out[i + j][0] + m[0], out[i + j][1] + m[1])
    return tuple(out)


def padd(f, g):
    return tuple((a[0] + b[0], a[1] + b[1]) for a, b in zip(f, g))


def rc(s):
    return "".join(COMP[c] for c in reversed(s))


def comp(s):
    return "".join(COMP[c] for c in s)


def main() -> int:
    ks = [int(a) for a in sys.argv[1:]] or [7, 8]
    fails = []
    for k in ks:
        kmers = ["".join(t) for t in itertools.product("ACGT", repeat=k)]
        classes = len({min(s, rc(s)) for s in kmers})
        orbits = len({min(s, s[::-1], comp(s), rc(s)) for s in kmers})
        burnside = (4 ** k + 4 ** ((k + 1) // 2) + 0 + (4 ** (k // 2) if k % 2 == 0 else 0)) // 4
        print(f"k={k}: {len(kmers)} k-mers, {classes} double-stranded classes, {orbits} orbits under "
              f"{{id, reverse, complement, reverse complement}} (Burnside {burnside})")
        if orbits != burnside:
            fails.append(f"k={k}: orbit count {orbits} != Burnside {burnside}")
        for name, (val, u) in ENC.items():
            prod, summ, pair = {}, set(), set()
            for s in kmers:
                f, r = poly(s, val), poly(rc(s), val)
                fr, fpr = pmul(f, r, u), padd(f, r)
                prod.setdefault(fr, []).append(s)
                summ.add(fpr)
                pair.add((fr, fpr))
            n_fr = len(prod)
            example = ""
            for group in prod.values():
                reps = {min(s, s[::-1], comp(s), rc(s)) for s in group}
                if len(reps) > 1:
                    a, b = sorted(group)[0], next(s for s in sorted(group) if min(s, s[::-1], comp(s), rc(s))
                                                   != min(group[0], group[0][::-1], comp(group[0]), rc(group[0])))
                    example = f"; e.g. {a} and {b} share F*R but are not related by reversal or complement"
                    break
            print(f"  {name}: F*R distinct {n_fr} (orbits - F*R = {orbits - n_fr} extra merges){example}")
            print(f"  {name}: F+R distinct {len(summ)}; (F+lam)(R+lam) as a polynomial in (x, lam) distinct {len(pair)}"
                  f" (classes {classes})")
            if len(pair) != classes:
                fails.append(f"k={k} {name}: (F+lam)(R+lam) separates {len(pair)} != {classes} classes")
            if n_fr > orbits:
                fails.append(f"k={k} {name}: F*R separates more than the orbits ({n_fr} > {orbits})")
    if fails:
        print("FAILED:", *fails, sep="\n  ")
        return 1
    print("all structural checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
