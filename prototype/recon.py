"""KP-8/KP-9/KP-10: partitioned characteristic-polynomial set reconciliation over F_{q^2}.

A set S of field elements (k-mer hashes) has characteristic polynomial chi_S(t) = prod_{s in S} (t - s).
For two sets A and B, chi_A / chi_B = chi_{A\\B} / chi_{B\\A} (Minsky, Trachtenberg & Zippel 2003). Evaluating both
at m shared points and solving a linear system recovers a numerator and a denominator whose roots include the two
differences, whenever |A\\B| + |B\\A| <= m.

Partitioning (Minsky & Trachtenberg): an independent hash splits the k-mer space into 2^bits buckets, so each bucket's
difference is small. A bucket whose difference exceeds m is detected at extra check points and split in two using the
next bit of the bucket hash, at the cost of one more round.

Decoding needs no root-finding: each side evaluates the recovered polynomial at its own bucket's hashes and keeps the
zeros. A spurious zero would need a random common factor of the recovered pair to vanish at one of our hashes, which
has probability about deg / q^2 per element, and the final check at independent points would catch it.

KP-10 power sums S_j = sum_s s^j (j = 1..m) are the linear, mergeable form of the same sketch. Newton's identities
turn the difference of two power-sum sketches into the series of chi_{A\\B}^rev / chi_{B\\A}^rev, and a Pade
approximation recovers the pair.
"""

from __future__ import annotations

from dataclasses import dataclass, field as dc_field

import numpy as np

from . import field as fq


# ── linear algebra over F_{q^2} ────────────────────────────────────────────────────────────────────────────────────

def solve(Mre, Mim, rre, rim):
    """Solve M x = r over F_{q^2} by Gauss-Jordan. Free variables are set to 0.

    Returns (x_re, x_im) or None when the system is inconsistent.
    """
    rows, cols = Mre.shape
    A = (np.concatenate([Mre, rre[:, None]], axis=1) % fq.Q, np.concatenate([Mim, rim[:, None]], axis=1) % fq.Q)
    A = (A[0].copy(), A[1].copy())
    pivots = []
    r = 0
    for c in range(cols):
        if r >= rows:
            break
        nz = np.nonzero((A[0][r:, c] != 0) | (A[1][r:, c] != 0))[0]
        if nz.size == 0:
            continue
        p = r + nz[0]
        if p != r:
            A[0][[r, p]] = A[0][[p, r]]
            A[1][[r, p]] = A[1][[p, r]]
        pinv = fq.inv((A[0][r, c], A[1][r, c]))
        A = (A[0].copy(), A[1].copy())
        row = fq.mul((A[0][r], A[1][r]), (np.full(cols + 1, pinv[0]), np.full(cols + 1, pinv[1])))
        A[0][r], A[1][r] = row
        factor = (A[0][:, c].copy(), A[1][:, c].copy())
        factor[0][r] = 0
        factor[1][r] = 0
        upd = fq.mul((factor[0][:, None], factor[1][:, None]), (row[0][None, :], row[1][None, :]))
        A = fq.sub(A, upd)
        pivots.append(c)
        r += 1
    # inconsistent if a zero row has a nonzero right-hand side
    if r < rows and np.any((A[0][r:, cols] != 0) | (A[1][r:, cols] != 0)):
        return None
    xre = np.zeros(cols, np.int64)
    xim = np.zeros(cols, np.int64)
    for i, c in enumerate(pivots):
        xre[c] = A[0][i, cols]
        xim[c] = A[1][i, cols]
    return xre, xim


def powers(z, n):
    """[z^0, ..., z^{n-1}] for a vector z of points: shape (len(z), n)."""
    re = np.ones((z[0].size, n), np.int64)
    im = np.zeros((z[0].size, n), np.int64)
    cur = (np.ones(z[0].size, np.int64), np.zeros(z[0].size, np.int64))
    for j in range(1, n):
        cur = fq.mul(cur, z)
        re[:, j], im[:, j] = cur
    return re, im


def poly_at(coeffs, x):
    """Evaluate a polynomial given as (re, im) coefficient arrays (low degree first) at points x."""
    acc = (np.zeros_like(x[0]), np.zeros_like(x[1]))
    for j in range(coeffs[0].size - 1, -1, -1):
        acc = fq.add(fq.mul(acc, x), (np.full_like(x[0], coeffs[0][j]), np.full_like(x[1], coeffs[1][j])))
    return acc


def chi_at(elems, z):
    """chi_S(z) = prod (z - s) for each point in z (vectors); elems is a small element array."""
    if elems[0].size == 0:
        return (np.ones(z[0].size, np.int64), np.zeros(z[0].size, np.int64))
    terms = fq.sub((z[0][:, None], z[1][:, None]), (elems[0][None, :], elems[1][None, :]))
    return fq.prod_rows(terms)


# ── sketches ───────────────────────────────────────────────────────────────────────────────────────────────────────

@dataclass
class Side:
    """One party's set: element hashes (re, im) and 64-bit bucket keys, sorted by key."""

    h: tuple
    key: np.ndarray
    codes: np.ndarray

    @classmethod
    def build(cls, codes, h, key):
        order = np.argsort(key, kind="stable")
        return cls((h[0][order], h[1][order]), key[order], codes[order])

    def bucket_slices(self, bits: int, buckets: np.ndarray):
        """(start, stop) index ranges of the given bucket ids at a resolution of `bits` bits."""
        b = (self.key >> np.uint64(64 - bits)).astype(np.int64) if bits else np.zeros(self.key.size, np.int64)
        lo = np.searchsorted(b, buckets, side="left")
        hi = np.searchsorted(b, buckets, side="right")
        return lo, hi


def padded(side: Side, lo, hi, pad):
    """Bucket elements as a padded (nb, L) element matrix; padding positions hold `pad` (re, im)."""
    sizes = hi - lo
    L = max(int(sizes.max()) if sizes.size else 1, 1)
    idx = lo[:, None] + np.arange(L)[None, :]
    valid = idx < hi[:, None]
    idx = np.where(valid, idx, 0)
    re = np.where(valid, side.h[0][idx], pad[0])
    im = np.where(valid, side.h[1][idx], pad[1])
    return (re, im), sizes


def eval_sketch(side: Side, bits: int, buckets: np.ndarray, z, chunk: int = 256):
    """chi_bucket(z_i) for every bucket and point: (nb, npts) arrays, plus bucket sizes."""
    lo, hi = side.bucket_slices(bits, buckets)
    out_re = np.empty((buckets.size, z[0].size), np.int64)
    out_im = np.empty_like(out_re)
    for s in range(0, buckets.size, chunk):
        (Hre, Him), _ = padded(side, lo[s:s + chunk], hi[s:s + chunk], (0, 0))
        valid = np.arange(Hre.shape[1])[None, :] < (hi[s:s + chunk] - lo[s:s + chunk])[:, None]
        for i in range(z[0].size):
            t = fq.sub((np.full_like(Hre, z[0][i]), np.full_like(Him, z[1][i])), (Hre, Him))
            t = (np.where(valid, t[0], 1), np.where(valid, t[1], 0))
            out_re[s:s + chunk, i], out_im[s:s + chunk, i] = fq.prod_rows(t)
    return (out_re, out_im), hi - lo


def power_sketch(side: Side, bits: int, buckets: np.ndarray, m: int, chunk: int = 256):
    """S_j = sum_s s^j for j = 1..m per bucket: (nb, m) arrays, plus bucket sizes."""
    lo, hi = side.bucket_slices(bits, buckets)
    out_re = np.empty((buckets.size, m), np.int64)
    out_im = np.empty_like(out_re)
    for s in range(0, buckets.size, chunk):
        H, _ = padded(side, lo[s:s + chunk], hi[s:s + chunk], (0, 0))
        pw = H
        for j in range(m):
            out_re[s:s + chunk, j] = pw[0].sum(axis=1) % fq.Q
            out_im[s:s + chunk, j] = pw[1].sum(axis=1) % fq.Q
            pw = fq.mul(pw, H)
    return (out_re, out_im), hi - lo


# ── decoding one bucket ────────────────────────────────────────────────────────────────────────────────────────────

def decode_eval(fA, fB, nA, nB, z, m, ownA, ownB):
    """Recover (flagsA, flagsB) for one bucket from evaluation sketches, or None on overflow.

    fA, fB: chi values at the m + c points (first m solve, the rest check). ownA, ownB: each side's bucket hashes.
    """
    f = fq.div(fA, fB)
    delta = int(nA) - int(nB)
    mbar = m if (m - delta) % 2 == 0 else m - 1
    if abs(delta) > mbar:
        return None
    dA, dB = (mbar + delta) // 2, (mbar - delta) // 2
    zs, fs = (z[0][:m], z[1][:m]), (f[0][:m], f[1][:m])
    Zp = powers(zs, max(dA, dB) + 1)
    # unknowns: p_0..p_{dA-1}, q_0..q_{dB-1};  sum p_j z^j - f sum q_j z^j = f z^dB - z^dA
    colsP = (Zp[0][:, :dA], Zp[1][:, :dA])
    colsQ = fq.neg(fq.mul((fs[0][:, None], fs[1][:, None]), (Zp[0][:, :dB], Zp[1][:, :dB])))
    Mre = np.concatenate([colsP[0], colsQ[0]], axis=1)
    Mim = np.concatenate([colsP[1], colsQ[1]], axis=1)
    rhs = fq.sub(fq.mul(fs, (Zp[0][:, dB], Zp[1][:, dB])), (Zp[0][:, dA], Zp[1][:, dA]))
    sol = solve(Mre, Mim, rhs[0], rhs[1])
    if sol is None:
        return None
    P = (np.append(sol[0][:dA], 1), np.append(sol[1][:dA], 0))
    Qp = (np.append(sol[0][dA:], 1), np.append(sol[1][dA:], 0))
    # check points: P/Q must equal f there
    zc, fc = (z[0][m:], z[1][m:]), (f[0][m:], f[1][m:])
    if not np.all(fq.eq(fq.mul(poly_at(P, zc), fq.inv(poly_at(Qp, zc))), fc)):
        return None
    flagsA = fq.is_zero(poly_at(P, ownA))
    flagsB = fq.is_zero(poly_at(Qp, ownB))
    # final check at the independent points: the flagged sets alone must explain f
    got = fq.div(chi_at(fq.take(ownA, flagsA), zc), chi_at(fq.take(ownB, flagsB), zc))
    if int(flagsA.sum()) - int(flagsB.sum()) != delta or not np.all(fq.eq(got, fc)):
        return None
    return flagsA, flagsB


def newton_series(S, n):
    """Coefficients r_0..r_n of exp(-sum_j S_j t^j / j), i.e. of prod(1 - a t) / prod(1 - b t) for power sums S."""
    r_re = np.zeros(n + 1, np.int64)
    r_im = np.zeros(n + 1, np.int64)
    r_re[0] = 1
    for k in range(1, n + 1):
        acc = (np.int64(0), np.int64(0))
        for j in range(1, k + 1):
            acc = fq.add(acc, fq.mul((S[0][j - 1], S[1][j - 1]), (r_re[k - j], r_im[k - j])))
        kinv = fq.inv((np.int64(k), np.int64(0)))
        val = fq.neg(fq.mul(acc, kinv))
        r_re[k], r_im[k] = val
    return r_re, r_im


def decode_power(SA, SB, nA, nB, m, c, ownA, ownB):
    """Recover (flagsA, flagsB) for one bucket from power-sum sketches S_1..S_{m+c}, or None on overflow."""
    S = fq.sub(SA, SB)
    delta = int(nA) - int(nB)
    mbar = m if (m - delta) % 2 == 0 else m - 1
    if abs(delta) > mbar:
        return None
    dA, dB = (mbar + delta) // 2, (mbar - delta) // 2
    r = newton_series((S[0][:m], S[1][:m]), mbar)
    # Pade: Q_rev (q_0 = 1, degree dB) with coefficients dA+1..dA+dB of Q_rev * R equal to zero
    if dB:
        Mre = np.empty((dB, dB), np.int64)
        Mim = np.empty_like(Mre)
        for i, n in enumerate(range(dA + 1, dA + dB + 1)):
            for j in range(1, dB + 1):
                Mre[i, j - 1], Mim[i, j - 1] = r[0][n - j], r[1][n - j]
        rhs = fq.neg((r[0][dA + 1:dA + dB + 1], r[1][dA + 1:dA + dB + 1]))
        sol = solve(Mre, Mim, rhs[0], rhs[1])
        if sol is None:
            return None
        qr = (np.concatenate([[1], sol[0]]), np.concatenate([[0], sol[1]]))
    else:
        qr = (np.array([1], np.int64), np.array([0], np.int64))
    pr_re = np.zeros(dA + 1, np.int64)
    pr_im = np.zeros(dA + 1, np.int64)
    for n in range(dA + 1):
        acc = (np.int64(0), np.int64(0))
        for j in range(0, min(n, dB) + 1):
            acc = fq.add(acc, fq.mul((qr[0][j], qr[1][j]), (r[0][n - j], r[1][n - j])))
        pr_re[n], pr_im[n] = acc
    # roots of the reversed polynomials are 1/a and 1/b
    flagsA = fq.is_zero(poly_at((pr_re, pr_im), fq.inv(ownA)))
    flagsB = fq.is_zero(poly_at(qr, fq.inv(ownB)))
    # check with the extra power sums S_{m+1}..S_{m+c}: the flagged sets alone must reproduce them
    fa, fb = fq.take(ownA, flagsA), fq.take(ownB, flagsB)
    for j in range(m, m + c):
        pa = fq.F(0)
        pb = fq.F(0)
        if fa[0].size:
            x = fa
            for _ in range(j):
                x = fq.mul(x, fa)
            pa = (x[0].sum() % fq.Q, x[1].sum() % fq.Q)
        if fb[0].size:
            x = fb
            for _ in range(j):
                x = fq.mul(x, fb)
            pb = (x[0].sum() % fq.Q, x[1].sum() % fq.Q)
        if not fq.eq(fq.sub(pa, pb), (S[0][j], S[1][j])):
            return None
    if int(flagsA.sum()) - int(flagsB.sum()) != delta:
        return None
    return flagsA, flagsB


# ── the protocol ───────────────────────────────────────────────────────────────────────────────────────────────────

@dataclass
class Result:
    a_minus_b: np.ndarray
    b_minus_a: np.ndarray
    rounds: int = 0
    buckets_decoded: int = 0
    splits: int = 0
    sketch_elements: int = 0  # field elements sent (evaluations or power sums), plus one size per bucket
    log: list = dc_field(default_factory=list)


def reconcile(A: Side, B: Side, bits: int, m: int, c: int, z, method: str = "eval",
              buckets: np.ndarray | None = None, max_bits: int = 62, max_rounds: int = 12) -> Result:
    """Reconcile A and B bucket by bucket, splitting any bucket whose difference overflows m.

    A correct decoder needs a split only for the rare bucket whose difference exceeds m, so a few rounds suffice.
    `max_rounds` turns a decoder that never succeeds into a fast error instead of a 62-level descent.
    """
    pending = np.arange(1 << bits, dtype=np.int64) if buckets is None else np.asarray(buckets, np.int64)
    a_out, b_out = [], []
    res = Result(np.empty(0, np.uint64), np.empty(0, np.uint64))
    while pending.size:
        res.rounds += 1
        if method == "eval":
            skA, nA = eval_sketch(A, bits, pending, z)
            skB, nB = eval_sketch(B, bits, pending, z)
        else:
            skA, nA = power_sketch(A, bits, pending, m + c)
            skB, nB = power_sketch(B, bits, pending, m + c)
        res.sketch_elements += pending.size * (m + c + 1)
        loA, hiA = A.bucket_slices(bits, pending)
        loB, hiB = B.bucket_slices(bits, pending)
        overflow = []
        for i, bkt in enumerate(pending):
            ownA = (A.h[0][loA[i]:hiA[i]], A.h[1][loA[i]:hiA[i]])
            ownB = (B.h[0][loB[i]:hiB[i]], B.h[1][loB[i]:hiB[i]])
            fa = (skA[0][i], skA[1][i])
            fb = (skB[0][i], skB[1][i])
            if method == "eval":
                out = decode_eval(fa, fb, nA[i], nB[i], z, m, ownA, ownB)
            else:
                out = decode_power(fa, fb, nA[i], nB[i], m, c, ownA, ownB)
            if out is None:
                overflow.append(bkt)
                continue
            fA, fB = out
            a_out.append(A.codes[loA[i]:hiA[i]][fA])
            b_out.append(B.codes[loB[i]:hiB[i]][fB])
            res.buckets_decoded += 1
        res.log.append((bits, int(pending.size), len(overflow)))
        if overflow and (bits + 1 > max_bits or res.rounds >= max_rounds):
            raise RuntimeError(f"bucket split limit reached after {res.rounds} rounds")
        res.splits += len(overflow)
        ov = np.asarray(overflow, np.int64)
        pending = np.sort(np.concatenate([2 * ov, 2 * ov + 1])) if ov.size else np.empty(0, np.int64)
        bits += 1
    res.a_minus_b = np.sort(np.concatenate(a_out)) if a_out else np.empty(0, np.uint64)
    res.b_minus_a = np.sort(np.concatenate(b_out)) if b_out else np.empty(0, np.uint64)
    return res
