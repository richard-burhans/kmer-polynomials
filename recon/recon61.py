"""KP-10 power-sum reconciliation over F_q, q = 2^61 - 1, for the Eisenstein KP-7 hash: one 8-byte element per entry.

The same protocol as prototype/recon.py (method "power"), specialised to a prime field:
  sketch   per bucket: power sums S_1..S_{m+c} of its elements and its size (vectorized, recon/fq61)
  decode   S = S^A - S^B; Newton's identities give r = exp(-sum S_j t^j / j) = P_rev / Q_rev with
           P_rev = prod_{a in A\\B} (1 - a t), Q_rev = prod_{b in B\\A} (1 - b t). A Pade step of degrees
           (dA, dB) = ((mbar + delta)/2, (mbar - delta)/2) recovers both; free variables are set to 0, as in the
           prototype. Each side flags its own elements a with P(a) = 0, where P is P_rev with its coefficients
           reversed (P(a) = a^dA P_rev(1/a), the same zeros for a != 0, and no modular inverse per element).
           The c extra sums must be reproduced by the flagged elements, and the flagged counts must differ by
           delta = |A_bucket| - |B_bucket|; otherwise the bucket overflows.
  split    an overflowing bucket is split on the next key bit and sent again (one more round).
Per-bucket decoding uses exact Python integers.
"""

from __future__ import annotations

from dataclasses import dataclass, field as dc_field

import numpy as np

from recon import fq61

Q = fq61.Q


@dataclass
class SideE:
    h: np.ndarray      # uint64 hashes in [0, q)
    key: np.ndarray    # uint64 bucket keys
    codes: np.ndarray  # uint64 k-mer codes

    @classmethod
    def build(cls, codes, h, key):
        order = np.argsort(key, kind="stable")
        return cls(h[order], key[order], codes[order])

    def bucket_slices(self, bits: int, buckets: np.ndarray):
        b = (self.key >> np.uint64(64 - bits)).astype(np.int64)
        return np.searchsorted(b, buckets, side="left"), np.searchsorted(b, buckets, side="right")


def power_sketch(side: SideE, bits: int, buckets: np.ndarray, m: int, chunk: int = 256):
    """(nb, m) uint64 power sums S_1..S_m per bucket, plus bucket sizes. Padding is 0, which adds nothing."""
    lo, hi = side.bucket_slices(bits, buckets)
    out = np.empty((buckets.size, m), np.uint64)
    for s in range(0, buckets.size, chunk):
        l, h = lo[s:s + chunk], hi[s:s + chunk]
        L = max(int((h - l).max()) if l.size else 1, 1)
        idx = l[:, None] + np.arange(L)[None, :]
        valid = idx < h[:, None]
        H = np.where(valid, side.h[np.where(valid, idx, 0)], np.uint64(0))
        pw = H
        for j in range(m):
            out[s:s + chunk, j] = fq61.summod(pw, axis=1)
            pw = fq61.mul(pw, H)
    return out, hi - lo


def _solve(M: list[list[int]], rhs: list[int]) -> list[int] | None:
    """Gauss-Jordan mod q; free variables 0; None if inconsistent."""
    n, cols = len(M), len(M[0]) if M else 0
    A = [row[:] + [r] for row, r in zip(M, rhs)]
    piv, r = [], 0
    for c in range(cols):
        p = next((i for i in range(r, n) if A[i][c] % Q), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = pow(A[r][c], Q - 2, Q)
        A[r] = [v * inv % Q for v in A[r]]
        for i in range(n):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(vi - f * vr) % Q for vi, vr in zip(A[i], A[r])]
        piv.append(c)
        r += 1
        if r == n:
            break
    if any(A[i][cols] % Q for i in range(r, n)):
        return None
    x = [0] * cols
    for i, c in enumerate(piv):
        x[c] = A[i][cols]
    return x


def _flags(desc_coeffs: list[int], own: np.ndarray) -> np.ndarray:
    """Zeros of sum desc_coeffs[i] t^(deg-i) at the elements `own` (Horner, vectorized)."""
    acc = np.zeros(own.size, np.uint64)
    for c in desc_coeffs:
        acc = fq61.add(fq61.mul(acc, own), np.uint64(c))
    return acc == 0


def _psum(elems: np.ndarray, j: int) -> int:
    return sum(pow(int(e), j, Q) for e in elems) % Q


def decode(SA, SB, nA: int, nB: int, m: int, c: int, ownA: np.ndarray, ownB: np.ndarray):
    S = [(int(a) - int(b)) % Q for a, b in zip(SA, SB)]  # S[j-1] = S_j
    delta = int(nA) - int(nB)
    mbar = m if (m - delta) % 2 == 0 else m - 1
    if abs(delta) > mbar:
        return None
    dA, dB = (mbar + delta) // 2, (mbar - delta) // 2
    r = [1] + [0] * mbar
    for k in range(1, mbar + 1):
        acc = sum(S[j - 1] * r[k - j] for j in range(1, k + 1)) % Q
        r[k] = (-acc * pow(k, Q - 2, Q)) % Q
    if dB:
        # r_i = 0 for i < 0: a negative Python index would wrap to the END of r (the bug fixed 2026-11-11; it hit
        # every bucket with dB > dA + 1, i.e. more differences on B's side)
        M = [[r[n - j] if n - j >= 0 else 0 for j in range(1, dB + 1)] for n in range(dA + 1, dA + dB + 1)]
        sol = _solve(M, [(-r[n]) % Q for n in range(dA + 1, dA + dB + 1)])
        if sol is None:
            return None
        qr = [1] + sol
    else:
        qr = [1]
    pr = [sum(qr[j] * r[n - j] for j in range(0, min(n, dB) + 1)) % Q for n in range(dA + 1)]
    flagsA = _flags(pr, ownA)   # pr = P_rev low first = P high first
    flagsB = _flags(qr, ownB)
    fa, fb = ownA[flagsA], ownB[flagsB]
    for j in range(m + 1, m + c + 1):
        if (_psum(fa, j) - _psum(fb, j)) % Q != S[j - 1]:
            return None
    if int(flagsA.sum()) - int(flagsB.sum()) != delta:
        return None
    return flagsA, flagsB


@dataclass
class Result:
    a_minus_b: np.ndarray
    b_minus_a: np.ndarray
    rounds: int = 0
    splits: int = 0
    sketch_elements: int = 0
    log: list = dc_field(default_factory=list)
    unresolved: int = 0  # overflowing buckets left when stop_at_limit ended the protocol early


def reconcile(A: SideE, B: SideE, bits: int, m: int, c: int, buckets: np.ndarray | None = None,
              max_bits: int = 62, max_rounds: int = 12, stop_at_limit: bool = False) -> Result:
    """stop_at_limit: at the round or bit limit, return what was decoded and count the still-overflowing buckets in
    Result.unresolved instead of raising (used to measure a control that cannot finish)."""
    pending = np.arange(1 << bits, dtype=np.int64) if buckets is None else np.asarray(buckets, np.int64)
    a_out, b_out = [], []
    res = Result(np.empty(0, np.uint64), np.empty(0, np.uint64))
    while pending.size:
        res.rounds += 1
        skA, nA = power_sketch(A, bits, pending, m + c)
        skB, nB = power_sketch(B, bits, pending, m + c)
        res.sketch_elements += pending.size * (m + c + 1)
        loA, hiA = A.bucket_slices(bits, pending)
        loB, hiB = B.bucket_slices(bits, pending)
        overflow = []
        for i, bkt in enumerate(pending):
            if nA[i] == nB[i] and np.array_equal(skA[i], skB[i]):
                continue  # identical sketches: no difference (with probability 1 - O(m/q))
            out = decode(skA[i], skB[i], nA[i], nB[i], m, c, A.h[loA[i]:hiA[i]], B.h[loB[i]:hiB[i]])
            if out is None:
                overflow.append(bkt)
                continue
            fA, fB = out
            a_out.append(A.codes[loA[i]:hiA[i]][fA])
            b_out.append(B.codes[loB[i]:hiB[i]][fB])
        res.log.append((bits, int(pending.size), len(overflow)))
        if overflow and (bits + 1 > max_bits or res.rounds >= max_rounds):
            if not stop_at_limit:
                raise RuntimeError(f"bucket split limit reached after {res.rounds} rounds")
            res.unresolved = len(overflow)
            break
        res.splits += len(overflow)
        ov = np.asarray(overflow, np.int64)
        pending = np.sort(np.concatenate([2 * ov, 2 * ov + 1])) if ov.size else np.empty(0, np.int64)
        bits += 1
    res.a_minus_b = np.sort(np.concatenate(a_out)) if a_out else np.empty(0, np.uint64)
    res.b_minus_a = np.sort(np.concatenate(b_out)) if b_out else np.empty(0, np.uint64)
    return res
