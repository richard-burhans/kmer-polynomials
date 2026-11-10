"""Vectorized arithmetic in F_q, q = 2^61 - 1, on uint64 numpy arrays, with no 128-bit types.

mul splits each operand into a 30-bit high and a 31-bit low limb, a = a1*2^31 + a0, so that
    a*b = a1 b1 2^62 + (a1 b0 + a0 b1) 2^31 + a0 b0,
and uses 2^61 = 1 (mod q): 2^62 = 2, and mid*2^31 = m1*2^61 + m0*2^31 = m1 + m0*2^31 with mid = m1*2^30 + m0.
Every partial product is < 2^62 and the total is < 2^63, so uint64 never overflows. Two Mersenne folds and one
conditional subtraction give the canonical value in [0, q).

summod adds many values < q without overflow by summing the high and low limbs separately.
Scalar decoding (Newton, Pade, Gauss-Jordan) uses Python integers, which are exact.
"""

from __future__ import annotations

import numpy as np

Q = (1 << 61) - 1
_Q = np.uint64(Q)
_M31 = np.uint64((1 << 31) - 1)
_M30 = np.uint64((1 << 30) - 1)
_S31, _S30, _S61, _ONE = np.uint64(31), np.uint64(30), np.uint64(61), np.uint64(1)


def red(s):
    """[0, 2^64) -> [0, q)."""
    r = (s & _Q) + (s >> _S61)
    r = (r & _Q) + (r >> _S61)
    return r - _Q * (r >= _Q).astype(np.uint64)


def mul(a, b):
    a = np.asarray(a, np.uint64)
    b = np.asarray(b, np.uint64)
    a1, a0 = a >> _S31, a & _M31
    b1, b0 = b >> _S31, b & _M31
    mid = a1 * b0 + a0 * b1
    s = ((a1 * b1) << _ONE) + (mid >> _S30) + ((mid & _M30) << _S31) + a0 * b0
    return red(s)


def add(a, b):
    s = np.asarray(a, np.uint64) + np.asarray(b, np.uint64)
    return s - _Q * (s >= _Q).astype(np.uint64)


def sub(a, b):
    a = np.asarray(a, np.uint64)
    b = np.asarray(b, np.uint64)
    return a + _Q * (a < b).astype(np.uint64) - b


def summod(x, axis=None):
    """sum(x) mod q for values < q, any count up to 2^33 per reduction, without overflow."""
    x = np.asarray(x, np.uint64)
    hi = red((x >> _S31).sum(axis=axis, dtype=np.uint64))
    lo = red((x & _M31).sum(axis=axis, dtype=np.uint64))
    return add(mul(hi, np.uint64(1 << 31)), lo)


def inv(a: int) -> int:
    return pow(int(a), Q - 2, Q)


def self_test(n: int = 200000) -> bool:
    rng = np.random.default_rng(11)
    a = rng.integers(0, Q, n, dtype=np.uint64)
    b = rng.integers(0, Q, n, dtype=np.uint64)
    edge = np.array([0, 1, 2, Q - 1, Q - 2, (1 << 60), (1 << 31) - 1, 1 << 31], np.uint64)
    a = np.concatenate([a, np.repeat(edge, edge.size)])
    b = np.concatenate([b, np.tile(edge, edge.size)])
    pm = mul(a, b)
    ok = all(int(pm[i]) == (int(a[i]) * int(b[i])) % Q for i in range(0, a.size, 97)) and \
        all(int(pm[i]) == (int(a[i]) * int(b[i])) % Q for i in range(a.size - edge.size ** 2, a.size))
    ok &= int(summod(a)) == sum(int(v) for v in a) % Q
    ok &= all(int(add(a, b)[i]) == (int(a[i]) + int(b[i])) % Q for i in range(0, a.size, 997))
    ok &= all(int(sub(a, b)[i]) == (int(a[i]) - int(b[i])) % Q for i in range(0, a.size, 997))
    return bool(ok)


if __name__ == "__main__":
    print("fq61 self-test:", "ok" if self_test() else "FAILED")
