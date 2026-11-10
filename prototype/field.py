"""Vectorized arithmetic in F_{q^2} = Z[i]/(q), q = 2^31 - 1.

q is a prime with q = 3 (mod 4), so -1 has no square root mod q and Z[i]/(q) is a field of q^2 elements. The k-mer
encoding A = 1, C = i, G = -i, T = -1 (IDEAS.md section 0) lives in it exactly.

An element is a pair (re, im) of int64 numpy arrays (or scalars) with entries in [0, q). Every product of two reduced
values is below 2^62, and every sum of two such products below 2^63, so int64 never overflows. The production choice
named in IDEAS.md KP-6, q = 2^61 - 1, needs 128-bit products and is not used here.
"""

from __future__ import annotations

import numpy as np

Q = (1 << 31) - 1
assert Q % 4 == 3


def F(re, im=0):
    """Make an element (pair of int64 arrays) from integers."""
    return (np.asarray(re, dtype=np.int64) % Q, np.asarray(im, dtype=np.int64) % Q)


def add(a, b):
    return ((a[0] + b[0]) % Q, (a[1] + b[1]) % Q)


def sub(a, b):
    return ((a[0] - b[0]) % Q, (a[1] - b[1]) % Q)


def neg(a):
    return ((-a[0]) % Q, (-a[1]) % Q)


def mul(a, b):
    """(ar + ai i)(br + bi i) = (ar br - ai bi) + (ar bi + ai br) i."""
    ar, ai = a
    br, bi = b
    re = (ar * br + ((Q - ai) % Q) * bi) % Q
    im = (ar * bi + ai * br) % Q
    return (re, im)


def is_zero(a):
    return (a[0] == 0) & (a[1] == 0)


def eq(a, b):
    return (a[0] == b[0]) & (a[1] == b[1])


def _powmod(base, e: int):
    """base^e mod Q, elementwise on an int64 array."""
    base = np.asarray(base, dtype=np.int64) % Q
    result = np.ones_like(base)
    while e:
        if e & 1:
            result = (result * base) % Q
        base = (base * base) % Q
        e >>= 1
    return result


def inv(a):
    """1/a = conj(a) / (ar^2 + ai^2). The norm is nonzero for a != 0 because q = 3 (mod 4). inv(0) returns 0."""
    ar, ai = a
    norm = (ar * ar + ai * ai) % Q
    ninv = _powmod(norm, Q - 2)
    return ((ar * ninv) % Q, (((Q - ai) % Q) * ninv) % Q)


def div(a, b):
    return mul(a, inv(b))


def random_elements(rng: np.random.Generator, n: int):
    return (rng.integers(0, Q, n, dtype=np.int64), rng.integers(0, Q, n, dtype=np.int64))


def take(a, idx):
    return (a[0][idx], a[1][idx])


def prod_rows(a):
    """Product along the last axis of a 2-D element array, by log-depth pairwise halving."""
    re, im = a
    while re.shape[-1] > 1:
        n = re.shape[-1]
        if n % 2:
            re = np.concatenate([re, np.ones(re.shape[:-1] + (1,), np.int64)], axis=-1)
            im = np.concatenate([im, np.zeros(im.shape[:-1] + (1,), np.int64)], axis=-1)
        re, im = mul((re[..., 0::2], im[..., 0::2]), (re[..., 1::2], im[..., 1::2]))
    return (re[..., 0], im[..., 0])


def poly_eval(coeffs, x):
    """Evaluate sum coeffs[j] x^j by Horner. coeffs is a list of scalar elements (low degree first); x an array."""
    acc = (np.zeros_like(x[0]), np.zeros_like(x[1]))
    for c in reversed(coeffs):
        acc = add(mul(acc, x), (np.full_like(x[0], c[0]), np.full_like(x[1], c[1])))
    return acc
