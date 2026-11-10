"""Canonical k-mers as uint64 codes, their F_{q^2} polynomial hashes, and an independent bucket hash."""

from __future__ import annotations

import numpy as np

from . import field as fq

K = 31
# A=0, C=1, G=2, T=3, anything else (N) = 4. Case-blind: soft-masking is ignored.
_LUT = np.full(256, 4, dtype=np.uint8)
for _ch, _v in zip(b"ACGTacgt", (0, 1, 2, 3, 0, 1, 2, 3)):
    _LUT[_ch] = _v
# Base value in F_{q^2}: A = 1, C = i, G = -i, T = -1 (IDEAS.md section 0).
_VAL_RE = np.array([1, 0, 0, fq.Q - 1], dtype=np.int64)
_VAL_IM = np.array([0, 1, fq.Q - 1, 0], dtype=np.int64)


def read_fasta(path) -> bytes:
    """All sequence bytes of a single-record FASTA."""
    with open(path, "rb") as fh:
        lines = fh.read().split(b"\n")
    return b"".join(line.strip() for line in lines if line and not line.startswith(b">"))


def encode(seq: bytes) -> np.ndarray:
    return _LUT[np.frombuffer(seq, dtype=np.uint8)]


def canonical_kmers(codes: np.ndarray, k: int = K) -> np.ndarray:
    """Canonical k-mer codes (min of forward and reverse-complement 2-bit codes) of every N-free window, unsorted.

    The forward code puts s_0 in the most significant digit. The reverse complement of s_0..s_{k-1} has (3 - s_{k-1})
    first, so its code is sum (3 - s_j) 4^j.
    """
    n = codes.size - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.uint64)
    bad = np.concatenate([[0], np.cumsum(codes == 4)])
    ok = (bad[k:] - bad[:-k]) == 0
    c = np.where(codes == 4, 0, codes).astype(np.uint64)
    fwd = np.zeros(n, dtype=np.uint64)
    rc = np.zeros(n, dtype=np.uint64)
    for t in range(k):
        window = c[t:t + n]
        fwd |= window << np.uint64(2 * (k - 1 - t))
        rc |= (np.uint64(3) - window) << np.uint64(2 * t)
    return np.minimum(fwd, rc)[ok]


def kmer_set(codes: np.ndarray, k: int = K) -> np.ndarray:
    return np.unique(canonical_kmers(codes, k))


def poly_hash(kmer_codes: np.ndarray, x0, k: int = K):
    """h(s) = p_s(x0) = sum_j val(s_j) x0^j over F_{q^2}, by Horner from s_{k-1} down to s_0."""
    x = (np.full(kmer_codes.size, x0[0], np.int64), np.full(kmer_codes.size, x0[1], np.int64))
    h = (np.zeros(kmer_codes.size, np.int64), np.zeros(kmer_codes.size, np.int64))
    for t in range(k):  # t = 0 is s_{k-1}, the least significant digit
        d = ((kmer_codes >> np.uint64(2 * t)) & np.uint64(3)).astype(np.intp)
        h = fq.add(fq.mul(h, x), (_VAL_RE[d], _VAL_IM[d]))
    return h


def mix64(x: np.ndarray) -> np.ndarray:
    """splitmix64 finalizer: a bucket hash independent of the polynomial hash."""
    with np.errstate(over="ignore"):
        z = x + np.uint64(0x9E3779B97F4A7C15)
        z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return z ^ (z >> np.uint64(31))
