"""The KP-7 strand-symmetric hash S(s) = (F + lam)(R + lam) over CM31 = F_{p^2}, p = 2^31 - 1, vectorized in numpy.

It is bit-identical to the C library hash/kp7m31.c (kp7m31_direct + kp7m31_sym): x and lam come from kp7m31_init
through ctypes, and `self-test` compares this port with the C at every k-mer of a random sample.

Element arithmetic is prototype/field.py (int64 pairs; every product of two reduced values is < 2^62, so no 128-bit
arithmetic is needed in this field). k-mers are the 2-bit canonical codes of prototype/kmers.py (A0 C1 G2 T3, s_0 in
the most significant digit).
    F = sum_j val(s_j) x^j       (prototype.kmers.poly_hash of the code)
    R = rc(p)(x) = F of rc(s)    (poly_hash of the reverse-complement code)
Because S is symmetric, hashing the canonical orientation or the other gives the same S.
"""

from __future__ import annotations

import ctypes
import pathlib
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from prototype import field as fq  # noqa: E402
from prototype.kmers import K, poly_hash  # noqa: E402

HASH = pathlib.Path(__file__).resolve().parents[1] / "hash"
LIB = pathlib.Path("/home/rcb112/project/claude/bench/data/kmer_poly_work/libkp7m31.so")
# Deterministic seeds for tests and benchmarks. params() without seeds draws x and lam from OS randomness, each
# independently (hash/kp7seed.h; round-1 adjudication L1 F1).
SEED_X, SEED_LAM = 0x5EED, 0x1A3D


class M31C(ctypes.Structure):
    _fields_ = [("re", ctypes.c_uint32), ("im", ctypes.c_uint32)]


class Params(ctypes.Structure):
    _fields_ = [("k", ctypes.c_uint), ("x", M31C), ("xinv", M31C), ("lam", M31C), ("xk1", M31C), ("xk", M31C),
                ("val", M31C * 4), ("tf", M31C * 16), ("tr", M31C * 16)]


def lib() -> ctypes.CDLL:
    """Build (if needed) and load the C library as a shared object. Scalar build: only its parameters and the scalar
    reference functions are used here."""
    src = HASH / "kp7m31.c"
    newest = max(src.stat().st_mtime, (HASH / "kp7seed.h").stat().st_mtime)
    if not LIB.exists() or LIB.stat().st_mtime < newest:
        subprocess.run(["gcc", "-O2", "-shared", "-fPIC", "-march=native", "-o", str(LIB), str(src)], check=True)
    so = ctypes.CDLL(str(LIB))
    so.kp7m31_init_seeds.argtypes = [ctypes.POINTER(Params), ctypes.c_uint, ctypes.POINTER(ctypes.c_uint64),
                                     ctypes.POINTER(ctypes.c_uint64)]
    so.kp7m31_direct.argtypes = [ctypes.POINTER(Params), ctypes.c_char_p, ctypes.POINTER(M31C), ctypes.POINTER(M31C)]
    so.kp7m31_direct.restype = ctypes.c_int
    so.kp7m31_sym.argtypes = [ctypes.POINTER(Params), M31C, M31C]
    so.kp7m31_sym.restype = ctypes.c_uint64
    return so


def params(k: int = K, seed_x: int | None = None, seed_lam: int | None = None) -> Params:
    """x and lam from two independent sources: a given seed is a deterministic stream for that parameter, None draws
    it from OS randomness (getrandom). Both None is the default for use; tests and benchmarks pass both."""
    p = Params()
    sx = None if seed_x is None else ctypes.byref(ctypes.c_uint64(seed_x))
    sl = None if seed_lam is None else ctypes.byref(ctypes.c_uint64(seed_lam))
    lib().kp7m31_init_seeds(ctypes.byref(p), k, sx, sl)
    return p


def rc_codes(codes: np.ndarray, k: int = K) -> np.ndarray:
    """Reverse-complement 2-bit codes (A0 C1 G2 T3: complement is 3 - c)."""
    out = np.zeros_like(codes)
    c = codes.copy()
    for _ in range(k):
        out = (out << np.uint64(2)) | (np.uint64(3) - (c & np.uint64(3)))
        c >>= np.uint64(2)
    return out


def sym_hash(codes: np.ndarray, p: Params, k: int = K, chunk: int = 4_000_000):
    """S(s) as an element (re, im) of int64 arrays in [0, p), for every k-mer code."""
    x = (int(p.x.re), int(p.x.im))
    lam = fq.F(int(p.lam.re), int(p.lam.im))
    re = np.empty(codes.size, np.int64)
    im = np.empty(codes.size, np.int64)
    for s in range(0, codes.size, chunk):
        c = codes[s:s + chunk]
        F = poly_hash(c, x, k)
        R = poly_hash(rc_codes(c, k), x, k)
        re[s:s + chunk], im[s:s + chunk] = fq.mul(fq.add(F, lam), fq.add(R, lam))
    return re, im


def pack(h) -> np.ndarray:
    """The C library's 62-bit packing: re | im << 31."""
    return h[0].astype(np.uint64) | (h[1].astype(np.uint64) << np.uint64(31))


def decode(code: int, k: int = K) -> bytes:
    return bytes(b"ACGT"[(code >> (2 * (k - 1 - t))) & 3] for t in range(k))


def self_test(n: int = 20000, seed_x: int = SEED_X, seed_lam: int = SEED_LAM) -> int:
    so, p = lib(), params(K, seed_x, seed_lam)
    rng = np.random.default_rng(7)
    codes = rng.integers(0, 1 << (2 * K), n, dtype=np.uint64)
    mine = pack(sym_hash(codes, p))
    bad = 0
    for i in range(n):
        F, R = M31C(), M31C()
        assert so.kp7m31_direct(ctypes.byref(p), decode(int(codes[i])), ctypes.byref(F), ctypes.byref(R)) == 1
        if so.kp7m31_sym(ctypes.byref(p), F, R) != int(mine[i]):
            bad += 1
    # strand symmetry
    sym_ok = np.array_equal(mine, pack(sym_hash(rc_codes(codes), p)))
    print(f"kp7hash self-test: {n - bad}/{n} equal to C kp7m31_sym; strand-symmetric: {sym_ok}")
    return 0 if bad == 0 and sym_ok else 1


if __name__ == "__main__":
    sys.exit(self_test())
