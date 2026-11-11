"""The Eisenstein KP-7 hash S(s) = (F + lam)(R + lam) over F_q, q = 2^61 - 1: ONE field element per k-mer.

Encoding A = 1, C = w, G = -w, T = -1, with w a primitive cube root of unity in F_q (q = 1 mod 3). Bit-identical to the
C library hash/kp7e.c (kp7e_direct + kp7e_sym): w, x and lam come from kp7e_init through ctypes, and `self-test`
compares this port with the C at every k-mer of a random sample. Collision bound: distinct double-stranded classes
collide with probability <= 2(k-1)/q (2.6e-17 at k = 31; results/kp7_eisenstein_2026-11-10.md), against
2(k-1)/p^2 for the Gaussian CM31 hash with x complex (recon/kp7hash.py, 1.3e-17).

k-mers are the canonical 2-bit codes of prototype/kmers.py (A0 C1 G2 T3, s_0 most significant), the same code order as
kp7e (complement 3 - c). F = sum_j val(s_j) x^j by Horner from s_{k-1}; R = F of the reverse-complement code.
"""

from __future__ import annotations

import ctypes
import pathlib
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from prototype.kmers import K  # noqa: E402
from recon import fq61  # noqa: E402
from recon.kp7hash import decode, rc_codes  # noqa: E402

HASH = pathlib.Path(__file__).resolve().parents[1] / "hash"
LIB = pathlib.Path("/home/rcb112/project/claude/bench/data/kmer_poly_work/libkp7e.so")
# Deterministic seeds for tests and benchmarks. params() without seeds draws x and lam from OS randomness, each
# independently (hash/kp7seed.h; round-1 adjudication L1 F1).
SEED_X, SEED_LAM = 0x5EED, 0x1A3D


class Params(ctypes.Structure):
    _fields_ = [("k", ctypes.c_uint), ("w", ctypes.c_uint64), ("val", ctypes.c_uint64 * 4), ("x", ctypes.c_uint64),
                ("xinv", ctypes.c_uint64), ("lam", ctypes.c_uint64), ("xk1", ctypes.c_uint64), ("xk", ctypes.c_uint64),
                ("tf", ctypes.c_uint64 * 16), ("tr", ctypes.c_uint64 * 16)]


def lib() -> ctypes.CDLL:
    src = HASH / "kp7e.c"
    newest = max(src.stat().st_mtime, (HASH / "kp7seed.h").stat().st_mtime)
    if not LIB.exists() or LIB.stat().st_mtime < newest:
        subprocess.run(["gcc", "-O2", "-shared", "-fPIC", "-march=native", "-o", str(LIB), str(src)], check=True)
    so = ctypes.CDLL(str(LIB))
    so.kp7e_init_seeds.argtypes = [ctypes.POINTER(Params), ctypes.c_uint, ctypes.POINTER(ctypes.c_uint64),
                                   ctypes.POINTER(ctypes.c_uint64)]
    so.kp7e_direct.argtypes = [ctypes.POINTER(Params), ctypes.c_char_p, ctypes.POINTER(ctypes.c_uint64),
                               ctypes.POINTER(ctypes.c_uint64)]
    so.kp7e_direct.restype = ctypes.c_int
    so.kp7e_sym.argtypes = [ctypes.POINTER(Params), ctypes.c_uint64, ctypes.c_uint64]
    so.kp7e_sym.restype = ctypes.c_uint64
    return so


def params(k: int = K, seed_x: int | None = None, seed_lam: int | None = None) -> Params:
    """x and lam from two independent sources: a given seed is a deterministic stream for that parameter, None draws
    it from OS randomness (getrandom). Both None is the default for use; tests and benchmarks pass both."""
    p = Params()
    sx = None if seed_x is None else ctypes.byref(ctypes.c_uint64(seed_x))
    sl = None if seed_lam is None else ctypes.byref(ctypes.c_uint64(seed_lam))
    lib().kp7e_init_seeds(ctypes.byref(p), k, sx, sl)
    return p


def poly(codes: np.ndarray, p: Params, k: int) -> np.ndarray:
    val = np.array([p.val[i] for i in range(4)], np.uint64)
    x = np.uint64(p.x)
    h = np.zeros(codes.size, np.uint64)
    for t in range(k):  # t = 0 is s_{k-1}, the least significant digit
        d = ((codes >> np.uint64(2 * t)) & np.uint64(3)).astype(np.intp)
        h = fq61.add(fq61.mul(h, x), val[d])
    return h


def sym_hash(codes: np.ndarray, p: Params, k: int = K, chunk: int = 4_000_000) -> np.ndarray:
    lam = np.uint64(p.lam)
    out = np.empty(codes.size, np.uint64)
    for s in range(0, codes.size, chunk):
        c = codes[s:s + chunk]
        out[s:s + chunk] = fq61.mul(fq61.add(poly(c, p, k), lam), fq61.add(poly(rc_codes(c, k), p, k), lam))
    return out


def self_test(n: int = 20000, seed_x: int = SEED_X, seed_lam: int = SEED_LAM) -> int:
    so, p = lib(), params(K, seed_x, seed_lam)
    rng = np.random.default_rng(7)
    codes = rng.integers(0, 1 << (2 * K), n, dtype=np.uint64)
    mine = sym_hash(codes, p)
    bad = 0
    for i in range(n):
        F, R = ctypes.c_uint64(), ctypes.c_uint64()
        assert so.kp7e_direct(ctypes.byref(p), decode(int(codes[i])), ctypes.byref(F), ctypes.byref(R)) == 1
        if so.kp7e_sym(ctypes.byref(p), F.value, R.value) != int(mine[i]):
            bad += 1
    sym_ok = np.array_equal(mine, sym_hash(rc_codes(codes), p))
    print(f"kp7e_hash self-test: {n - bad}/{n} equal to C kp7e_sym; strand-symmetric: {sym_ok}")
    return 0 if bad == 0 and sym_ok else 1


if __name__ == "__main__":
    sys.exit(self_test())
