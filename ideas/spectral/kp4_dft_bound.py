"""KP-4: does a truncated-DFT (Parseval) lower bound prune k-mer pairs at mismatch thresholds?

With the complex encoding, ||p - q||^2 = sum |a_t - b_t|^2 (2 per ordinary mismatch, 4 per complementary one), and by
Parseval (numpy's unnormalised DFT) ||p - q||^2 = (1/k) sum_j |P_j - Q_j|^2. Keeping any m coefficients gives a lower
bound LB_m <= ||p - q||^2. A pair with Hamming distance <= d has ||p - q||^2 <= 4d, so LB_m > 4d proves Hamming > d
(safe pruning). This script measures, on real k-mers:
  - the fraction of pairs each m prunes, against the fraction that COULD be pruned (true Hamming > d);
  - how much of ||p - q||^2 the first m low-|frequency| coefficients capture, against m/k for white noise (the Voss
    note: k-mer-length windows sit in the white-noise regime);
  - wall time per pair: LB_m on precomputed coefficients vs exact Hamming on 2-bit codes (XOR + popcount).
Two candidate sets: random pairs, and near neighbours (distinct k-mers sharing a `seed`-base prefix).

    kp4_dft_bound.py FASTA [--k 31] [--pairs 200000] [--seed-len 12] [--out results.json]
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from common import read_fasta_one

M_LIST = list(range(1, 9))
D_LIST = [1, 2, 3, 5]


def low_freq_order(k: int) -> np.ndarray:
    """0, 1, k-1, 2, k-2, ...: coefficients by increasing |frequency| (a complex signal has no conjugate symmetry)."""
    order = [0]
    for f in range(1, k // 2 + 1):
        order.append(f)
        if k - f != f:
            order.append(k - f)
    return np.array(order[:k])


def kmer_windows(codes: np.ndarray, k: int, starts: np.ndarray) -> np.ndarray:
    idx = starts[:, None] + np.arange(k)[None, :]
    return codes[idx]


def valid_starts(codes: np.ndarray, k: int) -> np.ndarray:
    bad = (codes == 4).astype(np.int32)
    c = np.concatenate([[0], np.cumsum(bad)])
    nbad = c[k:] - c[:-k]
    return np.nonzero(nbad == 0)[0]


def pack2(win: np.ndarray) -> np.ndarray:
    """2-bit pack (k <= 32) into uint64 for the popcount Hamming baseline."""
    w = win.astype(np.uint64)
    out = np.zeros(len(win), dtype=np.uint64)
    for t in range(win.shape[1]):
        out = (out << np.uint64(2)) | w[:, t]
    return out


def popcount64(x: np.ndarray) -> np.ndarray:
    return np.unpackbits(x.view(np.uint8).reshape(-1, 8), axis=1).sum(axis=1)


def hamming_packed(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    x = a ^ b
    lo = np.uint64(0x5555555555555555)
    y = (x | (x >> np.uint64(1))) & lo
    return popcount64(y)


def analyse(A: np.ndarray, Bw: np.ndarray, k: int, label: str) -> dict:
    vals = np.array([1, 1j, -1j, -1], dtype=np.complex128)
    a, b = vals[A], vals[Bw]
    ham = (A != Bw).sum(axis=1)
    l2 = (np.abs(a - b) ** 2).sum(axis=1)
    FA, FB = np.fft.fft(a, axis=1), np.fft.fft(b, axis=1)
    diff2 = np.abs(FA - FB) ** 2 / k  # per-coefficient share of ||p-q||^2; rows sum to l2
    assert np.allclose(diff2.sum(axis=1), l2), "Parseval check failed"
    res = {"label": label, "pairs": int(len(A)), "hamming_mean": float(ham.mean()),
           "hamming_quantiles": {q: float(np.quantile(ham, q)) for q in (0.01, 0.1, 0.5, 0.9)}, "rows": []}
    for name, order in (("first", np.arange(k)), ("lowfreq", low_freq_order(k))):
        cum = np.cumsum(diff2[:, order], axis=1)
        for m in M_LIST:
            lb = cum[:, m - 1]
            nz = l2 > 0
            share = float((lb[nz] / l2[nz]).mean())
            row = {"order": name, "m": m, "energy_share": share, "white_noise_share": m / k}
            for d in D_LIST:
                prunable = ham > d
                pruned = lb > 4 * d + 1e-9
                assert not np.any(pruned & ~prunable), "unsafe prune"  # the bound must never prune a true neighbour
                row[f"d{d}_pruned"] = float(pruned.mean())
                row[f"d{d}_prunable"] = float(prunable.mean())
                row[f"d{d}_recall_of_prunable"] = float(pruned.sum() / max(prunable.sum(), 1))
            res["rows"].append(row)
    return res


def timing(A: np.ndarray, Bw: np.ndarray, k: int, m: int = 8) -> dict:
    vals = np.array([1, 1j, -1j, -1], dtype=np.complex64)
    FA = np.fft.fft(vals[A], axis=1)[:, :m].astype(np.complex64)
    FB = np.fft.fft(vals[Bw], axis=1)[:, :m].astype(np.complex64)
    pa, pb = pack2(A), pack2(Bw)
    t0 = time.perf_counter()
    for _ in range(5):
        lb = (np.abs(FA - FB) ** 2).sum(axis=1) / k
    t1 = time.perf_counter()
    for _ in range(5):
        h = hamming_packed(pa, pb)
    t2 = time.perf_counter()
    assert np.array_equal(h, (A != Bw).sum(axis=1)), "packed Hamming check failed"
    n = len(A) * 5
    return {"m": m, "ns_per_pair_lb": (t1 - t0) / n * 1e9, "ns_per_pair_hamming": (t2 - t1) / n * 1e9,
            "bytes_per_kmer_lb": 8 * m, "bytes_per_kmer_packed": 8, "note": "numpy, vectorised; ratios, not absolute C costs"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta", type=pathlib.Path)
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--pairs", type=int, default=200000)
    ap.add_argument("--seed-len", type=int, default=12)
    ap.add_argument("--rng", type=int, default=1)
    ap.add_argument("--out", type=pathlib.Path)
    a = ap.parse_args()
    k, rng = a.k, np.random.default_rng(a.rng)
    name, codes, _ = read_fasta_one(a.fasta)
    vs = valid_starts(codes, k)
    out = {"fasta": str(a.fasta), "record": name, "k": k, "valid_kmers": int(len(vs))}
    # random pairs
    i, j = rng.choice(vs, a.pairs), rng.choice(vs, a.pairs)
    A, Bw = kmer_windows(codes, k, i), kmer_windows(codes, k, j)
    out["random"] = analyse(A, Bw, k, "random pairs")
    # near neighbours: distinct k-mers sharing a seed-length prefix, from a sample of 4M k-mers
    s = rng.choice(vs, min(len(vs), 4_000_000), replace=False)
    W = kmer_windows(codes, k, s)
    sl = a.seed_len
    key = np.zeros(len(W), dtype=np.uint64)
    for t in range(sl):
        key = (key << np.uint64(2)) | W[:, t].astype(np.uint64)
    order = np.argsort(key, kind="stable")
    key, W = key[order], W[order]
    same = key[1:] == key[:-1]
    distinct = np.any(W[1:] != W[:-1], axis=1)
    pick = np.nonzero(same & distinct)[0]
    if len(pick) > a.pairs:
        pick = rng.choice(pick, a.pairs, replace=False)
    out["seeded"] = analyse(W[pick], W[pick + 1], k, f"distinct k-mers sharing a {sl}-base prefix")
    out["timing_m8"] = timing(A[:100000], Bw[:100000], k, 8)
    js = json.dumps(out, indent=1)
    if a.out:
        a.out.write_text(js + "\n")
    print(js)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
