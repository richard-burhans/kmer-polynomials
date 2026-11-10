"""KP-9 experiments on KOMPa vs KOMPb chr7 (two haplotypes of one plant) and on a near-identical simulated pair.

    python -m prototype.run_kompa sets      # exact canonical 31-mer sets, |A|, |B|, D, hash collisions
    python -m prototype.run_kompa sample    # reconcile a random sample of buckets of the real haplotype pair
    python -m prototype.run_kompa sim N     # KOMPa chr7 vs itself with N simulated SNPs: full reconciliation

Work files go to bench/data/kmer_poly_work. Every recovered difference is checked element for element against the
exact numpy set difference.
"""

from __future__ import annotations

import json
import pathlib
import sys
import time

import numpy as np

from . import field as fq
from . import kmers, recon

W = pathlib.Path("/home/rcb112/project/claude/bench/data/kmer_poly_work")
M, C = 64, 4  # solve points and check points per bucket
SEED = 20261110


def load_set(name: str) -> np.ndarray:
    npy = W / f"{name}.k31.npy"
    if npy.exists():
        return np.load(npy)
    s = kmers.kmer_set(kmers.encode(kmers.read_fasta(W / f"{name}.fa")))
    np.save(npy, s)
    return s


def x0_and_points():
    rng = np.random.default_rng(SEED)
    x0 = fq.random_elements(rng, 1)
    return (x0[0][0], x0[1][0]), fq.random_elements(rng, M + C)


def side(codes, x0):
    return recon.Side.build(codes, kmers.poly_hash(codes, x0), kmers.mix64(codes))


def bits_for(D: int) -> int:
    """Smallest bucket count 2^bits giving an expected per-bucket difference of at most M/2."""
    bits = 0
    while D / (1 << bits) > M / 2:
        bits += 1
    return bits


def cmd_sets() -> dict:
    t = time.time()
    A, B = load_set("KOMPa.chr7"), load_set("KOMPb.chr7")
    t_sets = time.time() - t
    amb = np.setdiff1d(A, B, assume_unique=True)
    bma = np.setdiff1d(B, A, assume_unique=True)
    inter = A.size - amb.size
    x0, _ = x0_and_points()
    U = np.union1d(A, B)
    h = kmers.poly_hash(U, x0)
    distinct_h = np.unique(h[0] * (1 << 31) + h[1]).size
    out = {"A": int(A.size), "B": int(B.size), "A_minus_B": int(amb.size), "B_minus_A": int(bma.size),
           "intersection": int(inter), "D": int(amb.size + bma.size), "union": int(U.size),
           "distinct_hashes_of_union": int(distinct_h), "seconds_to_build_sets": round(t_sets, 1)}
    out["D_over_A"] = round(out["D"] / out["A"], 4)
    out["bits"] = bits_for(out["D"])
    out["sketch_elements"] = (1 << out["bits"]) * (M + C + 1)
    print(json.dumps(out, indent=1))
    (W / "sets_summary.json").write_text(json.dumps(out, indent=1))
    return out


def verify(res, A, B, keysA_sel, keysB_sel):
    """Compare recovered differences with exact set differences restricted to the selected elements."""
    want_amb = np.sort(np.setdiff1d(A[keysA_sel], B, assume_unique=True))
    want_bma = np.sort(np.setdiff1d(B[keysB_sel], A, assume_unique=True))
    return (np.array_equal(res.a_minus_b, want_amb) and np.array_equal(res.b_minus_a, want_bma),
            int(want_amb.size), int(want_bma.size))


def cmd_sample(n_buckets: int = 64, method: str = "eval") -> dict:
    A, B = load_set("KOMPa.chr7"), load_set("KOMPb.chr7")
    summ = json.loads((W / "sets_summary.json").read_text())
    bits = summ["bits"]
    x0, z = x0_and_points()
    rng = np.random.default_rng(SEED + 1)
    buckets = np.sort(rng.choice(1 << bits, n_buckets, replace=False))
    # restrict both sides to the sampled buckets (the protocol only touches those)
    kA = (kmers.mix64(A) >> np.uint64(64 - bits)).astype(np.int64)
    kB = (kmers.mix64(B) >> np.uint64(64 - bits)).astype(np.int64)
    selA, selB = np.isin(kA, buckets), np.isin(kB, buckets)
    sa, sb = side(A[selA], x0), side(B[selB], x0)
    t = time.time()
    res = recon.reconcile(sa, sb, bits=bits, m=M, c=C, z=z, method=method, buckets=buckets)
    dt = time.time() - t
    ok, na, nb = verify(res, A, B, selA, selB)
    out = {"pair": "KOMPa.chr7 vs KOMPb.chr7", "method": method, "bits": bits, "sampled_buckets": n_buckets,
           "elements_A": int(selA.sum()), "elements_B": int(selB.sum()), "diff_A_minus_B": na, "diff_B_minus_A": nb,
           "exact": ok, "rounds": res.rounds, "splits": res.splits, "log": res.log, "seconds": round(dt, 1),
           "sketch_elements": res.sketch_elements}
    print(json.dumps(out, indent=1))
    (W / f"sample_{method}.json").write_text(json.dumps(out, indent=1))
    return out


def mutate(seq: bytes, n_snp: int, rng) -> tuple[bytes, np.ndarray]:
    arr = np.frombuffer(seq, dtype=np.uint8).copy()
    upper = arr & 0xDF
    ok = np.nonzero(np.isin(upper, np.frombuffer(b"ACGT", np.uint8)))[0]
    pos = np.sort(rng.choice(ok, n_snp, replace=False))
    bases = np.frombuffer(b"ACGT", np.uint8)
    for p in pos:
        cur = upper[p]
        alt = rng.choice(bases[bases != cur])
        arr[p] = alt | (arr[p] & 0x20)  # keep the case
    return arr.tobytes(), pos


def cmd_sim(n_snp: int, method: str = "eval") -> dict:
    rng = np.random.default_rng(SEED + n_snp)
    seq = kmers.read_fasta(W / "KOMPa.chr7.fa")
    A = load_set("KOMPa.chr7")
    mut, _ = mutate(seq, n_snp, rng)
    B = kmers.kmer_set(kmers.encode(mut))
    amb = np.setdiff1d(A, B, assume_unique=True)
    bma = np.setdiff1d(B, A, assume_unique=True)
    D = int(amb.size + bma.size)
    bits = bits_for(max(D, 1))
    x0, z = x0_and_points()
    t = time.time()
    sa, sb = side(A, x0), side(B, x0)
    t_hash = time.time() - t
    t = time.time()
    res = recon.reconcile(sa, sb, bits=bits, m=M, c=C, z=z, method=method)
    dt = time.time() - t
    ok = np.array_equal(res.a_minus_b, np.sort(amb)) and np.array_equal(res.b_minus_a, np.sort(bma))
    out = {"pair": f"KOMPa.chr7 vs KOMPa.chr7 + {n_snp} SNPs", "method": method, "A": int(A.size), "B": int(B.size),
           "A_minus_B": int(amb.size), "B_minus_A": int(bma.size), "D": D, "D_over_A": round(D / A.size, 6),
           "bits": bits, "buckets_initial": 1 << bits, "rounds": res.rounds, "splits": res.splits, "log": res.log,
           "sketch_elements": res.sketch_elements, "sketch_bytes": res.sketch_elements * 8,
           "set_bytes": int(A.size) * 8, "exact": ok, "seconds_hash_and_sort": round(t_hash, 1),
           "seconds_reconcile": round(dt, 1)}
    print(json.dumps(out, indent=1))
    (W / f"sim_{n_snp}_{method}.json").write_text(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    cmd = sys.argv[1]
    method = sys.argv[3] if len(sys.argv) > 3 else "eval"
    if cmd == "sets":
        cmd_sets()
    elif cmd == "sample":
        cmd_sample(int(sys.argv[2]) if len(sys.argv) > 2 else 64, method)
    elif cmd == "sim":
        cmd_sim(int(sys.argv[2]), method)
