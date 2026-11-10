"""KP-10 end to end on a real genome pair: power-sum sketches of KP-7 strand-symmetric hashes, partitioned, decoded by
Newton's identities, verified element for element against exact set differences.

    python -m recon.kp10 run <genomeA> <genomeB> --bits B --m M --c C [--parts N] [--field gaussian|eisenstein]
    python -m recon.kp10 merge-test <genome> --bits B --m M [--field ...]  # sketch(A) == sketch(A1) + sketch(A2)

Two fields:
  gaussian    KP-7 over CM31 = F_{p^2}, p = 2^31 - 1, x complex (recon/kp7hash.py); prototype/recon.py decodes.
              One entry = 62 bits. Collision bound about 2k/p^2 = 1.4e-17 at k = 31.
  eisenstein  KP-7 over F_q, q = 2^61 - 1, A = 1, C = w, G = -w, T = -1 (recon/kp7e_hash.py); recon/recon61.py
              decodes. One entry = 61 bits. Collision bound about 2k/q = 2.7e-17 at k = 31.
Both entries fit 8 bytes. (A Gaussian sketch over F_{q^2} with q = 2^61 - 1 would need 16.)

Inputs are the partitioned canonical k-mer sets of recon/kmerset.py. Each element is the KP-7 hash S(s) over CM31
(recon/kp7hash.py, bit-identical to C kp7m31_sym). The bucket key is mix64(code), the same independent hash whose top
6 bits define kmerset's 64 partitions, so the buckets at resolution `bits` nest inside the partitions and every
partition can be reconciled on its own.

The protocol is prototype/recon.reconcile(method="power"): per bucket, power sums S_1..S_{m+c} and the bucket size.
Newton's identities give the series of chi_{A\\B}^rev / chi_{B\\A}^rev; a Pade step recovers both polynomials, and each
side flags its own elements that are roots. The c extra sums check the decode. An overflowing bucket (difference > m,
or a failed check) is split on the next key bit and sent again, at the cost of one more round.

Sketch size: every field element (a power sum or a bucket size) is counted at 8 bytes. A CM31 element is 62 bits.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from prototype import field as fq  # noqa: E402
from prototype.kmers import mix64  # noqa: E402
from prototype.recon import Side, power_sketch, reconcile  # noqa: E402
from recon import fq61, kmerset, recon61  # noqa: E402
from recon import kp7e_hash, kp7hash  # noqa: E402

RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"
PART_BITS = 6  # kmerset's 64 partitions


def side(codes: np.ndarray, p, fieldname: str = "gaussian"):
    if fieldname == "eisenstein":
        return recon61.SideE.build(codes, kp7e_hash.sym_hash(codes, p), mix64(codes))
    return Side.build(codes, kp7hash.sym_hash(codes, p), mix64(codes))


def field_params(fieldname: str, k: int):
    return kp7e_hash.params(k) if fieldname == "eisenstein" else kp7hash.params(k)


def part_buckets(part: int, bits: int) -> np.ndarray:
    span = 1 << (bits - PART_BITS)
    return np.arange(part * span, (part + 1) * span, dtype=np.int64)


def run(ga: str, gb: str, k: int, bits: int, m: int, c: int, nparts: int, fieldname: str = "gaussian") -> dict:
    assert bits >= PART_BITS
    p = field_params(fieldname, k)
    tot = {"field": fieldname, "A": ga, "B": gb, "k": k, "bits": bits, "m": m, "c": c, "parts": nparts,
           "D_recovered": 0, "exact": True,
           "sketch_elements": 0, "rounds_max": 0, "splits": 0, "seconds_hash": 0.0, "seconds_reconcile": 0.0,
           "elements_A": 0, "elements_B": 0}
    for part in range(nparts):
        a, b = kmerset.load_part(ga, k, part), kmerset.load_part(gb, k, part)
        t0 = time.time()
        A, B = side(a, p, fieldname), side(b, p, fieldname)
        t1 = time.time()
        if fieldname == "eisenstein":
            res = recon61.reconcile(A, B, bits, m, c, buckets=part_buckets(part, bits))
        else:
            res = reconcile(A, B, bits, m, c, None, method="power", buckets=part_buckets(part, bits))
        t2 = time.time()
        want_amb = np.setdiff1d(a, b, assume_unique=True)
        want_bma = np.setdiff1d(b, a, assume_unique=True)
        ok = np.array_equal(res.a_minus_b, want_amb) and np.array_equal(res.b_minus_a, want_bma)
        tot["exact"] &= bool(ok)
        tot["D_recovered"] += int(res.a_minus_b.size + res.b_minus_a.size)
        tot["sketch_elements"] += res.sketch_elements
        tot["rounds_max"] = max(tot["rounds_max"], res.rounds)
        tot["splits"] += res.splits
        tot["seconds_hash"] += t1 - t0
        tot["seconds_reconcile"] += t2 - t1
        tot["elements_A"] += int(a.size)
        tot["elements_B"] += int(b.size)
        print(f"part {part}: |A|={a.size} |B|={b.size} D={want_amb.size + want_bma.size} exact={ok} "
              f"rounds={res.rounds} splits={res.splits} hash {t1 - t0:.1f}s reconcile {t2 - t1:.1f}s", flush=True)
    tot["sketch_bytes"] = tot["sketch_elements"] * 8
    tot["bytes_per_difference"] = tot["sketch_bytes"] / max(tot["D_recovered"], 1)
    tot["seconds_hash"] = round(tot["seconds_hash"], 1)
    tot["seconds_reconcile"] = round(tot["seconds_reconcile"], 1)
    return tot


def merge_test(g: str, k: int, bits: int, m: int, part: int = 0, fieldname: str = "gaussian") -> bool:
    """Linearity: the power-sum sketch of a set equals the sum of the sketches of a split of it."""
    p = field_params(fieldname, k)
    a = kmerset.load_part(g, k, part)
    rng = np.random.default_rng(3)
    mask = rng.random(a.size) < 0.5
    bk = part_buckets(part, bits)
    if fieldname == "eisenstein":
        full, nf = recon61.power_sketch(side(a, p, fieldname), bits, bk, m)
        h1, n1 = recon61.power_sketch(side(a[mask], p, fieldname), bits, bk, m)
        h2, n2 = recon61.power_sketch(side(a[~mask], p, fieldname), bits, bk, m)
        ok = np.array_equal(full, fq61.add(h1, h2)) and np.array_equal(nf, n1 + n2)
    else:
        full, nf = power_sketch(side(a, p), bits, bk, m)
        h1, n1 = power_sketch(side(a[mask], p), bits, bk, m)
        h2, n2 = power_sketch(side(a[~mask], p), bits, bk, m)
        summed = ((h1[0] + h2[0]) % fq.Q, (h1[1] + h2[1]) % fq.Q)
        ok = np.array_equal(full[0], summed[0]) and np.array_equal(full[1], summed[1]) and np.array_equal(nf, n1 + n2)
    print(f"merge-test [{fieldname}] {g} part {part}: {a.size} elements in two halves of {mask.sum()} and {(~mask).sum()}; "
          f"sketch(A) == sketch(A1) + sketch(A2): {ok}")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "merge-test"])
    ap.add_argument("genomes", nargs="+")
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--bits", type=int, required=True)
    ap.add_argument("--m", type=int, default=16)
    ap.add_argument("--c", type=int, default=2)
    ap.add_argument("--parts", type=int, default=64)
    ap.add_argument("--field", choices=["gaussian", "eisenstein"], default="gaussian")
    ap.add_argument("--out", type=pathlib.Path)
    a = ap.parse_args()
    if a.cmd == "merge-test":
        return 0 if merge_test(a.genomes[0], a.k, a.bits, a.m, fieldname=a.field) else 1
    res = run(a.genomes[0], a.genomes[1], a.k, a.bits, a.m, a.c, a.parts, a.field)
    print(json.dumps(res))
    if a.out:
        a.out.write_text(json.dumps(res, indent=1) + "\n")
    return 0 if res["exact"] else 1


if __name__ == "__main__":
    sys.exit(main())
