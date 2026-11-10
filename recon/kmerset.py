"""Whole-genome canonical k-mer sets on disk, partitioned by an independent hash, and exact set differences.

    python -m recon.kmerset build <genome> [--k 31] [--parts 64]   # from hub-build/2bit/<genome>.2bit
    python -m recon.kmerset diff <genomeA> <genomeB>                # |A|, |B|, |A\\B|, |B\\A|, D; writes the differences

A genome of about 900 Mb has about 800 M distinct canonical 31-mers, 6.4 GB as uint64, which does not fit in memory
twice on this host (15 GB). So each sequence is streamed (twoBitToFa -seq), its canonical k-mers are split into
`parts` partitions by the top bits of mix64(code), and each partition is appended to its own file. A final pass
sorts and deduplicates each partition. Differences are then computed partition by partition.

Canonical = min(forward code, reverse-complement code), the same as prototype/kmers.canonical_kmers; soft-masking is
ignored (case-blind); windows containing N are skipped.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from prototype.kmers import canonical_kmers, encode, mix64  # noqa: E402

BENCH = pathlib.Path("/home/rcb112/project/claude/bench/data")
TWOBIT = BENCH / "hub-build" / "2bit"
KENT = BENCH / "kent-bin"
WORK = BENCH / "kmer_poly_work" / "kmersets"


def part_of(codes: np.ndarray, parts: int) -> np.ndarray:
    return (mix64(codes) >> np.uint64(64 - parts.bit_length() + 1)).astype(np.int64)


def gdir(genome: str, k: int) -> pathlib.Path:
    return WORK / f"{genome}.k{k}"


def build(genome: str, k: int, parts: int) -> dict:
    assert parts & (parts - 1) == 0, "parts must be a power of two"
    out = gdir(genome, k)
    if (out / "meta.json").exists():
        return json.loads((out / "meta.json").read_text())
    out.mkdir(parents=True, exist_ok=True)
    tb = TWOBIT / f"{genome}.2bit"
    names = [l.split("\t")[0] for l in subprocess.run([str(KENT / "twoBitInfo"), str(tb), "stdout"], check=True,
                                                       capture_output=True, text=True).stdout.splitlines()]
    raw = [open(out / f"raw_{p:03d}.u64", "wb") for p in range(parts)]
    t0, total = time.time(), 0
    for name in names:
        fa = subprocess.run([str(KENT / "twoBitToFa"), f"-seq={name}", str(tb), "stdout"], check=True,
                            capture_output=True).stdout
        seq = b"".join(fa.split(b"\n")[1:])
        codes = canonical_kmers(encode(seq), k)
        total += codes.size
        if codes.size:
            pid = part_of(codes, parts)
            order = np.argsort(pid, kind="stable")
            codes, pid = codes[order], pid[order]
            bounds = np.searchsorted(pid, np.arange(parts + 1))
            for p in range(parts):
                if bounds[p + 1] > bounds[p]:
                    codes[bounds[p]:bounds[p + 1]].tofile(raw[p])
    for f in raw:
        f.close()
    distinct = 0
    for p in range(parts):
        r = out / f"raw_{p:03d}.u64"
        u = np.unique(np.fromfile(r, dtype=np.uint64))
        u.tofile(out / f"part_{p:03d}.u64")
        r.unlink()
        distinct += u.size
    meta = {"genome": genome, "k": k, "parts": parts, "sequences": len(names), "kmers": total, "distinct": distinct,
            "seconds": round(time.time() - t0, 1)}
    (out / "meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    return meta


def load_part(genome: str, k: int, p: int) -> np.ndarray:
    return np.fromfile(gdir(genome, k) / f"part_{p:03d}.u64", dtype=np.uint64)


def diff(ga: str, gb: str, k: int) -> dict:
    ma, mb = (json.loads((gdir(g, k) / "meta.json").read_text()) for g in (ga, gb))
    assert ma["parts"] == mb["parts"]
    parts = ma["parts"]
    out = WORK / f"diff_{ga}__{gb}.k{k}"
    out.mkdir(parents=True, exist_ok=True)
    amb, bma = [], []
    for p in range(parts):
        a, b = load_part(ga, k, p), load_part(gb, k, p)
        amb.append(np.setdiff1d(a, b, assume_unique=True))
        bma.append(np.setdiff1d(b, a, assume_unique=True))
    amb_a, bma_a = np.concatenate(amb), np.concatenate(bma)
    amb_a.tofile(out / "a_minus_b.u64")
    bma_a.tofile(out / "b_minus_a.u64")
    res = {"A": ga, "B": gb, "k": k, "|A|": ma["distinct"], "|B|": mb["distinct"], "|A-B|": int(amb_a.size),
           "|B-A|": int(bma_a.size)}
    res["D"] = res["|A-B|"] + res["|B-A|"]
    res["D/|A|"] = res["D"] / res["|A|"]
    (out / "diff.json").write_text(json.dumps(res, indent=1) + "\n")
    return res


def self_test() -> int:
    # partitioning keeps every k-mer exactly once and in a deterministic partition
    rng = np.random.default_rng(1)
    codes = rng.integers(0, 1 << 62, 100000, dtype=np.uint64)
    pid = part_of(codes, 64)
    ok = pid.min() >= 0 and pid.max() < 64 and np.array_equal(pid, part_of(codes.copy(), 64))
    print("kmerset self-test:", "ok" if ok else "FAILED")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["build", "diff", "self-test"])
    ap.add_argument("genomes", nargs="*")
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--parts", type=int, default=64)
    a = ap.parse_args()
    if a.cmd == "self-test":
        return self_test()
    if a.cmd == "build":
        for g in a.genomes:
            print(json.dumps(build(g, a.k, a.parts)), flush=True)
    else:
        print(json.dumps(diff(a.genomes[0], a.genomes[1], a.k)), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
