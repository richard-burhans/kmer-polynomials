"""KP-10 on k-mers AS WRITTEN: no canonicalization, no 2-bit packing, bucket keys taken from the hash value itself.

    python -m recon.aswritten build <genome> --k K --mode sym|fwd --seed-x SX --seed-lam SL
    python -m recon.aswritten run <genomeA> <genomeB> --k K --mode sym|fwd --seed-x SX --seed-lam SL --bits B
                                  [--m 64] [--c 2] [--max-rounds R] [--out result.json]
    python -m recon.aswritten verify-canonical <genomeA> <genomeB> --seed-x SX --seed-lam SL   # k = 31 truth
    python -m recon.aswritten verify-kmc <genomeA> <genomeB> --k K --seed-x SX --seed-lam SL --amb F --bma F
    python -m recon.aswritten self-test

Why. recon/kp10.py reconciles min(fwd, rc) canonical codes with buckets keyed by mix64 of that code, so its orientation
invariance comes from canonicalization, and any element hash would have had it (round-1 completeness critic, gap 1).
Here each k-mer is hashed from the sequence bytes in the orientation it is written (hash/kp7e_stream.c), and:
  sym  element = S = (F + lam)(R + lam), strand-symmetric, so S(s) = S(rc s) without choosing a strand;
  fwd  element = F = sum s_t x^t, strand-specific: the control, whose difference counts orientation changes.
The element IS the hash value (the code column of SideE holds it too), and the bucket key is the value shifted to the
top of 64 bits, so partitions (top 6 bits) and buckets (top `bits` bits) are ranges of the value; nothing else is
hashed. Because the rolling hash reads bases, not a packed code, any k works (k > 32 included).

Each genome's values are streamed into 64 partitions by value range and deduplicated (the set of distinct values).
`run` reconciles partition by partition with recon/recon61.py (power sums S_1..S_{m+c} per bucket, Newton + Pade,
c check sums, split on overflow) and verifies every partition against the exact set difference of the VALUE sets.
`verify-*` then checks the recovered values against an independent truth computed on k-mers: the canonical 2-bit
difference of recon/kmerset.py (k <= 32) or a KMC 3 difference (any k), hashed with the same x and lam.

Bucket keys from S: a bucket's elements share their top bits, which changes nothing in the decode (power sums of
distinct field elements; each side flags only its own elements as roots, then c check sums and the count delta must
agree). The verification against independent truth is what shows no bucket decoded wrongly.
"""

from __future__ import annotations

import argparse
import ctypes
import json
import pathlib
import subprocess
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from recon import kmerset, kp7e_hash, recon61  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
STREAM = ROOT / "hash" / "kp7e_stream"
WORK = kmerset.BENCH / "kmer_poly_work" / "round1_owed" / "aswritten"
PART_BITS = 6


def stream_binary() -> pathlib.Path:
    subprocess.run(["make", "-s", "-C", str(ROOT / "hash"), "kp7e_stream"], check=True)
    return STREAM


def gdir(genome: str, k: int, mode: str, sx: int, sl: int) -> pathlib.Path:
    return WORK / f"{genome}.k{k}.{mode}.x{sx:x}.l{sl:x}"


PASS_PARTS = 4  # partitions per pass: a pass holds about 1/16 of a genome's k-mer values (~0.4 GB at 900 Mb)


def _stream_into(fasta_cmd: list[str] | None, fasta_path: pathlib.Path | None, out: pathlib.Path, k: int, mode: str,
                 sx: int, sl: int) -> dict:
    """Hash every k-mer as written and store each partition's distinct values. The genome is streamed once per pass
    of PASS_PARTS partitions (hashing is seconds; holding all values, or writing them all to disk, is not)."""
    out.mkdir(parents=True, exist_ok=True)
    binary = str(stream_binary())
    distinct, info = 0, None
    for lo in range(0, 1 << PART_BITS, PASS_PARTS):
        hi = lo + PASS_PARTS
        args = [binary, str(k), mode, str(sx), str(sl), str(PART_BITS), str(lo), str(hi)]
        if fasta_cmd:
            src = subprocess.Popen(fasta_cmd, stdout=subprocess.PIPE)
            r = subprocess.run(args, stdin=src.stdout, capture_output=True)
            src.stdout.close()
            assert src.wait() == 0, f"{fasta_cmd[0]} failed"
        else:
            with open(fasta_path, "rb") as fh:
                r = subprocess.run(args, stdin=fh, capture_output=True)
        assert r.returncode == 0, r.stderr.decode()
        pinfo = json.loads(r.stderr.decode())
        vals = np.frombuffer(r.stdout, dtype=np.uint64)
        assert vals.size == pinfo["written"]
        if info is None:
            info = {key: pinfo[key] for key in ("k", "mode", "x", "lam", "records", "bases", "kmers")}
        assert all(info[key] == pinfo[key] for key in info), "passes disagree"
        part = (vals >> np.uint64(61 - PART_BITS)).astype(np.int64)
        for p in range(lo, hi):
            u = np.unique(vals[part == p])
            u.tofile(out / f"part_{p:03d}.u64")
            distinct += int(u.size)
        del vals, part, r
    info["distinct"] = distinct
    info["passes"] = (1 << PART_BITS) // PASS_PARTS
    return info


def build(genome: str, k: int, mode: str, sx: int, sl: int) -> dict:
    out = gdir(genome, k, mode, sx, sl)
    if (out / "meta.json").exists():
        return json.loads((out / "meta.json").read_text())
    t0 = time.time()
    tb = kmerset.TWOBIT / f"{genome}.2bit"
    info = _stream_into([str(kmerset.KENT / "twoBitToFa"), str(tb), "stdout"], None, out, k, mode, sx, sl)
    info.update({"genome": genome, "seed_x": sx, "seed_lam": sl, "parts": 1 << PART_BITS,
                 "seconds": round(time.time() - t0, 1)})
    (out / "meta.json").write_text(json.dumps(info, indent=1) + "\n")
    return info


def load_part(genome: str, k: int, mode: str, sx: int, sl: int, p: int) -> np.ndarray:
    return np.fromfile(gdir(genome, k, mode, sx, sl) / f"part_{p:03d}.u64", dtype=np.uint64)


def side(values: np.ndarray) -> recon61.SideE:
    """Element = value = reported identity; key = the value at the top of 64 bits (values < 2^61)."""
    return recon61.SideE.build(values, values, values << np.uint64(64 - 61))


def part_buckets(part: int, bits: int) -> np.ndarray:
    span = 1 << (bits - PART_BITS)
    return np.arange(part * span, (part + 1) * span, dtype=np.int64)


def value_difference(ga: str, gb: str, k: int, mode: str, sx: int, sl: int) -> int:
    """|A \\ B| + |B \\ A| of the two VALUE sets (exact, partition by partition)."""
    d = 0
    for part in range(1 << PART_BITS):
        a, b = load_part(ga, k, mode, sx, sl, part), load_part(gb, k, mode, sx, sl, part)
        d += a.size + b.size - 2 * np.intersect1d(a, b, assume_unique=True).size
    return int(d)


def bits_for(D: int, m: int) -> int:
    """Smallest bucket count 2^bits (at least the 64 partitions) with an expected difference per bucket <= m/2: the
    rule prototype/run_kompa.py and the KP-10 runs used. It reads the true D, as those runs did."""
    bits = PART_BITS
    while D / (1 << bits) > m / 2:
        bits += 1
    return bits


def run(ga: str, gb: str, k: int, mode: str, sx: int, sl: int, bits: int | None, m: int, c: int,
        max_rounds: int = 12, parts: int = 1 << PART_BITS) -> dict:
    bits_rule = "given"
    if bits is None:
        bits, bits_rule = bits_for(value_difference(ga, gb, k, mode, sx, sl), m), "auto: bits_for(exact value D, m)"
    assert bits >= PART_BITS
    p = kp7e_hash.params(k, sx, sl)
    out = WORK / f"run_{ga}__{gb}.k{k}.{mode}.x{sx:x}.l{sl:x}.b{bits}.m{m}"
    out.mkdir(parents=True, exist_ok=True)
    tot = {"A": ga, "B": gb, "k": k, "mode": mode, "seed_x": sx, "seed_lam": sl, "x": int(p.x), "lam": int(p.lam),
           "bits": bits, "bits_rule": bits_rule, "m": m, "c": c, "max_rounds": max_rounds, "parts": parts, "bucket_key": "top bits of the value",
           "elements_A": 0, "elements_B": 0, "D_true_values": 0, "D_recovered": 0, "exact": True,
           "sketch_elements": 0, "rounds_max": 0, "splits": 0, "unresolved_buckets": 0, "seconds_reconcile": 0.0,
           "round_log": []}
    amb_all, bma_all = [], []
    for part in range(parts):
        a, b = load_part(ga, k, mode, sx, sl, part), load_part(gb, k, mode, sx, sl, part)
        t0 = time.time()
        res = recon61.reconcile(side(a), side(b), bits, m, c, buckets=part_buckets(part, bits), max_rounds=max_rounds,
                                stop_at_limit=True)
        dt = time.time() - t0
        want_amb = np.setdiff1d(a, b, assume_unique=True)
        want_bma = np.setdiff1d(b, a, assume_unique=True)
        ok = res.unresolved == 0 and np.array_equal(res.a_minus_b, want_amb) and np.array_equal(res.b_minus_a, want_bma)
        tot["exact"] &= bool(ok)
        tot["elements_A"] += int(a.size)
        tot["elements_B"] += int(b.size)
        tot["D_true_values"] += int(want_amb.size + want_bma.size)
        tot["D_recovered"] += int(res.a_minus_b.size + res.b_minus_a.size)
        tot["sketch_elements"] += res.sketch_elements
        tot["rounds_max"] = max(tot["rounds_max"], res.rounds)
        tot["splits"] += res.splits
        tot["unresolved_buckets"] += res.unresolved
        tot["seconds_reconcile"] += dt
        tot["round_log"].append(res.log)
        amb_all.append(res.a_minus_b)
        bma_all.append(res.b_minus_a)
        print(f"part {part}: |A|={a.size} |B|={b.size} D={want_amb.size + want_bma.size} "
              f"recovered={res.a_minus_b.size + res.b_minus_a.size} exact={ok} rounds={res.rounds} "
              f"splits={res.splits} unresolved={res.unresolved} {dt:.1f}s", flush=True)
    np.concatenate(amb_all).tofile(out / "a_minus_b.u64")
    np.concatenate(bma_all).tofile(out / "b_minus_a.u64")
    tot["sketch_bytes"] = tot["sketch_elements"] * 8
    tot["bytes_per_difference"] = tot["sketch_bytes"] / max(tot["D_recovered"], 1)
    tot["seconds_reconcile"] = round(tot["seconds_reconcile"], 1)
    tot["recovered_dir"] = str(out)
    (out / "result.json").write_text(json.dumps(tot, indent=1) + "\n")
    return tot


def _recovered(rundir: pathlib.Path) -> tuple[np.ndarray, np.ndarray]:
    return (np.fromfile(rundir / "a_minus_b.u64", dtype=np.uint64),
            np.fromfile(rundir / "b_minus_a.u64", dtype=np.uint64))


def _compare(name: str, rec_amb, rec_bma, true_amb_vals, true_bma_vals, extra: dict) -> dict:
    ta, tb = np.unique(true_amb_vals), np.unique(true_bma_vals)
    ra, rb = np.unique(rec_amb), np.unique(rec_bma)
    res = {"truth": name, **extra,
           "true_A_minus_B_kmers": int(len(true_amb_vals)), "true_B_minus_A_kmers": int(len(true_bma_vals)),
           "true_distinct_hash_values": int(ta.size + tb.size),
           "recovered_A_minus_B": int(ra.size), "recovered_B_minus_A": int(rb.size),
           "recovered_equals_truth": bool(np.array_equal(ra, ta) and np.array_equal(rb, tb)),
           "missed": int(np.setdiff1d(ta, ra).size + np.setdiff1d(tb, rb).size),
           "spurious": int(np.setdiff1d(ra, ta).size + np.setdiff1d(rb, tb).size)}
    return res


def verify_canonical(ga: str, gb: str, rundir: pathlib.Path, sx: int, sl: int, k: int = 31) -> dict:
    """Truth = the canonical 2-bit k-mer difference of recon/kmerset.py, hashed with the same x and lam."""
    d = kmerset.WORK / f"diff_{ga}__{gb}.k{k}"
    amb = np.fromfile(d / "a_minus_b.u64", dtype=np.uint64)
    bma = np.fromfile(d / "b_minus_a.u64", dtype=np.uint64)
    p = kp7e_hash.params(k, sx, sl)
    ra, rb = _recovered(rundir)
    return _compare(f"recon/kmerset.py canonical difference {d.name}", ra, rb, kp7e_hash.sym_hash(amb, p, k),
                    kp7e_hash.sym_hash(bma, p, k), {"rundir": str(rundir)})


def s_of_strings(kmers: list[bytes], k: int, sx: int, sl: int) -> np.ndarray:
    """S of k-mer strings by the C reference (kp7e_direct + kp7e_sym), any k."""
    so, p = kp7e_hash.lib(), kp7e_hash.params(k, sx, sl)
    out = np.empty(len(kmers), np.uint64)
    F, R = ctypes.c_uint64(), ctypes.c_uint64()
    for i, s in enumerate(kmers):
        assert len(s) == k and so.kp7e_direct(ctypes.byref(p), s, ctypes.byref(F), ctypes.byref(R)) == 1, s
        out[i] = so.kp7e_sym(ctypes.byref(p), F.value, R.value)
    return out


def read_kmc_dump(path: pathlib.Path) -> list[bytes]:
    with open(path, "rb") as fh:
        return [line.split(b"\t")[0] for line in fh if line.strip()]


def verify_kmc(rundir: pathlib.Path, k: int, sx: int, sl: int, amb_dump: pathlib.Path, bma_dump: pathlib.Path) -> dict:
    """Truth = KMC 3 canonical k-mer differences (kmc_tools simple ... kmers_subtract, then transform dump)."""
    amb, bma = read_kmc_dump(amb_dump), read_kmc_dump(bma_dump)
    ra, rb = _recovered(rundir)
    return _compare("KMC 3 canonical k-mer difference", ra, rb, s_of_strings(amb, k, sx, sl),
                    s_of_strings(bma, k, sx, sl), {"rundir": str(rundir), "amb_dump": str(amb_dump),
                                                   "bma_dump": str(bma_dump)})


def _values_of_fasta(fa: pathlib.Path, k: int, mode: str, sx: int, sl: int, tmp: pathlib.Path) -> np.ndarray:
    out = tmp / f"{fa.stem}.{mode}.k{k}"
    _stream_into(None, fa, out, k, mode, sx, sl)
    return np.concatenate([np.fromfile(out / f"part_{p:03d}.u64", dtype=np.uint64) for p in range(1 << PART_BITS)])


def self_test() -> int:
    """On a random sequence with N runs and soft-masking, and its reverse complement:
    sym values (as written) == numpy S of the canonical codes (k = 31); sym values are equal across orientation, fwd
    values are not; fwd values == numpy F of the forward codes; k = 63 sym values == the C reference on strings and
    are equal across orientation; every partition holds only its value range."""
    from prototype.kmers import canonical_kmers, encode
    rng = np.random.default_rng(17)
    seq = bytearray(rng.choice(np.frombuffer(b"ACGTacgt", np.uint8), 30000).tobytes())
    for s0 in (5000, 17000):
        seq[s0:s0 + 40] = b"N" * 40
    seq = bytes(seq)
    comp = bytes.maketrans(b"ACGTacgtN", b"TGCAtgcaN")
    sx, sl = 0x5EED, 0x1A3D
    ok = True
    with tempfile.TemporaryDirectory() as td:
        tmp = pathlib.Path(td)
        rcp = lambda b: b.translate(comp)[::-1]
        (tmp / "fw.fa").write_bytes(b">r1\n" + seq[:12000] + b"\n>r2\n" + seq[12000:] + b"\n")
        (tmp / "rc.fa").write_bytes(b">r2\n" + rcp(seq[12000:]) + b"\n>r1\n" + rcp(seq[:12000]) + b"\n")
        p31 = kp7e_hash.params(31, sx, sl)
        sym_fw = np.unique(_values_of_fasta(tmp / "fw.fa", 31, "sym", sx, sl, tmp))
        sym_rc = np.unique(_values_of_fasta(tmp / "rc.fa", 31, "sym", sx, sl, tmp))
        codes = np.concatenate([canonical_kmers(encode(seq[:12000]), 31), canonical_kmers(encode(seq[12000:]), 31)])
        want = np.unique(kp7e_hash.sym_hash(np.unique(codes), p31, 31))
        c1 = np.array_equal(sym_fw, want)
        c2 = np.array_equal(sym_fw, sym_rc)
        fwd_fw = np.unique(_values_of_fasta(tmp / "fw.fa", 31, "fwd", sx, sl, tmp))
        fwd_rc = np.unique(_values_of_fasta(tmp / "rc.fa", 31, "fwd", sx, sl, tmp))
        fcodes = _forward_codes(seq[:12000], 31) + _forward_codes(seq[12000:], 31)
        c3 = np.array_equal(fwd_fw, np.unique(kp7e_hash.poly(np.array(fcodes, np.uint64), p31, 31)))
        c4 = np.intersect1d(fwd_fw, fwd_rc).size < fwd_fw.size // 100
        sym63_fw = np.unique(_values_of_fasta(tmp / "fw.fa", 63, "sym", sx, sl, tmp))
        sym63_rc = np.unique(_values_of_fasta(tmp / "rc.fa", 63, "sym", sx, sl, tmp))
        up = seq.upper()
        strs = [part[i:i + 63] for part in (up[:12000], up[12000:]) for i in range(len(part) - 62)
                if b"N" not in part[i:i + 63]]
        c5 = np.array_equal(sym63_fw, np.unique(s_of_strings(strs, 63, sx, sl))) and np.array_equal(sym63_fw, sym63_rc)
        out = tmp / "fw.sym.k31"
        c6 = all(np.all((np.fromfile(out / f"part_{p:03d}.u64", dtype=np.uint64) >> np.uint64(61 - PART_BITS)) == p)
                 for p in range(1 << PART_BITS))
    for label, c in (("k=31 sym as written == S of canonical codes", c1), ("k=31 sym equal across orientation", c2),
                     ("k=31 fwd == F of forward codes", c3), ("k=31 fwd differs across orientation", c4),
                     ("k=63 sym == C reference on strings, equal across orientation", c5),
                     ("partitions are value ranges", c6)):
        print(f"  {'ok  ' if c else 'FAIL'} {label}")
        ok &= bool(c)
    print("aswritten self-test:", "ok" if ok else "FAILED")
    return 0 if ok else 1


def _forward_codes(s: bytes, k: int) -> list[int]:
    lut = {ord(ch): v for ch, v in zip("ACGTacgt", (0, 1, 2, 3, 0, 1, 2, 3))}
    out = []
    for i in range(len(s) - k + 1):
        w = s[i:i + k]
        if any(ch not in lut for ch in w):
            continue
        c = 0
        for ch in w:
            c = (c << 2) | lut[ch]
        out.append(c)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["build", "run", "verify-canonical", "verify-kmc", "self-test"])
    ap.add_argument("genomes", nargs="*")
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--mode", choices=["sym", "fwd"], default="sym")
    ap.add_argument("--seed-x", type=lambda v: int(v, 0), default=0x5EED)
    ap.add_argument("--seed-lam", type=lambda v: int(v, 0), default=0x1A3D)
    ap.add_argument("--bits", type=lambda v: None if v == "auto" else int(v), default=8, help="int or auto")
    ap.add_argument("--m", type=int, default=64)
    ap.add_argument("--c", type=int, default=2)
    ap.add_argument("--max-rounds", type=int, default=12)
    ap.add_argument("--rundir", type=pathlib.Path)
    ap.add_argument("--amb", type=pathlib.Path)
    ap.add_argument("--bma", type=pathlib.Path)
    ap.add_argument("--out", type=pathlib.Path)
    a = ap.parse_args()
    if a.cmd == "self-test":
        return self_test()
    if a.cmd == "build":
        for g in a.genomes:
            print(json.dumps(build(g, a.k, a.mode, a.seed_x, a.seed_lam)), flush=True)
        return 0
    if a.cmd == "run":
        res = run(a.genomes[0], a.genomes[1], a.k, a.mode, a.seed_x, a.seed_lam, a.bits, a.m, a.c, a.max_rounds)
    elif a.cmd == "verify-canonical":
        res = verify_canonical(a.genomes[0], a.genomes[1], a.rundir, a.seed_x, a.seed_lam, a.k)
    else:
        res = verify_kmc(a.rundir, a.k, a.seed_x, a.seed_lam, a.amb, a.bma)
    print(json.dumps(res, indent=1))
    if a.out:
        a.out.write_text(json.dumps(res, indent=1) + "\n")
    ok = res.get("exact", True) and res.get("recovered_equals_truth", True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
