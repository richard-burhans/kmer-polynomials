"""Shibuya, Belazzougui & Kucherov 2022 (km-peeler, pinned b8ec39b; recon/build_kmpeeler.sh) on a real genome pair.

    python -m recon.kmpeeler_run kmers <A> <B> --parts 0,1 --ratio 1.0 --trials 20      # IBLT of canonical k-mers
    python -m recon.kmpeeler_run ext <A> <B> --z 11 --ratio 1.0 --trials 5 [--seqs N]  # extended syncmers (Sec. 3.5)

kmers: the "straightforward" k-mer IBLT of Section 3.1/3.3 on the SAME canonical k-mer sets as recon/kp10.py (selected
kmerset partitions, written as text). n (the IBLT threshold, ibltseq build -n) = ratio x the true difference D, with
r = 3 hash functions and epsilon 0 (their defaults; m = n * c_r cells). Each trial uses a different seed (-s). A trial
succeeds when `list` returns exactly the true A\\B and B\\A.

ext: their proposed recover-all method (Section 3.5): extended syncmers of length 2k - z (ibltseq syncmers -k k -m z
-g k-z), deduplicated, IBLT difference, then the k-mers of the listed extended syncmers with those shared by both listed
sets removed. This runs on the FASTA AS WRITTEN (strand-specific, as km-peeler is). The result is compared with the
exact canonical and forward-strand k-mer differences: missing and spurious k-mers are reported. Needs a km-peeler
build with a 2k - z payload (build_kmpeeler.sh <dir> minimizers <k> <z>).

Sketch size is the IBLT file km-peeler writes.
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
from prototype.kmers import canonical_kmers, encode  # noqa: E402
from recon import kmerset  # noqa: E402

W = pathlib.Path("/home/rcb112/project/claude/bench/data/kmer_poly_work")
KMP_K = W / "km-peeler"          # payload k (configured syncmers k z: payload k)
KMP_EXT = W / "km-peeler-ext"    # payload 2k - z
ACGT = np.frombuffer(b"ACGT", dtype=np.uint8)


def to_text(codes: np.ndarray, path: pathlib.Path, k: int = 31) -> None:
    """Write k-mer codes as ACGT lines (s_0 = most significant digit), in chunks."""
    with open(path, "wb") as fh:
        for s in range(0, codes.size, 2_000_000):
            c = codes[s:s + 2_000_000]
            digits = ((c[:, None] >> (np.uint64(2) * np.arange(k - 1, -1, -1, dtype=np.uint64))) & np.uint64(3))
            rows = np.empty((c.size, k + 1), np.uint8)
            rows[:, :k] = ACGT[digits.astype(np.intp)]
            rows[:, k] = ord("\n")
            fh.write(rows.tobytes())


def from_text(lines: list[bytes], k: int = 31) -> np.ndarray:
    if not lines:
        return np.empty(0, np.uint64)
    arr = np.frombuffer(b"".join(lines), dtype=np.uint8).reshape(-1, k)
    lut = np.full(256, 0, np.uint64)
    lut[ord("C")], lut[ord("G")], lut[ord("T")] = 1, 2, 3
    v = lut[arr]
    out = np.zeros(arr.shape[0], np.uint64)
    for t in range(k):
        out = (out << np.uint64(2)) | v[:, t]
    return np.sort(out)


def ibltseq(binary: pathlib.Path, *args: str, stdin=None) -> bytes:
    return subprocess.run([str(binary / "ibltseq"), *args], check=True, capture_output=True, stdin=stdin).stdout


def build_diff_list(binary, a_txt, b_txt, n, seed, length, tmp):
    t0 = time.time()
    ibltseq(binary, "build", "-i", str(a_txt), "-n", str(n), "-r", "3", "-s", str(seed), "-l", str(length),
            "-o", str(tmp / "a.iblt"))
    ibltseq(binary, "build", "-i", str(b_txt), "-n", str(n), "-r", "3", "-s", str(seed), "-l", str(length),
            "-o", str(tmp / "b.iblt"))
    t1 = time.time()
    ibltseq(binary, "diff", "-i", str(tmp / "a.iblt"), "-j", str(tmp / "b.iblt"), "-o", str(tmp / "d.iblt"))
    out = ibltseq(binary, "list", "-i", str(tmp / "d.iblt")).splitlines()
    t2 = time.time()
    ia = [l[2:] for l in out if l.startswith(b"i,")]
    jb = [l[2:] for l in out if l.startswith(b"j,")]
    return ia, jb, (tmp / "a.iblt").stat().st_size, t1 - t0, t2 - t1


def run_kmers(ga, gb, k, parts, ratio, trials):
    tmp = W / f"kmp_kmers_{ga}__{gb}"
    tmp.mkdir(exist_ok=True)
    a = np.concatenate([kmerset.load_part(ga, k, p) for p in parts])
    b = np.concatenate([kmerset.load_part(gb, k, p) for p in parts])
    want_amb = np.sort(np.setdiff1d(a, b, assume_unique=True))
    want_bma = np.sort(np.setdiff1d(b, a, assume_unique=True))
    D = want_amb.size + want_bma.size
    to_text(a, tmp / "a.txt", k)
    to_text(b, tmp / "b.txt", k)
    n = max(int(round(ratio * D)), 1)
    ok = 0
    sizes, tb, td = [], [], []
    for t in range(trials):
        ia, jb, size, t_build, t_dec = build_diff_list(KMP_K, tmp / "a.txt", tmp / "b.txt", n, 1000 + t, k, tmp)
        got_a, got_b = from_text(ia, k), from_text(jb, k)
        success = np.array_equal(got_a, want_amb) and np.array_equal(got_b, want_bma)
        ok += success
        sizes.append(size); tb.append(t_build); td.append(t_dec)
        print(f"trial {t}: listed {got_a.size}+{got_b.size} of {D}; exact {success}; {size} bytes", flush=True)
    return {"mode": "kmers", "A": ga, "B": gb, "k": k, "parts": parts, "elements_A": int(a.size),
            "elements_B": int(b.size), "D": int(D), "n": n, "ratio": ratio, "trials": trials, "successes": ok,
            "failure_rate": 1 - ok / trials, "iblt_bytes": int(np.median(sizes)),
            "bytes_per_difference": float(np.median(sizes)) / D, "seconds_build_both": round(float(np.median(tb)), 2),
            "seconds_diff_list": round(float(np.median(td)), 2)}


def genome_fasta(g: str, nseqs: int | None, out: pathlib.Path) -> pathlib.Path:
    tb = kmerset.TWOBIT / f"{g}.2bit"
    if nseqs:
        names = [l.split("\t")[0] for l in subprocess.run([str(kmerset.KENT / "twoBitInfo"), str(tb), "stdout"],
                                                           check=True, capture_output=True, text=True).stdout.splitlines()]
        (out.parent / f"{g}.seqs").write_text("\n".join(names[:nseqs]) + "\n")
        subprocess.run([str(kmerset.KENT / "twoBitToFa"), f"-seqList={out.parent / (g + '.seqs')}", str(tb), str(out)],
                       check=True)
    else:
        subprocess.run([str(kmerset.KENT / "twoBitToFa"), str(tb), str(out)], check=True)
    return out


def run_ext(ga, gb, k, z, ratio, trials, nseqs):
    tmp = W / f"kmp_ext_{ga}__{gb}"
    tmp.mkdir(exist_ok=True)
    ext_txt = {}
    kmers_fwd, kmers_can = {}, {}
    for g in (ga, gb):
        fa = genome_fasta(g, nseqs, tmp / f"{g}.fa")
        raw = tmp / f"{g}.ext.raw"
        with open(raw, "wb") as fh:
            subprocess.run([str(KMP_EXT / "ibltseq"), "syncmers", "-k", str(k), "-m", str(z), "-g", str(k - z),
                            "-i", str(fa)], check=True, stdout=fh)
        ext_txt[g] = tmp / f"{g}.ext.txt"
        subprocess.run(f"LC_ALL=C sort -u -S 2G -T {tmp} {raw} | tr a-z A-Z | LC_ALL=C sort -u -S 2G -T {tmp} > {ext_txt[g]}",
                       shell=True, check=True)
        raw.unlink()
        # exact k-mer sets for scoring: forward strand (as written) and canonical
        recs, cur = [], []
        for l in open(fa, "rb"):
            if l.startswith(b">"):
                if cur:
                    recs.append(b"".join(cur))
                cur = []
            else:
                cur.append(l.strip())
        if cur:
            recs.append(b"".join(cur))
        # per record, so no k-mer spans two records
        kmers_can[g] = np.unique(np.concatenate([canonical_kmers(encode(r), k) for r in recs]))
        kmers_fwd[g] = None
    # exact element difference of extended-syncmer sets (for sizing n)
    D_ext = int(subprocess.run(f"LC_ALL=C comm -3 {ext_txt[ga]} {ext_txt[gb]} | wc -l", shell=True, check=True,
                               capture_output=True, text=True).stdout)
    want = np.union1d(np.setdiff1d(kmers_can[ga], kmers_can[gb]), np.setdiff1d(kmers_can[gb], kmers_can[ga]))
    n = max(int(round(ratio * D_ext)), 1)
    res_trials = []
    for t in range(trials):
        ia, jb, size, t_build, t_dec = build_diff_list(KMP_EXT, ext_txt[ga], ext_txt[gb], n, 2000 + t, 2 * k - z, tmp)
        listed = len(ia) + len(jb)
        # Section 3.5: k-mers of the listed extended syncmers, minus those shared by both listed sets (canonicalized
        # here for scoring against the canonical truth)
        def kmers_of(lines):
            if not lines:
                return np.empty(0, np.uint64)
            return np.unique(np.concatenate([canonical_kmers(encode(l.strip()), k) for l in lines]))
        ka, kb = kmers_of(ia), kmers_of(jb)
        got = np.union1d(np.setdiff1d(ka, kb), np.setdiff1d(kb, ka))
        missing = int(np.setdiff1d(want, got).size)
        spurious = int(np.setdiff1d(got, want).size)
        res_trials.append({"listed": listed, "complete": listed == D_ext, "missing": missing, "spurious": spurious,
                           "bytes": size, "seconds_build_both": round(t_build, 2), "seconds_diff_list": round(t_dec, 2)})
        print(f"trial {t}: listed {listed} of D_ext {D_ext}; k-mers missing {missing}, spurious {spurious}; {size} B",
              flush=True)
    return {"mode": "ext", "A": ga, "B": gb, "k": k, "z": z, "seqs": nseqs, "D_ext": D_ext, "D_kmers_canonical": int(want.size),
            "n": n, "ratio": ratio, "trials": res_trials}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["kmers", "ext"])
    ap.add_argument("A")
    ap.add_argument("B")
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--z", type=int, default=11)
    ap.add_argument("--parts", default="0")
    ap.add_argument("--ratio", type=float, default=1.0)
    ap.add_argument("--trials", type=int, default=20)
    ap.add_argument("--seqs", type=int)
    ap.add_argument("--out", type=pathlib.Path)
    a = ap.parse_args()
    if a.mode == "kmers":
        res = run_kmers(a.A, a.B, a.k, [int(x) for x in a.parts.split(",")], a.ratio, a.trials)
    else:
        res = run_ext(a.A, a.B, a.k, a.z, a.ratio, a.trials, a.seqs)
    print(json.dumps(res))
    if a.out:
        a.out.write_text(json.dumps(res, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
