"""KP-13: all-offset FFT correlation as a diagonal finder, tested against real chains.

The correlation of a target window T with a query segment Q at offset o is sum_t conj(T_t) Q_{t+o}. With the complex
encoding its real part scores +1 per match, -1 per complementary base and 0 otherwise (IDEAS.md KP-1, KP-13). One FFT
of length >= |T| + |Q| gives every offset. A real alignment along diagonal d (query position minus target position)
shows as a peak at the offset for d.

Test pair: cs10 NC_044378.1 (target) against KOMPa.chr7 (query), whose chains are known (the 12-of-12 vgp rerun,
bench/data/multiz-lowmem/cache/KOMPa.over.chain). For each sampled target window:
  - choose the strand and the median diagonal carrying the most chain-aligned bases in the window;
  - correlate the window against a query segment of +/- `--flank` around that diagonal, on that strand (for a '-'
    chain, the reverse complement of the query, whose coordinates the chain uses);
  - call peaks (local maxima of the robust z-score of the real part, non-maximum suppression +/- `--tol`);
  - score them against truth: the window's chain blocks on that strand, clustered by diagonal (+/- `--tol`).
    Precision is the share of peaks within `--tol` of any chain diagonal; recall is the share of chain-diagonal
    clusters carrying >= `--min-weight` aligned bases that have a peak within `--tol`.
Variants: `cplx` (match - complement); `cplx_mask` (soft-masked bases zeroed); `ind4` (four 0/1 indicator sequences,
the exact match count, minus its median). Negative control: the same window against a same-sized query segment that
no chain block of the window maps to. Cost: FFT time per window, and lastz on the same window and segment.

    kp13_fft_diag.py TARGET_FA QUERY_FA CHAIN [--windows 40] [--win 100000] [--flank 500000] [--lastz /usr/bin/lastz]
                     [--genome-search 3] [--out results.json]
"""
from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import tempfile
import time

import numpy as np

from common import read_chains, read_fasta_one, revcomp, to_complex

ZS = (6.0, 10.0)


def nextpow2(n: int) -> int:
    return 1 << (int(n - 1).bit_length())


def xcorr_real(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Real part of sum_t conj(a_t) b_{t+o} for o = 0 .. len(b) - len(a) (full overlap)."""
    n = nextpow2(len(a) + len(b))
    fa = np.fft.fft(a.astype(np.complex64), n)
    fb = np.fft.fft(b.astype(np.complex64), n)
    r = np.fft.ifft(np.conj(fa) * fb)
    return r.real[: len(b) - len(a) + 1].astype(np.float64)


def xcorr_ind4(ta: np.ndarray, qb: np.ndarray) -> np.ndarray:
    """Exact match count per offset: sum over bases of the indicator correlations (real FFTs)."""
    n = nextpow2(len(ta) + len(qb))
    tot = None
    for c in range(4):
        x = (ta == c).astype(np.float32)
        y = (qb == c).astype(np.float32)
        r = np.fft.irfft(np.conj(np.fft.rfft(x, n)) * np.fft.rfft(y, n), n)[: len(qb) - len(ta) + 1]
        tot = r if tot is None else tot + r
    return np.rint(tot)


def robust_z(s: np.ndarray) -> np.ndarray:
    med = np.median(s)
    mad = np.median(np.abs(s - med)) * 1.4826
    return (s - med) / (mad if mad > 0 else 1.0)


def peaks(z: np.ndarray, thr: float, tol: int, cap: int = 200) -> list[int]:
    idx = np.nonzero(z >= thr)[0]
    if len(idx) == 0:
        return []
    order = idx[np.argsort(-z[idx])]
    taken: list[int] = []
    for o in order:
        if all(abs(o - p) > tol for p in taken):
            taken.append(int(o))
            if len(taken) >= cap:
                break
    return taken


def window_truth(chains, ws: int, we: int):
    """{strand: [(diagonal, overlap_bases)]} for chain blocks overlapping [ws, we)."""
    out = {"+": [], "-": []}
    for ch in chains:
        for t, q, size in ch["blocks"]:
            lo, hi = max(t, ws), min(t + size, we)
            if hi > lo:
                out[ch["strand"]].append((q - t, hi - lo))
    return out


def cluster(diags, tol: int):
    if not diags:
        return []
    diags = sorted(diags)
    cl = [[diags[0][0], diags[0][1], diags[0][0]]]  # [weighted-mean-ish anchor, weight, last]
    for d, w in diags[1:]:
        if d - cl[-1][2] <= tol:
            cl[-1][1] += w
            cl[-1][2] = d
        else:
            cl.append([d, w, d])
    return [((c[0] + c[2]) // 2, c[1], c[0], c[2]) for c in cl]  # (centre, weight, lo, hi)


def near(d: int, clusters, tol: int) -> bool:
    return any(c[2] - tol <= d <= c[3] + tol for c in clusters)


def run_lastz(lastz: str, T: str, Q: str, strand: str, tmp: pathlib.Path) -> tuple[float, int]:
    tf, qf = tmp / "t.fa", tmp / "q.fa"
    tf.write_text(">t\n" + T + "\n")
    qf.write_text(">q\n" + Q + "\n")
    t0 = time.perf_counter()
    r = subprocess.run([lastz, str(tf), str(qf), "--strand=plus", "--format=general:score", "--output=/dev/stdout"],
                       capture_output=True, text=True)
    dt = time.perf_counter() - t0
    n = max(len([l for l in r.stdout.splitlines() if l and not l.startswith("#")]), 0)
    return dt, n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("target", type=pathlib.Path)
    ap.add_argument("query", type=pathlib.Path)
    ap.add_argument("chain", type=pathlib.Path)
    ap.add_argument("--windows", type=int, default=40)
    ap.add_argument("--win", type=int, default=100_000)
    ap.add_argument("--flank", type=int, default=500_000)
    ap.add_argument("--tol", type=int, default=50)
    ap.add_argument("--min-weight", type=int, default=1000)
    ap.add_argument("--min-covered", type=int, default=5000)
    ap.add_argument("--lastz", default="/usr/bin/lastz")
    ap.add_argument("--lastz-windows", type=int, default=5)
    ap.add_argument("--genome-search", type=int, default=3)
    ap.add_argument("--rng", type=int, default=7)
    ap.add_argument("--out", type=pathlib.Path)
    a = ap.parse_args()
    rng = np.random.default_rng(a.rng)
    tname, tc, tl = read_fasta_one(a.target)
    qname, qc, ql = read_fasta_one(a.query)
    qrc, qlrc = revcomp(qc, ql)
    Qs = {"+": (qc, ql), "-": (qrc, qlrc)}
    chains = read_chains(a.chain, tname, qname)
    letters = np.frombuffer(b"ACGTN", dtype=np.uint8)
    out = {"target": tname, "query": qname, "chains": len(chains), "win": a.win, "flank": a.flank, "tol": a.tol,
           "min_weight": a.min_weight, "windows": []}
    starts = np.linspace(0, len(tc) - a.win, a.windows * 3).astype(int)
    rng.shuffle(starts)
    tmp = pathlib.Path(tempfile.mkdtemp(dir="/home/rcb112/project/claude/bench/data/kmer_poly_work/spectral"))
    done = lz = 0
    for ws in starts:
        if done >= a.windows:
            break
        we = ws + a.win
        truth = window_truth(chains, ws, we)
        cov = {s: sum(w for _, w in truth[s]) for s in "+-"}
        strand = max(cov, key=cov.get)
        if cov[strand] < a.min_covered:
            continue
        ds = np.repeat([d for d, _ in truth[strand]], [w for _, w in truth[strand]])
        dmed = int(np.median(ds))
        qcodes, qlow = Qs[strand]
        s0 = max(0, ws + dmed - a.flank)
        s1 = min(len(qcodes), ws + dmed + a.win + a.flank)
        if s1 - s0 < a.win * 2:
            continue
        T, Tl = tc[ws:we], tl[ws:we]
        Q, Ql = qcodes[s0:s1], qlow[s0:s1]
        clusters = cluster([(d, w) for d, w in truth[strand] if s0 - ws <= d <= s1 - we], a.tol)
        heavy = [c for c in clusters if c[1] >= a.min_weight]
        rec = {"ws": int(ws), "strand": strand, "covered": int(cov[strand]), "clusters": len(clusters),
               "heavy_clusters": len(heavy), "heavy_weight": int(sum(c[1] for c in heavy)), "variants": {}}
        # negative control: a far segment of the same length on the same strand with no chain diagonal of this window
        for _ in range(20):
            n0 = int(rng.integers(0, len(qcodes) - (s1 - s0)))
            if not any(n0 - ws <= d <= n0 + (s1 - s0) - we for d, _ in truth[strand]):
                break
        N, Nl = qcodes[n0:n0 + (s1 - s0)], qlow[n0:n0 + (s1 - s0)]
        for var in ("cplx", "cplx_mask", "ind4"):
            t0 = time.perf_counter()
            if var == "ind4":
                sc, scn = xcorr_ind4(T, Q), xcorr_ind4(T, N)
            else:
                m = var == "cplx_mask"
                tv = to_complex(T, Tl, m)
                sc, scn = xcorr_real(tv, to_complex(Q, Ql, m)), xcorr_real(tv, to_complex(N, Nl, m))
            dt = (time.perf_counter() - t0) / 2
            z, zn = robust_z(sc), robust_z(scn)
            vr = {"seconds_per_correlation": dt, "max_z": float(z.max()), "neg_max_z": float(zn.max())}
            for thr in ZS:
                pk = peaks(z, thr, a.tol)
                diag = [s0 + o - ws for o in pk]
                hits = [d for d in diag if near(d, clusters, a.tol)]
                found = [c for c in heavy if any(c[2] - a.tol <= d <= c[3] + a.tol for d in diag)]
                vr[f"z{thr:g}"] = {"peaks": len(pk), "precision": len(hits) / len(pk) if pk else None,
                                   "recall_heavy": len(found) / len(heavy) if heavy else None,
                                   "recall_heavy_weight": (sum(c[1] for c in found) / rec["heavy_weight"])
                                   if heavy else None,
                                   "neg_peaks": len(peaks(zn, thr, a.tol))}
            rec["variants"][var] = vr
        if lz < a.lastz_windows:
            Ts = letters[T].tobytes().decode()
            Qq = letters[Q].tobytes().decode()
            dt, nal = run_lastz(a.lastz, Ts, Qq, strand, tmp)
            rec["lastz"] = {"seconds": dt, "alignments": nal}
            lz += 1
        out["windows"].append(rec)
        done += 1
        print(f"window {done} ws={ws} strand={strand} cov={cov[strand]} heavy={len(heavy)} "
              f"cplx z6 P={rec['variants']['cplx']['z6']['precision']} R={rec['variants']['cplx']['z6']['recall_heavy']}",
              flush=True)
    # whole-chromosome search: a window against the WHOLE query, both strands, overlap-save in 4 Mb chunks
    gs = []
    for rec in out["windows"][: a.genome_search]:
        ws, we = rec["ws"], rec["ws"] + a.win
        truth = window_truth(chains, ws, we)
        tv = to_complex(tc[ws:we], tl[ws:we], False)
        g = {"ws": ws, "strands": {}}
        t0 = time.perf_counter()
        for strand in "+-":
            qcodes, qlow = Qs[strand]
            step = 4_000_000
            scores = []
            for c0 in range(0, len(qcodes) - a.win, step):
                seg = to_complex(qcodes[c0:c0 + step + a.win], qlow[c0:c0 + step + a.win], False)
                scores.append(xcorr_real(tv, seg)[:step])
            sc = np.concatenate(scores)
            z = robust_z(sc)
            pk = peaks(z, 10.0, a.tol, cap=50)
            cl = cluster(truth[strand], a.tol)
            diag = [o - ws for o in pk]
            g["strands"][strand] = {"peaks_z10": len(pk), "on_chain_diagonal": sum(near(d, cl, a.tol) for d in diag),
                                    "max_z": float(z.max()), "chain_bases": int(sum(w for _, w in truth[strand]))}
        g["seconds_both_strands"] = time.perf_counter() - t0
        gs.append(g)
        print(f"genome search ws={ws}: {g}", flush=True)
    out["genome_search"] = gs
    js = json.dumps(out, indent=1)
    if a.out:
        a.out.write_text(js + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
