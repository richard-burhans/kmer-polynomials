"""KP-1 Eisenstein revisit: substitution spectrum of real cannabis differences, scored by the two encodings.

Reads pairwise synNet MAFs whose reference is cs10 NC_044378.1 (bench/data/wfi_laila/NC_044378.1/bare/<label>.maf,
the WF-I inputs) and counts aligned substitutions by class:
  transition (A-G, C-T), transversion A-C/G-T ("tv13": distance 3 under Eisenstein), transversion A-T/C-G
  ("complement": distance 4 under both encodings).
Pairs reported: cs10 vs each query, and the two haplotypes of one plant against each other (KOMPa vs KOMPb), compared
only at cs10 positions where both align (projected through cs10 coordinates). Case-blind; gaps and N skipped. With
--upper, only columns where every compared base is uppercase (outside soft-masked repeats).

Squared distance per class: quarter map (A=1, C=i, G=-i, T=-1) 2 / 2 / 4; Eisenstein (A=1, C=w, G=-w, T=-1) 1 / 3 / 4.

    kp1_tstv.py MAF_DIR --ref-len 71238074 --labels KOMPa KOMPb [--upper] [--out out.json]
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

CODE = np.full(256, 255, dtype=np.uint8)
for ch, v in zip("ACGTacgt", (0, 1, 2, 3, 0, 1, 2, 3)):
    CODE[ord(ch)] = v
# class of an unordered pair of distinct codes A0 C1 G2 T3
CLS = np.full((4, 4), -1, dtype=np.int8)
for x, y, c in ((0, 2, 0), (1, 3, 0), (0, 1, 1), (2, 3, 1), (0, 3, 2), (1, 2, 2)):
    CLS[x, y] = CLS[y, x] = c
NAMES = ("transition", "tv13", "complement")
COST = {"quarter": (2, 2, 4), "eisenstein": (1, 3, 4)}


def project(maf: pathlib.Path, ref_len: int):
    """Arrays over cs10 coordinates: query base code (255 = unaligned/gap/N) and uppercase flags for ref and query."""
    q = np.full(ref_len, 255, dtype=np.uint8)
    up = np.zeros(ref_len, dtype=bool)
    rbase = np.full(ref_len, 255, dtype=np.uint8)
    rows = []
    with open(maf, "rb") as fh:
        for line in fh:
            if line.startswith(b"s "):
                rows.append(line.split())
            elif not line.strip() and rows:
                ref, qry = rows[0], rows[1]
                start = int(ref[2])
                rt = np.frombuffer(ref[6], dtype=np.uint8)
                qt = np.frombuffer(qry[6], dtype=np.uint8)
                keep = rt != ord("-")
                pos = start + np.arange(int(keep.sum()))
                rk, qk = rt[keep], qt[keep]
                rbase[pos] = CODE[rk]
                q[pos] = CODE[qk]  # '-' maps to 255
                up[pos] = (rk < ord("a")) & (qk < ord("a")) & (qk != ord("-"))
                rows = []
    return rbase, q, up


def spectrum(a: np.ndarray, b: np.ndarray, ok: np.ndarray) -> dict:
    m = ok & (a < 4) & (b < 4)
    diff = m & (a != b)
    cls = CLS[a[diff], b[diff]]
    counts = [int((cls == c).sum()) for c in range(3)]
    n = sum(counts)
    out = {"compared": int(m.sum()), "substitutions": n, "counts": dict(zip(NAMES, counts))}
    if n:
        ts = counts[0]
        tv = counts[1] + counts[2]
        out["ts_tv"] = ts / tv if tv else None
        out["fraction"] = {k: v / n for k, v in zip(NAMES, counts)}
        for enc, cost in COST.items():
            out[f"mean_cost_{enc}"] = sum(c * w for c, w in zip(counts, cost)) / n
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("maf_dir", type=pathlib.Path)
    ap.add_argument("--ref-len", type=int, required=True)
    ap.add_argument("--labels", nargs=2, default=["KOMPa", "KOMPb"])
    ap.add_argument("--out", type=pathlib.Path)
    a = ap.parse_args()
    la, lb = a.labels
    ra, qa, ua = project(a.maf_dir / f"{la}.maf", a.ref_len)
    rb, qb, ub = project(a.maf_dir / f"{lb}.maf", a.ref_len)
    ref = np.where(ra < 4, ra, rb)
    res = {"ref_len": a.ref_len, "labels": [la, lb], "costs": COST}
    for upper in (False, True):
        tag = "uppercase_only" if upper else "all_columns"
        res[tag] = {
            f"cs10_vs_{la}": spectrum(ref, qa, ua if upper else np.ones_like(ua)),
            f"cs10_vs_{lb}": spectrum(ref, qb, ub if upper else np.ones_like(ub)),
            f"{la}_vs_{lb}_at_shared_cs10_positions": spectrum(qa, qb, (ua & ub) if upper else np.ones_like(ua)),
        }
    js = json.dumps(res, indent=1)
    if a.out:
        a.out.write_text(js + "\n")
    print(js)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
