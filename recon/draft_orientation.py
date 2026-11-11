"""How much of a draft assembly is written in the other orientation from its final assembly, and what that does to a
strand-specific k-mer method.

    python -m recon.draft_orientation orient <final> <draft> [--out JSON]     # contigs and bases same / flipped
    python -m recon.draft_orientation extdiff <final> <draft> [--out JSON]    # extended-syncmer difference, as written
                                                                              # and after canonicalization

orient (round-1 attack L7 F10, adjudication L7 F10). 31-mers are sampled by a strand-independent content filter
(mix64 of the canonical code, 1 in 1,024; reverse-complement palindromes dropped). A sampled k-mer is kept only if it
occurs exactly once in the final assembly, and its strand there (forward code == canonical code or not) is recorded.
Each draft contig is then called
    same      if >= 90% of its kept sampled k-mers are on the same strand as in the final assembly,
    flipped   if >= 90% are on the other strand,
    mixed     otherwise, and
    unplaced  if it has no kept sampled k-mer,
and the bases of each class are summed. Adapted from the attacker's orient.py, which held whole chromosomes in numpy;
here records are encoded in overlapping chunks so memory stays small.

extdiff. km-peeler's extended syncmers (recon/kmpeeler_run.py ext: ibltseq syncmers -k 31 -m 11 -g 20, deduplicated,
uppercased) are compared as written (LC_ALL=C comm -3) and after each string is replaced by min(s, rc(s))
(recon/canon_lines.py, then sort -u). Which positions are syncmers and the rightward extension both depend on the
strand, so canonicalizing the strings cannot remove the orientation effect.
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
from prototype.kmers import encode, mix64  # noqa: E402
from recon import kmerset  # noqa: E402

K = 31
SAMPLE_MASK = 1023
CHUNK = 4_000_000
WORK = kmerset.BENCH / "kmer_poly_work"
OWED = WORK / "round1_owed" / "orientation"


def records(genome: str):
    """(name, sequence) for every sequence of a hub genome, via twoBitToFa."""
    tb = kmerset.TWOBIT / f"{genome}.2bit"
    names = [l.split("\t")[0] for l in subprocess.run([str(kmerset.KENT / "twoBitInfo"), str(tb), "stdout"],
                                                       check=True, capture_output=True, text=True).stdout.splitlines()]
    for name in names:
        fa = subprocess.run([str(kmerset.KENT / "twoBitToFa"), f"-seq={name}", str(tb), "stdout"], check=True,
                            capture_output=True).stdout
        yield name, b"".join(fa.split(b"\n")[1:])


def sampled(seq: bytes) -> tuple[np.ndarray, np.ndarray]:
    """Canonical codes of the sampled k-mers of seq and whether each is written on its canonical (forward) strand."""
    cans, fwds = [], []
    for s0 in range(0, max(len(seq) - K + 1, 0), CHUNK):
        c = encode(seq[s0:s0 + CHUNK + K - 1])
        n = c.size - K + 1
        if n <= 0:
            continue
        bad = np.concatenate([[0], np.cumsum(c == 4)])
        ok = (bad[K:] - bad[:-K]) == 0
        cc = np.where(c == 4, 0, c).astype(np.uint64)
        fwd = np.zeros(n, np.uint64)
        rev = np.zeros(n, np.uint64)
        for j in range(K):
            fwd = (fwd << np.uint64(2)) | cc[j:j + n]
            rev |= (np.uint64(3) - cc[j:j + n]) << np.uint64(2 * j)
        can = np.minimum(fwd, rev)
        keep = ok & ((mix64(can) & np.uint64(SAMPLE_MASK)) == 0) & (fwd != rev)
        cans.append(can[keep])
        fwds.append((fwd == can)[keep])
    if not cans:
        return np.empty(0, np.uint64), np.empty(0, bool)
    return np.concatenate(cans), np.concatenate(fwds)


def orient(final: str, draft: str) -> dict:
    t0 = time.time()
    cs, fs = [], []
    for _, s in records(final):
        c, f = sampled(s)
        cs.append(c)
        fs.append(f)
    can, isf = np.concatenate(cs), np.concatenate(fs)
    u, idx, cnt = np.unique(can, return_index=True, return_counts=True)
    single = cnt == 1
    fin_codes, fin_fwd = u[single], isf[idx[single]]
    counts = {c: 0 for c in ("same", "flipped", "mixed", "unplaced")}
    bases = dict(counts)
    for _, s in records(draft):
        c, f = sampled(s)
        pos = np.searchsorted(fin_codes, c)
        pos = np.minimum(pos, max(fin_codes.size - 1, 0))
        hit = fin_codes[pos] == c
        agree = int(np.sum(hit & (fin_fwd[pos] == f)))
        dis = int(np.sum(hit & (fin_fwd[pos] != f)))
        cls = ("unplaced" if agree + dis == 0 else "same" if agree >= 0.9 * (agree + dis)
               else "flipped" if dis >= 0.9 * (agree + dis) else "mixed")
        counts[cls] += 1
        bases[cls] += len(s)
    tot = sum(bases.values())
    return {"final": final, "draft": draft, "k": K, "sample": f"mix64(canonical) & {SAMPLE_MASK} == 0, no palindromes",
            "final_sampled_kmers": int(can.size), "final_single_copy_sampled": int(fin_codes.size),
            "draft_contigs": counts, "draft_bases": bases, "draft_total_bases": tot,
            "fraction_of_bases": {c: round(b / tot, 4) for c, b in bases.items()},
            "seconds": round(time.time() - t0, 1)}


def _sh(cmd: str) -> str:
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True, executable="/bin/bash").stdout


def extdiff(final: str, draft: str) -> dict:
    src = WORK / f"kmp_ext_{final}__{draft}"
    OWED.mkdir(parents=True, exist_ok=True)
    canon = pathlib.Path(__file__).resolve().parent / "canon_lines.py"
    t0 = time.time()
    res = {"final": final, "draft": draft, "source": str(src)}
    for g in (final, draft):
        txt = src / f"{g}.ext.txt"
        _sh(f"LC_ALL=C sort -c -u {txt}")  # as written: must already be sorted and unique
        res[f"ext_{g}"] = int(_sh(f"wc -l < {txt}"))
        out = OWED / f"{g}.ext.canon.txt"
        if not out.exists():
            _sh(f"{sys.executable} {canon} < {txt} | LC_ALL=C sort -u -S 1G -T {OWED} > {out}.tmp && mv {out}.tmp {out}")
        res[f"ext_canon_{g}"] = int(_sh(f"wc -l < {out}"))
    a, b = src / f"{final}.ext.txt", src / f"{draft}.ext.txt"
    res["D_ext_as_written"] = int(_sh(f"LC_ALL=C comm -3 {a} {b} | wc -l"))
    ca, cb = OWED / f"{final}.ext.canon.txt", OWED / f"{draft}.ext.canon.txt"
    res["D_ext_canonicalized"] = int(_sh(f"LC_ALL=C comm -3 {ca} {cb} | wc -l"))
    res["seconds"] = round(time.time() - t0, 1)
    return res


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["orient", "extdiff"])
    ap.add_argument("final")
    ap.add_argument("draft")
    ap.add_argument("--out", type=pathlib.Path)
    a = ap.parse_args()
    res = orient(a.final, a.draft) if a.cmd == "orient" else extdiff(a.final, a.draft)
    print(json.dumps(res, indent=1))
    if a.out:
        a.out.write_text(json.dumps(res, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
