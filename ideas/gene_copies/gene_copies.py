#!/usr/bin/env python3
"""Find the copies and pseudogenes of a reference gene family in a genome: a k-mer scan against lastz and miniprot.

Pilot family: the cannabinoid oxidocyclases (THCAS, CBDAS, CBCAS). All three are single-exon genes, so there are no
exons to treat separately; a copy is scored as full or partial along the reference, and as intact or pseudogene by
its open reading frame (k-mer, lastz, truth) or by miniprot's frameshift and stop-codon calls.

    gene_copies.py verify-refs                         # re-fetch every accession from NCBI and compare with refs/
    gene_copies.py stage   --genome G --src PATH       # .fa/.fa.gz/.2bit -> DATA/G/genome.fa (cs10 is linked)
    gene_copies.py kmer    --genome G [--jobs 6]       # k-mer scans: k = 11, 15, 21, 31, 15 with 1 mismatch, spaced
    gene_copies.py lastz   --genome G [--jobs 8]       # lastz DNA baseline, per sequence, in parallel
    gene_copies.py miniprot --genome G [--jobs 16]     # miniprot protein-to-genome baseline
    gene_copies.py blast   --genome G [--jobs 16]      # blastn (dc-megablast) + tblastn: the truth search, not scored
    gene_copies.py truth   --genome cs10               # truth set: blast + RefSeq family genes + Grassa 2021 Table S9
    gene_copies.py score   --genome cs10               # recall/precision per method and identity bin -> DATA/cs10/score.json
    gene_copies.py census  --genome G                  # per-genome copy counts by class, every method
    gene_copies.py self-test

Large outputs go to DATA (bench/data/gene_copies/<genome>/). Methods write loci_<method>.tsv in one shared schema.
"""
from __future__ import annotations

import argparse
import gzip
import json
import multiprocessing
import os
import pathlib
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
BENCH = pathlib.Path("/home/rcb112/project/claude/bench")
DATA = BENCH / "data" / "gene_copies"
KENT = BENCH / "data" / "kent-bin"
MINIPROT = DATA / "env" / "bin" / "miniprot"
CG = BENCH / "repos" / "cannabis-genome"
CS10_FA = BENCH / "data" / "kmer_poly_work" / "cs10.fa"
CS10_GFF = BENCH / "data" / "cs10" / "GCF_900626175.2_cs10_genomic.gff.gz"
FAMILY_GENES = CG / "data" / "processed" / "cs10_tandem_genes.tsv"
GRASSA_S9 = CG / "data" / "processed" / "grassa_s9_refseq_overlap.tsv"
REFS_NT = HERE / "refs" / "synthases_cds.fa"
REFS_AA = HERE / "refs" / "synthases_protein.fa"
MANIFEST = HERE / "refs" / "manifest.json"

#: seed variants: name -> (span, care positions or None for contiguous, mismatches allowed in the reference ball)
PALINDROMIC_MASK = "1110111011101110111011101110111"  # ideas/algebra/test_spaced.c, ideas/eisenstein (KP-14)
SEEDS = {
    "k11": (11, None, 0),
    "k15": (15, None, 0),
    "k21": (21, None, 0),
    "k31": (31, None, 0),
    "k15d1": (15, None, 1),
    "sp31w24": (31, PALINDROMIC_MASK, 0),
}
#: a k-mer locus is CALLED when its hits on one diagonal band cover at least this many reference bases (fixed
#: before scoring; 150 bp = 9% of a 1,635-bp synthase, enough to call a partial copy)
MIN_COVER = 150
MIN_HITS = 3
CLUSTER_GAP = 1000   # bp between consecutive hits of one cluster
DIAG_BAND = 150      # bp of diagonal drift (indels) allowed inside one segment
JOIN_GAP = 15000     # collinear segments of one reference closer than this are one copy (an insertion splits it)
FULL_FRACTION = 0.90  # reference span covered for a "full" copy
ORF_FRACTION = 0.95   # longest ATG..stop ORF, as a fraction of the reference protein, for "intact"
ORF_MAX = 1.15
BINS = [(">=95", 95.0, 101.0), ("90-95", 90.0, 95.0), ("80-90", 80.0, 90.0), ("<80", -1.0, 80.0)]
#: the <80 bin split further: blastn finds a DNA hit at 70-80%, under 70%, or none (tblastn/RefSeq only)
SUBBINS = [("70-80", 70.0, 80.0), ("<70", 0.0, 70.0), ("protein-only", None, None)]

LUT = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    LUT[_b] = _i
    LUT[_b + 32] = _i
COMP = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")
CODON = {a + b + c: aa for (a, b, c), aa in zip(((x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"),
                                                  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}


# ---------------------------------------------------------------- sequence I/O
def read_fasta(path: pathlib.Path):
    """Yield (name, sequence bytes) from a FASTA file, gzipped or not."""
    op = gzip.open if str(path).endswith(".gz") else open
    name, chunks = None, []
    with op(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(chunks)
                name, chunks = line[1:].split()[0].decode(), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, b"".join(chunks)


def revcomp(s: bytes) -> bytes:
    return s.translate(COMP)[::-1]


def translate(s: bytes) -> str:
    t = s.decode().upper()
    return "".join(CODON.get(t[i:i + 3], "X") for i in range(0, len(t) - 2, 3))


def load_refs() -> dict[str, bytes]:
    return {n: s.upper() for n, s in read_fasta(REFS_NT)}


# ---------------------------------------------------------------- k-mer keys
def seed_positions(span: int, mask: str | None) -> np.ndarray:
    if mask is None:
        return np.arange(span)
    assert len(mask) == span
    return np.array([i for i, c in enumerate(mask) if c == "1"])


def window_keys(codes: np.ndarray, span: int, pos: np.ndarray):
    """Canonical keys of every window: (key, forward_is_canonical, valid). Pure.

    Forward key packs codes[i + p] over care positions p; the reverse-complement key packs the complement of
    codes[i + span - 1 - p], i.e. the same mask read on the other strand. min() of the two is identical for a window and
    its reverse complement for ANY mask (the packed-code twin of ideas/eisenstein/spaced_e.h's (F + lam)(R + lam)).
    A window with a non-ACGT base anywhere in its span is invalid, as in spaced_e.h."""
    n = len(codes) - span + 1
    if n <= 0:
        e = np.zeros(0, np.uint64)
        return e, np.zeros(0, bool), np.zeros(0, bool)
    cv = np.where(codes < 4, codes, 0).astype(np.uint64)
    c3 = np.uint64(3) - cv
    fw = np.zeros(n, np.uint64)
    rc = np.zeros(n, np.uint64)
    two = np.uint64(2)
    for p in pos:
        fw <<= two
        fw |= cv[p:p + n]
        q = span - 1 - p
        rc <<= two
        rc |= c3[q:q + n]
    bad = np.concatenate(([0], np.cumsum(codes == 4, dtype=np.int64)))
    valid = (bad[span:span + n] - bad[:n]) == 0
    fwd = fw <= rc
    return np.minimum(fw, rc), fwd, valid


def ref_index(refs: dict[str, bytes], seed: str):
    """Sorted arrays (key, ref id, ref pos, ref-forward-is-canonical) over every reference window, with the
    Hamming ball of radius d around each forward window when the seed allows mismatches. Pure."""
    span, mask, d = SEEDS[seed]
    pos = seed_positions(span, mask)
    w = len(pos)
    keys, rid, rpos, ror = [], [], [], []
    for i, (name, s) in enumerate(refs.items()):
        codes = LUT[np.frombuffer(s, np.uint8)]
        cv = codes.astype(np.uint64)
        n = len(codes) - span + 1
        fw = np.zeros(n, np.uint64)
        for p in pos:
            fw = (fw << np.uint64(2)) | cv[p:p + n]
        ball = [fw]
        if d >= 1:
            for j in range(w):
                sh = np.uint64(2 * (w - 1 - j))
                for sub in (1, 2, 3):
                    ball.append(fw ^ (np.uint64(sub) << sh))
        assert d <= 1, "radius > 1 not implemented"
        for f in ball:
            rcf = rc_key(f, w)
            keys.append(np.minimum(f, rcf)); ror.append(f <= rcf)
            rid.append(np.full(n, i, np.int32)); rpos.append(np.arange(n, dtype=np.int32))
    k = np.concatenate(keys); o = np.argsort(k, kind="stable")
    return k[o], np.concatenate(rid)[o], np.concatenate(rpos)[o], np.concatenate(ror)[o]


def rc_key(f: np.ndarray, w: int) -> np.ndarray:
    """Reverse complement of packed keys of weight w (2 bits per base). For a palindromic mask this is the key of the
    reverse-complement window; the reference side is built only from such masks or contiguous seeds. Pure."""
    out = np.zeros_like(f)
    x = f.copy()
    for _ in range(w):
        out = (out << np.uint64(2)) | (np.uint64(3) - (x & np.uint64(3)))
        x >>= np.uint64(2)
    return out


def scan_sequence(args):
    """All reference-k-mer hits on one genome sequence: arrays (gpos, rid, rpos, same_strand). Chunked to bound RAM."""
    path, name, seed, idx = args
    K, R, P, O = idx
    span, mask, _ = SEEDS[seed]
    pos = seed_positions(span, mask)
    seq = fetch_seq(path, name)
    codes = LUT[np.frombuffer(seq, np.uint8)]
    out = []
    step = 8_000_000
    for a in range(0, len(codes), step):
        c = codes[a:a + step + span - 1]
        key, fwd, valid = window_keys(c, span, pos)
        j = np.searchsorted(K, key)
        j[j == len(K)] = len(K) - 1
        hit = np.nonzero((K[j] == key) & valid)[0]
        if not len(hit):
            continue
        hk = key[hit]
        lo = np.searchsorted(K, hk, "left"); hi = np.searchsorted(K, hk, "right")
        cnt = hi - lo
        g = np.repeat(hit, cnt)
        first = np.repeat(lo - np.concatenate(([0], np.cumsum(cnt)[:-1])), cnt)
        jj = first + np.arange(cnt.sum())
        out.append((g + a, R[jj], P[jj], np.repeat(fwd[hit], cnt) == O[jj]))
    if not out:
        z = np.zeros(0, np.int64)
        return name, (z, z.astype(np.int32), z.astype(np.int32), z.astype(bool))
    return name, tuple(np.concatenate([o[i] for o in out]) for i in range(4))


def fetch_seq(fa: pathlib.Path, name: str) -> bytes:
    """One record by its .fai offsets (no full-file read per worker)."""
    for line in pathlib.Path(str(fa) + ".fai").read_text().splitlines():
        f = line.split("\t")
        if f[0] == name:
            n, off, lb, lw = int(f[1]), int(f[2]), int(f[3]), int(f[4])
            with open(fa, "rb") as fh:
                fh.seek(off)
                raw = fh.read(n + (n // lb + 1) * (lw - lb))
            return raw.replace(b"\n", b"").replace(b"\r", b"")[:n]
    raise KeyError(name)


def segments_from_hits(g, rid, rpos, same, span):
    """Chain hits by diagonal into segments (rid, strand, gs, ge, rs, re, cover, nhits). Pure.

    Per reference and strand: sort by genome position, cut where consecutive hits are > CLUSTER_GAP apart; inside a
    cluster sort by diagonal (g - r on +, g + r on -) and cut where the diagonal jumps by > DIAG_BAND."""
    segs = []
    for r in np.unique(rid):
        for strand, sel in (("+", same), ("-", ~same)):
            m = (rid == r) & sel
            if not m.any():
                continue
            gg, rr = g[m], rpos[m].astype(np.int64)
            o = np.argsort(gg, kind="stable"); gg, rr = gg[o], rr[o]
            cuts = np.nonzero(np.diff(gg) > CLUSTER_GAP)[0] + 1
            for cg, cr in zip(np.split(gg, cuts), np.split(rr, cuts)):
                dg = cg - cr if strand == "+" else cg + cr
                o2 = np.argsort(dg, kind="stable")
                dcut = np.nonzero(np.diff(dg[o2]) > DIAG_BAND)[0] + 1
                for part in np.split(o2, dcut):
                    if len(part) < MIN_HITS:
                        continue
                    sg, sr = cg[part], cr[part]
                    cover = union_len(np.sort(sr), span)
                    segs.append((int(r), strand, int(sg.min()), int(sg.max()) + span, int(sr.min()),
                                 int(sr.max()) + span, cover, len(part)))
    return segs


def union_len(starts_sorted: np.ndarray, span: int) -> int:
    """Bases covered by intervals [s, s + span) over sorted starts. Pure."""
    if not len(starts_sorted):
        return 0
    gaps = np.diff(starts_sorted)
    return int(span + np.minimum(gaps, span).sum())


# ---------------------------------------------------------------- loci, shared by every method
LOCUS_COLS = ["chrom", "start", "end", "strand", "ref", "ref_start", "ref_end", "ref_cover", "identity", "n_segments",
              "full", "status", "orf_aa", "detail"]


def join_collinear(segs):
    """Merge segments (chrom, rid, strand, gs, ge, rs, re, cover, ident, detail) of one reference and strand that are
    collinear along the reference and closer than JOIN_GAP in the genome. Pure."""
    out = []
    key = lambda s: (s[0], s[1], s[2], s[3])
    for s in sorted(segs, key=key):
        if out:
            p = out[-1]
            same = p[0] == s[0] and p[1] == s[1] and p[2] == s[2]
            near = s[3] - p[4] <= JOIN_GAP
            colinear = (s[5] >= p[6] - 100) if s[2] == "+" else (s[6] <= p[5] + 100)
            if same and near and colinear:
                w1, w2 = p[7], s[7]
                ident = (p[8] * w1 + s[8] * w2) / (w1 + w2) if (p[8] is not None and s[8] is not None) else (p[8] or s[8])
                out[-1] = (p[0], p[1], p[2], p[3], max(p[4], s[4]), min(p[5], s[5]), max(p[6], s[6]), p[7] + s[7],
                           ident, p[9] + s[9], (p[10] or []) + (s[10] or []))
                continue
        out.append(tuple(s[:9]) + (s[9] if len(s) > 9 else 1, list(s[10]) if len(s) > 10 and s[10] else []))
    return out


def merge_refs(chains, refnames):
    """Merge chains of different references that overlap in the genome into one locus; the reference with the most
    covered bases names it. Pure. Returns dicts in LOCUS_COLS order (status filled later)."""
    chains = sorted(chains, key=lambda c: (c[0], c[3]))
    loci = []
    for c in chains:
        if loci and loci[-1]["chrom"] == c[0] and c[3] < loci[-1]["end"]:
            L = loci[-1]
            L["end"] = max(L["end"], c[4])
            L["members"].append(c)
            continue
        loci.append({"chrom": c[0], "start": c[3], "end": c[4], "members": [c]})
    out = []
    for L in loci:
        best = max(L["members"], key=lambda c: c[7] * (c[8] if c[8] is not None else 100.0))
        covers = {}
        for c in L["members"]:
            covers[refnames[c[1]]] = max(covers.get(refnames[c[1]], 0), c[7])
        out.append({"chrom": L["chrom"], "start": L["start"], "end": L["end"], "strand": best[2],
                    "ref": refnames[best[1]], "ref_start": best[5], "ref_end": best[6], "ref_cover": best[7],
                    "identity": None if best[8] is None else round(best[8], 2), "n_segments": best[9],
                    "detail": ";".join(f"{k}:{v}" for k, v in sorted(covers.items())) + (
                        (";" + ",".join(best[10])) if best[10] else "")})
    return out


def longest_orf(seq: bytes) -> int:
    """Length in codons of the longest ATG..stop ORF on the given strand (stop excluded). Pure."""
    best = 0
    s = seq.upper()
    for f in range(3):
        aa = translate(s[f:])
        i = 0
        while True:
            m = aa.find("M", i)
            if m < 0:
                break
            stop = aa.find("*", m)
            if stop < 0:
                best = max(best, len(aa) - m)
                break
            best = max(best, stop - m)
            i = stop + 1
    return best


def classify(loci, genome: dict[str, bytes] | None, refs: dict[str, bytes], use_orf=True):
    """Fill full (reference span >= FULL_FRACTION) and, where the method has no frameshift calls, status by the ORF test
    on the implied gene region (the locus extended to both reference ends, plus 300 bp). Mutates and returns loci."""
    for L in loci:
        rl = len(refs[L["ref"]])
        L["full"] = (L["ref_end"] - L["ref_start"]) >= FULL_FRACTION * rl
        if use_orf and genome is not None:
            if L["strand"] == "+":
                a, b = L["start"] - L["ref_start"] - 300, L["end"] + (rl - L["ref_end"]) + 300
            else:
                a, b = L["start"] - (rl - L["ref_end"]) - 300, L["end"] + L["ref_start"] + 300
            s = genome[L["chrom"]][max(0, a):b]
            if L["strand"] == "-":
                s = revcomp(s)
            L["orf_aa"] = longest_orf(s)
            # an ORF much LONGER than the reference is read-through into an insertion, not an intact copy
            L["status"] = "intact" if ORF_FRACTION * (rl // 3 - 1) <= L["orf_aa"] <= ORF_MAX * (rl // 3 - 1) \
                else "pseudogene"
    return loci


def write_loci(path: pathlib.Path, loci):
    with open(path, "w") as fh:
        fh.write("\t".join(LOCUS_COLS) + "\n")
        for L in loci:
            fh.write("\t".join("" if L.get(c) is None else str(L.get(c)) for c in LOCUS_COLS) + "\n")


def read_loci(path: pathlib.Path):
    rows = []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            d = dict(zip(hdr, line.rstrip("\n").split("\t")))
            for c in ("start", "end", "ref_start", "ref_end", "ref_cover", "n_segments"):
                d[c] = int(d[c])
            d["identity"] = float(d["identity"]) if d["identity"] else None
            d["full"] = d["full"] == "True"
            rows.append(d)
    return rows


# ---------------------------------------------------------------- per-genome staging and methods
def gdir(genome: str) -> pathlib.Path:
    d = DATA / genome
    d.mkdir(parents=True, exist_ok=True)
    return d


def genome_fa(genome: str) -> pathlib.Path:
    p = gdir(genome) / "genome.fa"
    if not p.exists():
        raise SystemExit(f"⛔ {p} missing: run stage first")
    return p


def seq_names(fa: pathlib.Path) -> list[tuple[str, int]]:
    fai = pathlib.Path(str(fa) + ".fai")
    if not fai.exists():
        subprocess.run(["samtools", "faidx", str(fa)], check=True)
    return [(l.split("\t")[0], int(l.split("\t")[1])) for l in fai.read_text().splitlines()]


def load_genome(fa: pathlib.Path) -> dict[str, bytes]:
    return dict(read_fasta(fa))


def cpu() -> float:
    """User + system CPU seconds of this process and every reaped child (pool workers, lastz, miniprot)."""
    import resource
    a, b = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    return a.ru_utime + a.ru_stime + b.ru_utime + b.ru_stime


def timing(genome: str, method: str, seconds: float, threads: int, extra: dict | None = None, cpu_s: float | None = None):
    p = gdir(genome) / "timings.json"
    t = json.loads(p.read_text()) if p.exists() else {}
    t[method] = {"wall_s": round(seconds, 1), "threads": threads, **(extra or {})}
    if cpu_s is not None:
        t[method]["cpu_s"] = round(cpu_s, 1)
    p.write_text(json.dumps(t, indent=1, sort_keys=True))


def cmd_stage(a):
    out = gdir(a.genome) / "genome.fa"
    src = pathlib.Path(a.src)
    if out.exists():
        print(f"{out} exists"); return
    if src.suffix == ".2bit":
        subprocess.run([str(KENT / "twoBitToFa"), str(src), str(out)], check=True)
    elif src.suffix == ".gz":
        with gzip.open(src, "rb") as i, open(out, "wb") as o:
            while b := i.read(1 << 24):
                o.write(b)
    else:
        out.symlink_to(src)
    seq_names(out)
    print(f"staged {out}")


def cmd_kmer(a):
    fa = genome_fa(a.genome)
    refs = load_refs(); names = list(refs)
    seqs = [n for n, L in seq_names(fa) if L >= 1000]
    genome = None
    for seed in (a.seeds or list(SEEDS)):
        t0, c0 = time.time(), cpu()
        idx = ref_index(refs, seed)
        span = SEEDS[seed][0]
        segs = []
        nhits = 0
        with ProcessPoolExecutor(a.jobs, mp_context=multiprocessing.get_context("fork")) as ex:
            for name, (g, rid, rpos, same) in ex.map(scan_sequence, [(fa, n, seed, idx) for n in seqs]):
                nhits += len(g)
                for s in segments_from_hits(g, rid, rpos, same, span):
                    if s[6] >= MIN_COVER:
                        segs.append((name, s[0], s[1], s[2], s[3], s[4], s[5], s[6], None, 1, []))
        loci = merge_refs(join_collinear(segs), names)
        t_scan, c_scan = time.time() - t0, cpu() - c0
        if genome is None:
            genome = load_genome(fa)
        classify(loci, genome, refs)
        write_loci(gdir(a.genome) / f"loci_{seed}.tsv", loci)
        timing(a.genome, seed, t_scan, a.jobs, {"hits": int(nhits), "ref_keys": int(len(idx[0])), "loci": len(loci)}, c_scan)
        print(f"{seed}: {nhits} hits, {len(loci)} loci, {t_scan:.1f} s scan")


LASTZ_FMT = "general:name1,strand1,zstart1,end1,name2,strand2,zstart2+,end2+,size2,nmatch,nmismatch,score"


def lastz_one(args):
    fa, name, out = args
    pathlib.Path(str(out) + ".name").write_text(name + "\n")
    cmd = ["lastz", f"{fa}[multiple,subset={out}.name,unmask]", f"{REFS_NT}[multiple,unmask]", "--strand=both",
           "--ambiguous=iupac", f"--format={LASTZ_FMT}", f"--output={out}"]
    subprocess.run(cmd, check=True)
    return out


def cmd_lastz(a):
    fa = genome_fa(a.genome)
    d = gdir(a.genome) / "lastz_parts"; d.mkdir(exist_ok=True)
    seqs = [n for n, L in seq_names(fa) if L >= 1000]
    t0, c0 = time.time(), cpu()
    if a.parse_only:
        parts = [d / "whole.tsv"] if (d / "whole.tsv").exists() else [d / f"{i}.tsv" for i in range(len(seqs))]
    elif a.jobs == 1:
        # one lastz over the whole genome: the honest single-thread cost (the parallel driver below re-parses the
        # whole FASTA in every per-sequence job, which inflates its CPU time)
        out = d / "whole.tsv"
        subprocess.run(["lastz", f"{fa}[multiple,unmask]", f"{REFS_NT}[multiple,unmask]", "--strand=both",
                        "--ambiguous=iupac", f"--format={LASTZ_FMT}", f"--output={out}"], check=True)
        parts = [out]
    else:
        with ProcessPoolExecutor(a.jobs, mp_context=multiprocessing.get_context("fork")) as ex:
            parts = list(ex.map(lastz_one, [(fa, n, d / f"{i}.tsv") for i, n in enumerate(seqs)]))
    elapsed, c_el = time.time() - t0, cpu() - c0
    rows = []
    for p in parts:
        for line in open(p):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            rows.append(f)
    with open(gdir(a.genome) / "lastz.tsv", "w") as fh:
        fh.write("#" + LASTZ_FMT + "\n")
        for f in rows:
            fh.write("\t".join(f) + "\n")
    refs = load_refs(); names = list(refs)
    segs = []
    for f in rows:
        chrom, _, gs, ge, ref, rstrand, rs, re_, _, nm, nmm, _ = f
        cover = int(re_) - int(rs)
        if cover < MIN_COVER:
            continue
        ident = 100.0 * int(nm) / max(1, int(nm) + int(nmm))
        segs.append((chrom, names.index(ref), "+" if rstrand == "+" else "-", int(gs), int(ge), int(rs), int(re_),
                     cover, ident, 1, []))
    loci = merge_refs(join_collinear(segs), names)
    classify(loci, load_genome(fa), refs)
    write_loci(gdir(a.genome) / "loci_lastz.tsv", loci)
    if not a.parse_only:
        timing(a.genome, "lastz", elapsed, a.jobs, {"loci": len(loci)}, c_el)
    print(f"lastz: {len(rows)} HSPs, {len(loci)} loci, {elapsed:.1f} s")


def cigar_counts(cg: str) -> dict[str, int]:
    import re
    c: dict[str, int] = {}
    for n, op in re.findall(r"(\d+)([MIDFGNUV])", cg):
        c[op] = c.get(op, 0) + 1
    return c


PART_BP = 150_000_000


def split_genome(fa: pathlib.Path, d: pathlib.Path) -> list[pathlib.Path]:
    """Write the genome as FASTA parts of at most PART_BP (a longer sequence gets a part of its own)."""
    done = d / ".done"
    if done.exists():
        return sorted(d.glob("part*.fa"))
    d.mkdir(exist_ok=True)
    parts, size, fh = [], PART_BP + 1, None
    for name, seq in read_fasta(fa):
        if size + len(seq) > PART_BP:
            if fh:
                fh.close()
            parts.append(d / f"part{len(parts):03d}.fa"); fh = open(parts[-1], "wb"); size = 0
        fh.write(b">" + name.encode() + b"\n" + seq + b"\n"); size += len(seq)
    fh.close()
    done.touch()
    return parts


def cmd_miniprot(a):
    fa = genome_fa(a.genome)
    out = gdir(a.genome) / "miniprot.paf"
    # the whole-genome index does not fit in this sandbox's 15 GB, so the genome is split into parts of <= PART_BP
    # and searched part by part; the elapsed time covers all parts (the split itself is excluded)
    parts = split_genome(fa, gdir(a.genome) / "mp_parts")
    t0, c0 = time.time(), cpu()
    with open(out, "a" if a.parse_only else "w") as fh:
        for part in ([] if a.parse_only else parts):
            subprocess.run([str(MINIPROT), "-t", str(a.jobs), "-N", "200", "-p", "0.05", "--outs=0.05", "--outc=0.05",
                            "--outn=1000", str(part), str(REFS_AA)], stdout=fh, check=True)
    elapsed, c_el = time.time() - t0, cpu() - c0
    refs = load_refs(); names = list(refs)
    segs = []
    for line in open(out):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12 or f[0].startswith("#"):
            continue
        q, qlen, qs, qe, strand, chrom, _, ts, te, nmatch, alen = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], \
            f[6], int(f[7]), int(f[8]), int(f[9]), int(f[10])
        tags = dict((t[:2], t[5:]) for t in f[12:])
        cov_nt = 3 * (qe - qs)
        if cov_nt < MIN_COVER:
            continue
        cc = cigar_counts(tags.get("cg", ""))
        fs, st = int(tags.get("fs", 0)), int(tags.get("st", 0))
        introns = cc.get("N", 0) + cc.get("U", 0) + cc.get("V", 0)
        ident = 100.0 * nmatch / alen if alen else None  # PAF col 10/11: matching / aligned nucleotides
        detail = [f"fs={fs}", f"st={st}", f"introns={introns}"]
        segs.append((chrom, names.index(q), strand, ts, te, 3 * qs, 3 * qe, cov_nt, ident, 1, detail))
    loci = merge_refs(join_collinear(segs), names)
    classify(loci, None, refs, use_orf=False)
    for L in loci:
        broken = sum(int(kv.split("=")[1]) for kv in L["detail"].split(";")[-1].split(",") if "=" in kv) > 0
        L["status"] = "pseudogene" if (broken or not L["full"]) else "intact"
    write_loci(gdir(a.genome) / "loci_miniprot.tsv", loci)
    if not a.parse_only:
        timing(a.genome, "miniprot", elapsed, a.jobs, {"loci": len(loci)}, c_el)
    print(f"miniprot: {len(loci)} loci, {elapsed:.1f} s")


def seqid(x: str) -> str:
    """blast with -parse_seqids reports 'ref|NC_044374.1|'; the genome names the sequence 'NC_044374.1'. Pure."""
    return x.split("|")[1] if "|" in x else x


BLAST_FMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident"


def cmd_blast(a):
    fa = genome_fa(a.genome)
    d = gdir(a.genome)
    db = d / "blastdb" / "genome"
    if a.parse_only:
        blast_loci(a.genome); return
    t0 = time.time()
    if not pathlib.Path(str(db) + ".njs").exists() and not list(db.parent.glob("genome*.nsq")):
        db.parent.mkdir(exist_ok=True)
        subprocess.run(["makeblastdb", "-in", str(fa), "-dbtype", "nucl", "-out", str(db), "-parse_seqids"],
                       check=True, stdout=subprocess.DEVNULL)
    t_db = time.time() - t0
    t1 = time.time()
    subprocess.run(["blastn", "-task", "dc-megablast", "-query", str(REFS_NT), "-db", str(db), "-evalue", "1e-10",
                    "-max_target_seqs", "5000", "-max_hsps", "5000", "-num_threads", str(a.jobs),
                    "-outfmt", BLAST_FMT, "-out", str(d / "blastn.tsv")], check=True)
    t_bn = time.time() - t1
    t2 = time.time()
    subprocess.run(["tblastn", "-query", str(REFS_AA), "-db", str(db), "-evalue", "1e-10", "-max_target_seqs",
                    "5000", "-max_hsps", "5000", "-num_threads", str(a.jobs), "-outfmt", BLAST_FMT,
                    "-out", str(d / "tblastn.tsv")], check=True)
    t_tbn = time.time() - t2
    timing(a.genome, "blastn", t_bn + t_db, a.jobs, {"makeblastdb_s": round(t_db, 1)})
    timing(a.genome, "tblastn", t_tbn, a.jobs)
    blast_loci(a.genome)


def blast_segments(path, names, protein=False):
    segs = []
    for line in open(path):
        f = line.split("\t")
        q, s, pid, length = f[0], f[1], float(f[2]), int(f[3])
        qs, qe, ss, se, nid = int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[12])
        strand = "+" if ss <= se else "-"
        gs, ge = min(ss, se) - 1, max(ss, se)
        rs, re_ = (qs - 1, qe) if not protein else (3 * (qs - 1), 3 * qe)
        if re_ - rs < (60 if protein else 100):
            continue
        segs.append((seqid(s), names.index(q), strand, gs, ge, rs, re_, re_ - rs, pid, 1, []))
    return segs


def blast_loci(genome):
    """Per-method loci for blastn and tblastn (used for the truth set, and as an extra reference point)."""
    d = gdir(genome); refs = load_refs(); names = list(refs)
    g = load_genome(genome_fa(genome))
    for m, prot in (("blastn", False), ("tblastn", True)):
        loci = merge_refs(join_collinear(blast_segments(d / f"{m}.tsv", names, prot)), names)
        classify(loci, g, refs)
        write_loci(d / f"loci_{m}.tsv", loci)
        print(f"{m}: {len(loci)} loci")


# ---------------------------------------------------------------- truth set (cs10)
def dna_identity(genome_path: pathlib.Path, names):
    """Per (chrom, strand) interval list of blastn HSPs with identity, for identity lookup. Pure-ish."""
    rows = []
    for line in open(genome_path):
        f = line.split("\t")
        ss, se = int(f[8]), int(f[9])
        rows.append((seqid(f[1]), min(ss, se) - 1, max(ss, se), f[0], int(f[12]), int(f[3]), int(f[6]) - 1, int(f[7]),
                     "+" if ss <= se else "-"))
    return rows


def cmd_truth(a):
    """Truth loci for cs10: union of blastn/tblastn loci, RefSeq oxidocyclase/BBE-like genes and Grassa S9, merged by
    overlap; each gets DNA identity to every reference (blastn identities summed over its HSPs) and an ORF status."""
    d = gdir(a.genome); refs = load_refs(); names = list(refs)
    g = load_genome(genome_fa(a.genome))
    cand = []  # (chrom, start, end, strand or ".", source)
    for m in ("blastn", "tblastn"):
        for L in read_loci(d / f"loci_{m}.tsv"):
            if L["ref_cover"] >= MIN_COVER:
                cand.append((L["chrom"], L["start"], L["end"], L["strand"], m))
    for line in open(FAMILY_GENES):
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[8] == "oxidocyclase/BBE-like":
            cand.append((f[2], int(f[3]) - 1, int(f[4]), f[5], "refseq:" + f[1]))
    for i, line in enumerate(open(GRASSA_S9)):
        if i == 0:
            continue
        f = line.rstrip("\n").split("\t")
        cand.append(("NC_044378.1", int(f[0]) - 1, int(f[1]), ".", f"grassaS9:{f[4]}:{f[5]}"))
    cand.sort()
    merged = []
    for c in cand:
        if merged and merged[-1]["chrom"] == c[0] and c[1] < merged[-1]["end"]:
            M = merged[-1]; M["end"] = max(M["end"], c[2]); M["sources"].append(c[4])
            if M["strand"] == "." and c[3] != ".":
                M["strand"] = c[3]
        else:
            merged.append({"chrom": c[0], "start": c[1], "end": c[2], "strand": c[3], "sources": [c[4]]})
    hsps = dna_identity(d / "blastn.tsv", names)
    tbl = dna_identity(d / "tblastn.tsv", names)
    mp = {(L["chrom"], L["start"]): L for L in read_loci(d / "loci_miniprot.tsv")}
    mp_loci = list(mp.values())
    out = []
    for i, M in enumerate(merged):
        ident, cover, rspan, strands = {}, {}, {}, {}
        for (c, s, e, q, nid, ln, qs, qe, st) in hsps:
            if c == M["chrom"] and s < M["end"] and e > M["start"]:
                x = ident.setdefault(q, [0, 0]); x[0] += nid; x[1] += ln
                strands.setdefault(q, {}).setdefault(st, 0); strands[q][st] += qe - qs
                cover[q] = cover.get(q, 0) + (qe - qs)
                lo, hi = rspan.get(q, (10**9, -1)); rspan[q] = (min(lo, qs), max(hi, qe))
        aa = {}
        for (c, s, e, q, nid, ln, qs, qe, _) in tbl:
            if c == M["chrom"] and s < M["end"] and e > M["start"]:
                x = aa.setdefault(q, [0, 0]); x[0] += nid; x[1] += ln
        pid = {q: 100.0 * v[0] / v[1] for q, v in ident.items()}
        # the closest reference: most identical aligned bases (span x identity), so a 1,638-bp THCAS hit at 95.9%
        # does not beat a 1,635-bp CBCAS hit at 99.8% on length alone
        best = max(pid, key=lambda q: min(cover[q], len(refs[q])) * pid[q]) if pid else None
        aapid = {q: 100.0 * v[0] / v[1] for q, v in aa.items()}
        T = {"id": f"T{i + 1:03d}", "chrom": M["chrom"], "start": M["start"], "end": M["end"], "strand": M["strand"],
             "sources": ",".join(sorted(set(M["sources"]))),
             "best_ref": best, "dna_identity": round(pid[best], 2) if best else None,
             "dna_ref_span": (rspan[best][1] - rspan[best][0]) if best else 0,
             **{f"pid_{q}": round(pid[q], 2) for q in pid}, **{f"aa_{q}": round(aapid[q], 1) for q in aapid}}
        # status: ORF test on the implied gene region of the best DNA hit; partial copies are pseudogenes
        if best:
            rl = len(refs[best])
            lo, hi = rspan[best]
            T["full"] = (hi - lo) >= FULL_FRACTION * rl
            T["strand"] = max(strands[best], key=strands[best].get)
            L = {"chrom": M["chrom"], "start": M["start"], "end": M["end"], "strand": T["strand"],
                 "ref": best, "ref_start": lo, "ref_end": hi}
            classify([L], g, refs)
            T["orf_aa"] = L["orf_aa"]
            T["status"] = L["status"]
        else:
            T["full"] = None; T["status"] = "n/a (protein-level only)"
        T["bin"] = bin_of(T["dna_identity"])
        mm = [L for L in mp_loci if L["chrom"] == M["chrom"] and L["start"] < M["end"] and L["end"] > M["start"]]
        T["miniprot_detail"] = mm[0]["detail"] if mm else ""
        out.append(T)
    cols = ["id", "chrom", "start", "end", "strand", "best_ref", "dna_identity", "dna_ref_span", "bin", "full",
            "status", "orf_aa", "pid_THCAS", "pid_CBDAS", "pid_CBCAS", "aa_THCAS", "aa_CBDAS", "aa_CBCAS",
            "miniprot_detail", "sources"]
    with open(d / "truth.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for T in out:
            fh.write("\t".join("" if T.get(c) is None else str(T.get(c)) for c in cols) + "\n")
    print(f"truth: {len(out)} loci -> {d / 'truth.tsv'}")
    for b, *_ in BINS:
        print(b, sum(1 for T in out if T["bin"] == b))


def bin_of(pid):
    if pid is None:
        return "<80"
    for b, lo, hi in BINS:
        if lo <= pid < hi:
            return b
    return "<80"


def read_tsv(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh]


def overlaps(a, b, frac=0.5):
    """True if intervals overlap by at least frac of the shorter one. Pure."""
    o = min(a[1], b[1]) - max(a[0], b[0])
    return o > 0 and o >= frac * min(a[1] - a[0], b[1] - b[0])


METHODS = ["k11", "k15", "k21", "k31", "k15d1", "sp31w24", "lastz", "miniprot", "blastn", "tblastn"]


def cmd_score(a):
    d = gdir(a.genome)
    truth = read_tsv(d / "truth.tsv")
    for T in truth:
        T["start"], T["end"] = int(T["start"]), int(T["end"])
    tim = json.loads((d / "timings.json").read_text())
    res = {}
    for m in METHODS:
        p = d / f"loci_{m}.tsv"
        if not p.exists():
            continue
        pred = read_loci(p)
        matched_truth = {}
        fp = []
        orient_ok = orient_n = 0
        for L in pred:
            # one truth locus per predicted locus (the largest overlap), so a locus spliced across two copies
            # earns one, not two
            hits = [T for T in truth if T["chrom"] == L["chrom"] and overlaps((L["start"], L["end"]), (T["start"], T["end"]))]
            if not hits:
                fp.append(L); continue
            T = max(hits, key=lambda T: min(L["end"], T["end"]) - max(L["start"], T["start"]))
            matched_truth.setdefault(T["id"], []).append(L)
        r = {"loci": len(pred), "tp_loci": len(pred) - len(fp), "fp_loci": len(fp),
             "precision": round((len(pred) - len(fp)) / len(pred), 3) if pred else None,
             "fp": [f"{L['chrom']}:{L['start']}-{L['end']}{L['strand']} {L['ref']} cover={L['ref_cover']}" for L in fp],
             "time": tim.get(m), "bins": {}}
        status_agree = status_n = full_agree = 0
        for b, *_ in BINS:
            tb = [T for T in truth if T["bin"] == b]
            found = [T for T in tb if T["id"] in matched_truth]
            r["bins"][b] = {"truth": len(tb), "found": len(found),
                            "recall": round(len(found) / len(tb), 3) if tb else None,
                            "missed": [T["id"] for T in tb if T["id"] not in matched_truth]}
        for b, lo, hi in SUBBINS:
            tb = [T for T in truth if T["bin"] == "<80" and (
                (lo is None and not T["dna_identity"]) or (lo is not None and T["dna_identity"]
                                                         and lo <= float(T["dna_identity"]) < hi))]
            r["bins"]["<80:" + b] = {"truth": len(tb), "found": sum(1 for T in tb if T["id"] in matched_truth)}
        for T in truth:
            Ls = matched_truth.get(T["id"])
            if not Ls:
                continue
            L = max(Ls, key=lambda x: x["ref_cover"])
            if T["strand"] in "+-":
                orient_n += 1; orient_ok += L["strand"] == T["strand"]
            # intact/full agreement only where the truth status is meaningful: the synthase clade (>= 80% DNA
            # identity, full-length DNA alignments). Below 80% the truth locus comes from partial HSPs.
            if T["status"] in ("intact", "pseudogene") and T["bin"] != "<80":
                status_n += 1; status_agree += L.get("status") == T["status"]
                full_agree += (L["full"] == (T["full"] == "True"))
        r["orientation_agree"] = f"{orient_ok}/{orient_n}"
        r["status_agree"] = f"{status_agree}/{status_n}"
        r["full_agree"] = f"{full_agree}/{status_n}"
        res[m] = r
    (d / "score.json").write_text(json.dumps(res, indent=1))
    allb = [b for b, *_ in BINS] + ["<80:" + b for b, *_ in SUBBINS]
    print(f"{'method':9} {'loci':>5} {'prec':>6} " + " ".join(f"{b:>11}" for b in allb) + "  orient status full  wall_s cpu_s")
    for m, r in res.items():
        cells = " ".join(f"{r['bins'][b]['found']:>5}/{r['bins'][b]['truth']:<5}" for b in allb)
        t = f"{r['time']['wall_s']} {r['time'].get('cpu_s', '')}" if r["time"] else ""
        print(f"{m:9} {r['loci']:>5} {r['precision']!s:>6} {cells}  {r['orientation_agree']:>6} {r['status_agree']:>6} "
              f"{r['full_agree']:>6} {t}")


def cmd_census(a):
    """Per-genome copy counts. The synthase clade is a lastz locus at >= 80% DNA identity to its closest reference
    (cs10 has nothing between 78% and 90%); every other method's loci are split by overlap with those. Functional-
    looking copies: lastz >= 98% identity, full and intact (ORF test), by reference."""
    d = gdir(a.genome)
    lz = read_loci(d / "loci_lastz.tsv")
    clade = [L for L in lz if L["identity"] is not None and L["identity"] >= 80]
    in_clade = lambda L: any(C["chrom"] == L["chrom"] and overlaps((L["start"], L["end"]), (C["start"], C["end"]))
                             for C in clade)
    out = {"functional_lastz": {}, "methods": {}}
    for L in clade:
        if L["identity"] >= 98 and L["full"] and L["status"] == "intact":
            out["functional_lastz"][L["ref"]] = out["functional_lastz"].get(L["ref"], 0) + 1
    for m in METHODS:
        p = d / f"loci_{m}.tsv"
        if not p.exists():
            continue
        loci = read_loci(p)
        cl = [L for L in loci if in_clade(L)]
        out["methods"][m] = {
            "loci": len(loci), "clade": len(cl), "outside_clade": len(loci) - len(cl),
            "clade_full_intact": sum(1 for L in cl if L["full"] and L["status"] == "intact"),
            "clade_full_pseudogene": sum(1 for L in cl if L["full"] and L["status"] != "intact"),
            "clade_partial": sum(1 for L in cl if not L["full"])}
    (d / "census.json").write_text(json.dumps(out, indent=1))
    print(a.genome, json.dumps(out["functional_lastz"]))
    for m, r in out["methods"].items():
        print(f"  {m:9} " + " ".join(f"{k}={v}" for k, v in r.items()))


# ---------------------------------------------------------------- spike-in (fills the bins cs10 leaves empty)
SPIKE_LEVELS = [97, 92, 87, 82, 77, 72, 67]
SPIKE_REPS = 3
SPIKE_STEP = 250_000


def diverge(seq: bytes, identity: float, rng) -> bytes:
    """A copy at the nominal identity: (1 - identity) x length mutation events, 85% substitutions (61% of them
    transitions, the cs10 vs KOMP spectrum in results/ideas_eisenstein_2026-11-10.md) and 15% indels of 1-3 bp. Pure
    given rng."""
    s = bytearray(seq)
    n = int(round((1 - identity / 100.0) * len(s)))
    sites = sorted(rng.choice(len(s) - 4, n, replace=False), reverse=True)
    ts = {ord("A"): b"G", ord("G"): b"A", ord("C"): b"T", ord("T"): b"C"}
    tv = {ord("A"): b"CT", ord("G"): b"CT", ord("C"): b"AG", ord("T"): b"AG"}
    for i in sites:  # right to left, so earlier indices stay valid
        r = rng.random()
        if r < 0.85:
            s[i] = ts[s[i]][0] if rng.random() < 0.61 else tv[s[i]][int(rng.integers(0, 2))]
        else:
            L = int(rng.integers(1, 4))
            if r < 0.925:
                del s[i:i + L]
            else:
                s[i:i] = bytes(rng.choice(np.frombuffer(b"ACGT", np.uint8), L))
    return bytes(s)


def cmd_spike(a):
    """Plant copies of each reference at SPIKE_LEVELS identity, both strands, into cs10 NC_044370.1 0-20 Mb (no family
    member there: the nearest truth locus is at 29.6 Mb) -> DATA/spike/genome.fa and truth.tsv."""
    rng = np.random.default_rng(20261111)
    refs = load_refs()
    host = bytearray(fetch_seq(DATA / "cs10" / "genome.fa", "NC_044370.1")[:20_000_000].upper())
    rows, at = [], 100_000
    for lvl in SPIKE_LEVELS:
        for name, ref in refs.items():
            for rep_ in range(SPIKE_REPS):
                c = diverge(ref, lvl, rng)
                strand = "+" if rng.random() < 0.5 else "-"
                if strand == "-":
                    c = revcomp(c)
                host[at:at + len(c)] = c
                rows.append({"id": f"S{len(rows) + 1:03d}", "chrom": "spike", "start": at, "end": at + len(c),
                             "strand": strand, "best_ref": name, "dna_identity": lvl, "dna_ref_span": len(ref),
                             "bin": bin_of(lvl), "full": True, "status": "n/a (spike)", "sources": "spike"})
                at += SPIKE_STEP
    d = gdir("spike")
    (d / "genome.fa").write_bytes(b">spike\n" + bytes(host) + b"\n")
    for f in d.glob("genome.fa.fai"):
        f.unlink()
    seq_names(d / "genome.fa")
    cols = list(rows[0])
    with open(d / "truth.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"spike: {len(rows)} copies planted -> {d}")


# ---------------------------------------------------------------- references
def cmd_verify_refs(a):
    """Re-fetch every accession in refs/manifest.json from NCBI efetch and compare sequence and SHA-256."""
    import hashlib
    import re
    import urllib.request
    man = json.loads(MANIFEST.read_text())
    nt = dict(read_fasta(REFS_NT)); aa = dict(read_fasta(REFS_AA))
    ok = True
    for name, e in man["references"].items():
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&retmode=text"
               f"&tool=gene_copies_pilot&id={e['nucleotide']}")
        gb = urllib.request.urlopen(url, timeout=60).read().decode()
        seq = re.sub(r"[^acgt]", "", gb.split("ORIGIN")[1]).upper()
        a0, b0 = e["cds"]
        cds = seq[a0 - 1:b0].encode()
        same = cds == nt[name].upper()
        prot = translate(cds).rstrip("*")
        same_aa = prot == aa[name].decode()
        if e.get("protein"):
            url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&retmode=text"
                   f"&tool=gene_copies_pilot&id={e['protein']}")
            pf = urllib.request.urlopen(url, timeout=60).read().decode()
            same_aa = same_aa and "".join(pf.splitlines()[1:]) == prot
        sha = hashlib.sha256(cds).hexdigest()
        good = same and same_aa and sha == e["sha256_cds"]
        ok &= good
        print(f"{name}: {e['nucleotide']} CDS {a0}..{b0} {'OK' if good else 'MISMATCH'} sha256 {sha[:16]}")
        time.sleep(0.4)
    if not ok:
        raise SystemExit("⛔ a reference does not match NCBI")


# ---------------------------------------------------------------- self-test
def mutate(seq: bytes, rate: float, rng) -> bytes:
    s = bytearray(seq)
    for i in np.nonzero(rng.random(len(s)) < rate)[0]:
        s[i] = b"ACGT"[(b"ACGT".index(s[i]) + int(rng.integers(1, 4))) % 4]
    return bytes(s)


def cmd_self_test(a):
    rng = np.random.default_rng(7)
    rnd = lambda n: rng.choice(np.frombuffer(b"ACGT", np.uint8), n).tobytes()
    # strand symmetry of the canonical keys, contiguous and spaced
    s = rnd(500)
    for seed in SEEDS:
        span, mask, _ = SEEDS[seed]
        pos = seed_positions(span, mask)
        k1, f1, v1 = window_keys(LUT[np.frombuffer(s, np.uint8)], span, pos)
        k2, f2, v2 = window_keys(LUT[np.frombuffer(revcomp(s), np.uint8)], span, pos)
        assert np.array_equal(k1, k2[::-1]), seed
        assert v1.all() and v2.all()
    # rc_key of a reference key equals the genome rc key for a palindromic mask
    pos = seed_positions(31, PALINDROMIC_MASK)
    w = s[:31]
    kf, _, _ = window_keys(LUT[np.frombuffer(w, np.uint8)], 31, pos)
    fw = 0
    for p in pos:
        fw = (fw << 2) | int(LUT[w[p]])
    assert min(fw, int(rc_key(np.array([fw], np.uint64), len(pos))[0])) == int(kf[0])
    # N handling
    k, f, v = window_keys(LUT[np.frombuffer(b"ACGTNACGTACGTAC", np.uint8)], 4, seed_positions(4, None))
    assert list(v) == [True] + [False] * 4 + [True] * 7
    # planted copies: a 1,635-bp reference at 100%, 92%, 85% (+), and 97% (-) in random sequence
    ref = rnd(1635)
    refs = {"R": ref}
    plant = [(5000, 0.0, "+"), (20000, 0.08, "+"), (40000, 0.15, "+"), (60000, 0.03, "-")]
    g = bytearray(rnd(80000))
    for at, rate, strand in plant:
        c = mutate(ref, rate, rng)
        if strand == "-":
            c = revcomp(c)
        g[at:at + len(c)] = c
    g = bytes(g)
    with tempfile.TemporaryDirectory() as t:
        fa = pathlib.Path(t) / "g.fa"
        fa.write_text(">chrT\n" + g.decode() + "\n")
        seq_names(fa)
        for seed, expect in (("k11", 4), ("k15", 4), ("k15d1", 4), ("k31", 2)):
            idx = ref_index(refs, seed)
            _, (gp, rid, rpos, same) = scan_sequence((fa, "chrT", seed, idx))
            segs = [s for s in segments_from_hits(gp, rid, rpos, same, SEEDS[seed][0]) if s[6] >= MIN_COVER]
            chains = join_collinear([("chrT",) + s[:7] + (None, 1, []) for s in segs])
            loci = merge_refs([c[:1] + c[1:] for c in chains], ["R"])
            found = {(L["start"] // 1000, L["strand"]) for L in loci}
            want = {(at // 1000, st) for at, rate, st in plant[:expect]} if seed != "k31" else None
            if want is not None:
                assert want <= found, (seed, found, want, [(L["start"], L["end"], L["ref_cover"]) for L in loci], [s[2:8] for s in segments_from_hits(gp, rid, rpos, same, SEEDS[seed][0])][:40])
            assert all(L["strand"] == dict((at // 1000, st) for at, _, st in plant)[L["start"] // 1000]
                       for L in loci if L["start"] // 1000 in {at // 1000 for at, _, _ in plant}), seed
            classify(loci, {"chrT": g}, refs)
    # ORF test: a frameshift breaks it
    thcas = load_refs()["THCAS"]
    assert longest_orf(thcas) == 545 and longest_orf(revcomp(thcas)) < 545
    assert longest_orf(thcas[:800] + thcas[801:]) < ORF_FRACTION * 544
    # collinear join: an insertion splitting a copy is rejoined
    segs = [("c", 0, "+", 1000, 1800, 0, 800, 800, 99.0, 1, []), ("c", 0, "+", 6800, 7635, 800, 1635, 835, 98.0, 1, [])]
    j = join_collinear(segs)
    assert len(j) == 1 and j[0][5] == 0 and j[0][6] == 1635
    print("self-test OK")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    for c in ("stage", "kmer", "lastz", "miniprot", "blast", "truth", "score", "census"):
        p = sp.add_parser(c)
        p.add_argument("--genome", required=True)
        p.add_argument("--jobs", type=int, default=8)
        if c == "stage":
            p.add_argument("--src", required=True)
        if c in ("blast", "miniprot", "lastz"):
            p.add_argument("--parse-only", action="store_true", help="rebuild loci from the existing raw output")
        if c == "kmer":
            p.add_argument("--seeds", nargs="*")
    sp.add_parser("verify-refs"); sp.add_parser("self-test"); sp.add_parser("spike")
    a = ap.parse_args()
    {"stage": cmd_stage, "kmer": cmd_kmer, "lastz": cmd_lastz, "miniprot": cmd_miniprot, "blast": cmd_blast,
     "truth": cmd_truth, "score": cmd_score, "census": cmd_census, "verify-refs": cmd_verify_refs, "spike": cmd_spike,
     "self-test": cmd_self_test}[a.cmd](a)


if __name__ == "__main__":
    main()
