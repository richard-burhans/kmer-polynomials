#!/usr/bin/env python3
"""Turn one hash/run_kept_laila.sh job into cycles per input k-mer, per-kept costs and speed-ups.

    python3 hash/kept_table.py <laila_perf results dir> [--min]      markdown on stdout

--min uses each method's MINIMUM over repeats (and the minimum load run) instead of the median (and the mean load):
interference from other jobs on laila only ADDS cycles, so on a busy host the minimum is the cross-check.

Reads run.log (each run_N's command and its RESULT line) and counters.tsv. The mean of the job's __load_only__ runs
is subtracted from every method run, and the result divided by the input k-mer count. A method's value is the MEDIAN
over its repeats; spread is (max - min) / median. Checks, before printing anything:
  - every run of a method reports the same kept count and XOR (determinism);
  - every sampler at the same rate keeps the same number of k-mers (they share the gate);
  - samplers of the same hash family at the same rate give the same XOR;
  - the check run ended with "every evaluator identical".
Per kept k-mer: the slope between rates 1/16 and 1/4 (as manuscript/paper_bench_table.py), and, where the method's
rate-0 base was run, (cycles(rate) - cycles(0)) / rate at each rate. Break-even against batched ntHash: the rate r with
intercept + slope * r = ntHash, from the 1/16-1/4 line.
"""
import collections
import pathlib
import re
import statistics
import sys

KMERS = 736_547_799  # 31-mers of cs10 (results/paper_bench_2026-11-11.md)
RATES = {"0": 0.0, "1/1000": 1 / 1000, "1/100": 1 / 100, "1/16": 1 / 16, "1/4": 1 / 4}
RUN = re.compile(r"^=== (run_\d+): (.*)$")
ONLY = re.compile(r'only=([^"]*)"')
RESULT = re.compile(r"^RESULT (.*)\tkept=(\d+)\txor=([0-9a-f]+)\tsecs=([0-9.]+)$")
FAMILY = {"lib-e": "e", "ec": "e", "ep": "e", "gp": "g", "lib-g": "g", "mp": "m", "mc": "m", "lib-m": "m",
          "nt": "nt", "uh": "uh"}


def load(d: pathlib.Path):
    cmds, results, cur, check_ok = {}, {}, None, False
    for line in (d / "run.log").read_text().splitlines():
        m = RUN.match(line)
        if m:
            cur = m.group(1)
            cmds[cur] = m.group(2)
            continue
        m = RESULT.match(line)
        if m and cur:
            results[cur] = (m.group(1), int(m.group(2)), m.group(3), float(m.group(4)))
        if "check: every evaluator identical to its reference" in line:
            check_ok = True
    ctr = collections.defaultdict(dict)
    for line in (d / "counters.tsv").read_text().splitlines()[1:]:
        run, ev, val = line.split("\t")[:3]
        ctr[run][ev] = float(val)
    return cmds, results, ctr, check_ok


def parse(method: str):
    """(ev, B, collect, rate, batch) or None for a non-sampler."""
    f = method.split()
    if len(f) >= 4 and f[0] in FAMILY and f[3] in RATES:
        return f[0], int(f[1]), f[2], f[3], (f[4] if len(f) > 4 else "b512")
    return None


def main() -> int:
    d = pathlib.Path(sys.argv[1])
    use_min = "--min" in sys.argv[2:]
    cmds, results, ctr, check_ok = load(d)
    errors = []
    if not check_ok:
        errors.append("the check run did not report every evaluator identical")
    by = collections.defaultdict(list)
    for run, cmd in cmds.items():
        m = ONLY.search(cmd)
        if m:
            by[m.group(1)].append(run)
    load_runs = by.pop("__load_only__")
    agg_load = min if use_min else statistics.mean
    base = {ev: agg_load(ctr[r][ev] for r in load_runs) for ev in ("cycles", "instructions", "branch-misses")}
    load_spread = (max(ctr[r]["cycles"] for r in load_runs) - min(ctr[r]["cycles"] for r in load_runs)) / base["cycles"]
    rows = {}
    for meth, runs in by.items():
        outs = {(results[r][1], results[r][2]) for r in runs if r in results}
        if len(outs) != 1 or len(runs) != sum(r in results for r in runs):
            errors.append(f"{meth}: runs disagree or lack a RESULT: {outs}")
        cyc = [(ctr[r]["cycles"] - base["cycles"]) / KMERS for r in runs]
        ins = [(ctr[r]["instructions"] - base["instructions"]) / KMERS for r in runs]
        bm = [(ctr[r]["branch-misses"] - base["branch-misses"]) / KMERS for r in runs]
        med = min(cyc) if use_min else statistics.median(cyc)
        rows[meth] = dict(cyc=med, spread=(max(cyc) - min(cyc)) / med, n=len(runs), ins=statistics.median(ins),
                          bm=statistics.median(bm), kept=results[runs[0]][1], xor=results[runs[0]][2])
    # shared gate and families
    kept_at, xor_at = collections.defaultdict(set), collections.defaultdict(set)
    for meth, row in rows.items():
        p = parse(meth)
        if p:
            kept_at[p[3]].add(row["kept"])
            xor_at[(FAMILY[p[0]], p[3])].add(row["xor"])
    for rate, ks in kept_at.items():
        if len(ks) != 1:
            errors.append(f"rate {rate}: samplers keep different counts {ks}")
    for key, xs in xor_at.items():
        if len(xs) != 1:
            errors.append(f"family {key}: different XORs {xs}")
    if errors:
        print("\n".join("ERROR: " + e for e in errors))
        return 1

    nt = rows["ntHash batched"]["cyc"]
    print(f"Statistic: {'MINIMUM' if use_min else 'median'} over repeats. Load baseline: {len(load_runs)} runs, {'min' if use_min else 'mean'} {base['cycles']:.4g} cycles, spread {load_spread:.1%}. "
          f"Batched ntHash: {nt:.2f} cycles per k-mer.\n")
    print("| method | kept | cycles/k-mer | spread | n | instructions/k-mer | branch misses/k-mer | vs batched ntHash |")
    print("|---|---|---|---|---|---|---|---|")
    for meth in sorted(rows, key=lambda m: (parse(m) is not None, parse(m)[0:3] if parse(m) else (m,),
                                            RATES[parse(m)[3]] if parse(m) else 0, parse(m)[4] if parse(m) else "")):
        r = rows[meth]
        print(f"| {meth} | {r['kept']:,} | {r['cyc']:.2f} | {r['spread']:.1%} | {r['n']} | {r['ins']:.1f} | "
              f"{r['bm']:.3f} | {nt / r['cyc']:.2f}x |")
    # per-kept costs
    print("\n| evaluator | collect | batch | per kept, slope 1/16 to 1/4 | intercept | break-even vs ntHash | "
          "per kept at 1/1000 | at 1/100 | (base: rate-0 run of) |")
    print("|---|---|---|---|---|---|---|---|---|")
    groups = collections.defaultdict(dict)
    for meth, r in rows.items():
        p = parse(meth)
        if p:
            groups[(p[0], p[1], p[2], p[4])][p[3]] = r["cyc"]
    for key in sorted(groups):
        g = groups[key]
        ev, B, col, batch = key
        slope = icpt = be = None
        if "1/16" in g and "1/4" in g:
            slope = (g["1/4"] - g["1/16"]) / (1 / 4 - 1 / 16)
            icpt = g["1/16"] - slope / 16
            be = (nt - icpt) / slope if slope > 0 else None
        # base: this method's own rate-0 run, else the shared batched base of the same collection
        b0, src = g.get("0"), f"{ev} {B} {col}"
        if b0 is None:
            for alt in (("lib-e", 8, col, "b512"), ("ec", 8, col, "b512")):
                if alt in groups and "0" in groups[alt] and col == "bat":
                    b0, src = groups[alt]["0"], f"{alt[0]} {alt[1]} {col}"
                    break
        sp = {r: (g[r] - b0) / RATES[r] for r in ("1/1000", "1/100") if r in g and b0 is not None}
        fmt = lambda v, f="{:.1f}": "—" if v is None else f.format(v)
        print(f"| {ev} {B} | {col} | {batch} | {fmt(slope)} | {fmt(icpt, '{:.2f}')} | "
              f"{'—' if be is None else ('1/%.1f' % (1 / be) if be > 0 else 'never')} | {fmt(sp.get('1/1000'))} | "
              f"{fmt(sp.get('1/100'))} | {src if b0 is not None else '—'} |")
    return 0


if __name__ == "__main__":
    sys.exit(main())
