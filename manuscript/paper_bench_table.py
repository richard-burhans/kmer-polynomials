#!/usr/bin/env python3
"""Turn one run_paper_bench.sh laila job into the paper's speed numbers.

    python3 manuscript/paper_bench_table.py <laila_perf results dir>   markdown on stdout

Reads run.log (the command of each run_N) and counters.tsv (its counters). For each binary, the mean of its
__load_only__ runs is subtracted from each method run, and the result divided by the input k-mer count. A method's
value is the mean over its repeats; the spread column is (max - min) / mean over those repeats. The per-kept cost and
break-even rates come from the batched samplers' slope between 1/16 and 1/4, as in results/kp7_eisenstein_2026-11-10.md.
"""
import collections
import pathlib
import re
import statistics
import sys

KMERS = 736_547_799  # 31-mers of cs10 (bench_e and bench print it; checked below)
RUN = re.compile(r'^=== (run_\d+): (.*)$')
WHAT = re.compile(r'\./(bench_e|bench|bench_spaced_e) .*?(?:"only=([^"]*)")?$')


def load(d: pathlib.Path):
    cmds, logs, cur = {}, collections.defaultdict(list), None
    for line in (d / "run.log").read_text().splitlines():
        m = RUN.match(line)
        if m:
            cur = m.group(1)
            cmds[cur] = m.group(2)
        elif cur:
            logs[cur].append(line)
    ctr = collections.defaultdict(dict)
    for line in (d / "counters.tsv").read_text().splitlines()[1:]:
        run, ev, val = line.split("\t")[:3]
        ctr[run][ev] = float(val)
    return cmds, logs, ctr


def main() -> int:
    d = pathlib.Path(sys.argv[1])
    cmds, logs, ctr = load(d)
    for run, lines in logs.items():  # every counted bench run must have hashed the expected k-mers
        for ln in lines:
            m = re.search(r"k-mers (\d+)\s+median", ln)
            if m and int(m.group(1)) != KMERS:
                sys.exit(f"{run}: {m.group(1)} k-mers, expected {KMERS}")
    per = collections.defaultdict(list)  # (binary, method) -> [run]
    spaced = None
    for run, cmd in cmds.items():
        m = WHAT.search(cmd)
        if not m:
            sys.exit(f"cannot parse {run}: {cmd}")
        if m.group(1) == "bench_spaced_e":
            spaced = run
        else:
            per[(m.group(1), m.group(2))].append(run)
    base = {b: {ev: statistics.mean(ctr[r][ev] for r in per[(b, "__load_only__")]) for ev in ("cycles", "instructions")}
            for b in sorted({b for b, meth in per if meth == "__load_only__"})}
    rows = {}
    for (b, meth), runs in per.items():
        if meth == "__load_only__":
            continue
        cyc = [(ctr[r]["cycles"] - base[b]["cycles"]) / KMERS for r in runs]
        ins = [(ctr[r]["instructions"] - base[b]["instructions"]) / KMERS for r in runs]
        bm = statistics.mean(ctr[r]["branch-misses"] / ctr[r]["branches"] for r in runs)
        rows[meth] = (statistics.mean(cyc), (max(cyc) - min(cyc)) / statistics.mean(cyc), statistics.mean(ins), bm,
                      len(runs), b)
    ref = rows["ntHash batched (ntbatch)"][0]
    print(f"Load baselines (mean of each binary's __load_only__ runs): " + "; ".join(
        f"{b} {v['cycles'] / 1e9:.2f}e9 cycles, {v['instructions'] / 1e9:.2f}e9 instructions "
        f"(n = {len(per[(b, '__load_only__')])}, spread "
        f"{(max(ctr[r]['cycles'] for r in per[(b, '__load_only__')]) - min(ctr[r]['cycles'] for r in per[(b, '__load_only__')])) / v['cycles']:.1%})"
        for b, v in base.items()) + "\n")
    print("| method | binary | cycles/k-mer | spread | instructions/k-mer | branch-miss rate | vs batched ntHash |")
    print("|---|---|---|---|---|---|---|")
    for meth, (c, sp, i, bm, n, b) in rows.items():
        print(f"| {meth} | {b} | {c:.2f} | {sp:.1%} | {i:.1f} | {bm:.2%} | {ref / c:.2f}× |")
    for fam in ("e", "m31"):
        lo, hi = rows[f"{fam} batched 1/16"][0], rows[f"{fam} batched 1/4"][0]
        slope = (hi - lo) / (0.25 - 0.0625)
        b0 = lo - slope * 0.0625
        print(f"\n{fam} batched: {slope:.1f} cycles per kept k-mer, base {b0:.2f}; beats batched ntHash ({ref:.2f}) "
              f"below 1/{slope / (ref - b0):.1f}; beats e sym_seq L2 ({rows['e sym_seq L2'][0]:.2f}) below "
              f"1/{slope / (rows['e sym_seq L2'][0] - b0):.1f}")
    if spaced:
        print(f"\nSpaced seeds ({spaced}, its own wall-clock rates):\n```")
        print("\n".join(ln for ln in logs[spaced] if ln.strip()))
        print("```")
    return 0


if __name__ == "__main__":
    sys.exit(main())
