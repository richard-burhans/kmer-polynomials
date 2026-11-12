#!/usr/bin/env python3
"""Turn one run_paper_bench.sh laila job into the paper's speed numbers.

    python3 manuscript/paper_bench_table.py <laila_perf results dir>   markdown on stdout

Reads run.log (the command of each run_N) and counters.tsv (its counters). For each binary, the median of its
__load_only__ runs is subtracted from each method run, and the result divided by the input k-mer count. A method's
value is the MEDIAN over its repeats, its spread (max - min) / median, and its minimum is shown too. Each method is
compared with batched ntHash timed by the same binary (bench, which has none, uses bench_e's). The per-kept cost and
break-even rates come from each batched sampler's slope between 1/16 and 1/4. It refuses to print if any run exited
non-zero, if bench_kept's identity check did not pass, or if repeats of a bench_kept method kept different k-mers.
"""
import collections
import pathlib
import re
import statistics
import sys

KMERS = 736_547_799  # 31-mers of cs10
RUN = re.compile(r'^=== (run_\d+): (.*)$')
WHAT = re.compile(r'\./(bench_e|bench_kept|bench|bench_spaced_e) .*?(?:"only=([^"]*)"|( check)|(\d+ \d+))?$')
NT = {"bench_e": "ntHash batched (ntbatch)", "bench_kept": "ntHash batched", "bench": "ntHash batched (ntbatch)"}
NT_BIN = {"bench_e": "bench_e", "bench_kept": "bench_kept", "bench": "bench_e"}
SLOPES = [("bench_e", "e batched"), ("bench_kept", "lib-e 8 bat"), ("bench_kept", "nt 8 bat"), ("bench_kept", "uh 0 bat")]


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
    meta = dict(ln.split("\t", 1) for ln in (d / "meta.tsv").read_text().splitlines() if "\t" in ln)
    return cmds, logs, ctr, meta


def compute(d: pathlib.Path):
    """Validate one job's raw outputs and return (load baselines, per-method rows, spaced-seed run, run logs).

    rows maps (binary, method) to (median cycles/k-mer, spread, min, instructions/k-mer, branch-miss rate, n).
    """
    cmds, logs, ctr, meta = load(d)
    bad = {k: v for k, v in meta.items() if (k.startswith("pcount_exit_") or k == "build_exit") and v.strip() != "0"}
    if bad:
        sys.exit(f"non-zero exits: {bad}")
    per = collections.defaultdict(list)  # (binary, method) -> [run]
    spaced = check = None
    kept = collections.defaultdict(set)  # bench_kept method -> {(kept, xor)}
    for run, cmd in cmds.items():
        m = WHAT.search(cmd)
        if not m:
            sys.exit(f"cannot parse {run}: {cmd}")
        if m.group(1) == "bench_spaced_e":
            spaced = run
        elif m.group(3):
            check = run
        else:
            per[(m.group(1), m.group(2))].append(run)
            for ln in logs[run]:
                r = re.match(r"RESULT (.*)\tkept=(\d+)\txor=([0-9a-f]+)", ln)
                if r:
                    kept[r.group(1)].add((int(r.group(2)), r.group(3)))
    if check is None or not any("every evaluator identical" in ln for ln in logs[check]):
        sys.exit("bench_kept check missing or failed")
    for meth, s in kept.items():
        if len(s) != 1:
            sys.exit(f"{meth}: repeats differ {s}")
    by_rate = collections.defaultdict(set)
    for meth, s in kept.items():
        if meth not in ("ntHash batched", "uh all", "__load_only__"):
            by_rate[meth.split()[-1]].add(next(iter(s))[0])
    for rate, s in by_rate.items():
        if len(s) != 1:
            sys.exit(f"samplers at {rate} keep different counts {s}")
    base = {}
    for b in sorted({b for b, meth in per if meth == "__load_only__"}):
        cyc = [ctr[r]["cycles"] for r in per[(b, "__load_only__")]]
        base[b] = (statistics.median(cyc), statistics.median(ctr[r]["instructions"] for r in per[(b, "__load_only__")]),
                   len(cyc), (max(cyc) - min(cyc)) / statistics.median(cyc))
    rows = {}
    for (b, meth), runs in per.items():
        if meth == "__load_only__":
            continue
        cyc = [(ctr[r]["cycles"] - base[b][0]) / KMERS for r in runs]
        ins = [(ctr[r]["instructions"] - base[b][1]) / KMERS for r in runs]
        bm = statistics.median(ctr[r]["branch-misses"] / ctr[r]["branches"] for r in runs)
        med = statistics.median(cyc)
        rows[(b, meth)] = (med, (max(cyc) - min(cyc)) / med, min(cyc), statistics.median(ins), bm, len(runs))
    return base, rows, spaced, logs


def slope_fit(rows, b: str, fam: str):
    """Cycles per kept k-mer and intercept of a batched sampler, from its 1/16 and 1/4 medians."""
    lo, hi = rows[(b, f"{fam} 1/16")][0], rows[(b, f"{fam} 1/4")][0]
    slope = (hi - lo) / (0.25 - 0.0625)
    return slope, lo - slope * 0.0625


def main() -> int:
    base, rows, spaced, logs = compute(pathlib.Path(sys.argv[1]))
    print("Load baselines (median of each binary's __load_only__ runs): " + "; ".join(
        f"{b} {c / 1e9:.2f}e9 cycles (n = {n}, spread {sp:.1%})" for b, (c, _, n, sp) in base.items()) + "\n")
    print("| method | binary | cycles/k-mer (median) | spread | min | instructions/k-mer | branch-miss rate | n "
          "| vs batched ntHash (same binary) |")
    print("|---|---|---|---|---|---|---|---|---|")
    for (b, meth), (c, sp, mn, i, bm, n) in rows.items():
        ref = rows[(NT_BIN[b], NT[b])][0]
        print(f"| {meth} | {b} | {c:.2f} | {sp:.1%} | {mn:.2f} | {i:.1f} | {bm:.2%} | {n} | {ref / c:.2f}× |")
    print()
    for b, fam in SLOPES:
        slope, b0 = slope_fit(rows, b, fam)
        ref = rows[(b, NT[b])][0]
        be = f"1/{slope / (ref - b0):.1f}" if ref > b0 else "never"
        print(f"- {b} `{fam}`: {slope:.1f} cycles per kept k-mer, intercept {b0:.2f}; beats batched ntHash "
              f"({ref:.2f}) below {be}")
    if spaced:
        print(f"\nSpaced seeds ({spaced}, its own wall-clock medians):\n```")
        print("\n".join(ln for ln in logs[spaced] if ln.strip()))
        print("```")
    return 0


if __name__ == "__main__":
    sys.exit(main())
