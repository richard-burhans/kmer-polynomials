"""Summarize a laila run of recon/kp10_vs_iblt_job.sh (fetched by scripts/laila_perf.py) into TSV tables.

    python3 recon/kp10_vs_iblt_summarize.py <laila_perf result dir> <out dir>

Reads run.log (the "@@ <kind> <json>" lines and the pcount run headers) and counters.tsv (pcount's wall_s and
task-clock per run). Writes one TSV per table into <out dir>; every row carries `src`, the run.log line(s) it comes
from, so each number in the write-up can cite file:line. Prints the tables as Markdown.

Success criteria (decided before the run, recon/kp10_vs_iblt_job.sh):
  KP-10   `exact`: the recovered k-mers (canonicalized only for scoring) equal the truth, nothing unresolved.
  IBLT    `cancelled_exact`: km-peeler's listing equals the truth after k-mers listed on BOTH sides are dropped (the
          generous reading); `exact` (raw listing) is reported too. An IBLT trial is "detected" as failed when
          `ibltseq list` exits non-zero; a trial that exits 0 yet is not cancelled-exact is an UNDETECTED failure.
Failure-rate upper bounds are exact (Clopper-Pearson) one-sided 95%: 1 - 0.05^(1/n) when no failure is seen.
"""

from __future__ import annotations

import json
import pathlib
import statistics
import sys
from collections import defaultdict


def cp_upper(fail: int, n: int, alpha: float = 0.05) -> float:
    """One-sided upper (1 - alpha) Clopper-Pearson bound on a failure probability, by bisection on the binomial CDF."""
    if n == 0:
        return 1.0
    if fail == n:
        return 1.0
    from math import comb

    def cdf(p: float) -> float:
        return sum(comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(fail + 1))

    lo, hi = fail / n, 1.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if cdf(mid) > alpha:
            lo = mid
        else:
            hi = mid
    return hi


def load(run_dir: pathlib.Path):
    recs = defaultdict(list)  # kind -> [(line_no, obj)]
    for no, line in enumerate((run_dir / "run.log").read_text().splitlines(), 1):
        if line.startswith("@@ "):
            _, kind, js = line.split(" ", 2)
            recs[kind].append((no, json.loads(js)))
    runs = {}
    for no, line in enumerate((run_dir / "run.log").read_text().splitlines(), 1):
        if line.startswith("=== run_"):
            name, cmd = line[4:].split(": ", 1)
            runs[name] = {"cmd": cmd.strip(), "line": no}
    counters = defaultdict(dict)
    ctsv = run_dir / "counters.tsv"
    if ctsv.exists():
        for no, line in enumerate(ctsv.read_text().splitlines(), 1):
            f = line.split("\t")
            if no == 1 or len(f) < 3:
                continue
            counters[f[0]][f[1]] = (f[5] if len(f) > 5 and f[5] else f[2], no)
    return recs, runs, counters


def write(out: pathlib.Path, name: str, header: list[str], rows: list[list]) -> None:
    with open(out / f"{name}.tsv", "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"\n### {name}  ({out.name}/{name}.tsv)\n")
    print("| " + " | ".join(header) + " |")
    print("|" + "---|" * len(header))
    for r in rows:
        print("| " + " | ".join(str(x) for x in r) + " |")


def lines(xs) -> str:
    xs = sorted(xs)
    return f"run.log:{xs[0]}" if len(xs) == 1 else f"run.log:{xs[0]}-{xs[-1]} ({len(xs)} lines)"


def main() -> int:
    run_dir, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    recs, runs, counters = load(run_dir)

    # 1. sets and truth
    rows = []
    for no, d in recs["kset"]:
        o = d["out"]
        rows.append([d["genome"], o["records"], o["bases"], o["kmers"], o["distinct_canonical"], d["wall"], d["user"],
                     f"run.log:{no}"])
    write(out, "sets", ["genome", "records", "bases", "kmer_windows", "distinct_canonical_31mers", "kset_wall_s",
                        "kset_user_s", "src"], rows)
    rows = []
    diffs = {d["pair"]: (no, d["out"]) for no, d in recs["diff"]}
    for no, d in recs["kmc"]:
        p, s = d["pair"], d["summary"]
        dn, dd = diffs[p]
        rows.append([p, dd["distinct_A"], dd["distinct_B"], dd["A_minus_B"], dd["B_minus_A"], dd["D"],
                     f"{dd['D'] / dd['distinct_A']:.3e}", s["distinct_A"], s["distinct_B"], s["A_minus_B"],
                     s["B_minus_A"], d["same_amb"]["identical"] and d["same_bma"]["identical"], s["tool"],
                     f"run.log:{dn},{no}"])
    write(out, "truth", ["pair", "kmtools_A", "kmtools_B", "A_minus_B", "B_minus_A", "D", "D_over_A", "kmc_A", "kmc_B",
                         "kmc_A_minus_B", "kmc_B_minus_A", "kmc_difference_identical_element_for_element", "kmc",
                         "src"], rows)
    D = {p: dd["D"] for p, (_, dd) in diffs.items()}

    # 2. timed runs (one thread each) with pcount's wall and task-clock for the whole run
    def pc(cmd_suffix: str):
        for name, r in runs.items():
            if r["cmd"].endswith(cmd_suffix):
                c = counters.get(name, {})
                w, wl = c.get("wall_s", ("", 0))
                t, _ = c.get("task-clock", ("", 0))
                return name, w, f"{int(t) / 1e9:.1f}" if t else "", wl
        return "", "", "", 0

    rows = []
    for no, d in recs["kp10_timed"]:
        p = d["label"].split("_", 1)[1]
        name, w, t, cl = pc(f"kp10_timed {p}")
        cpu_sketch = d["cpu_hash_A"] + d["cpu_dedup_A"] + d["cpu_sketch_A"]
        rows.append([p, "KP-10 (kp10c, 1 thread)", d["sketch_bytes"], f"{d['bytes_per_true_difference']:.2f}", d["rounds"],
                     d["exact"], f"{d['cpu_hash_A']:.1f}+{d['cpu_dedup_A']:.1f}+{d['cpu_sketch_A']:.1f}={cpu_sketch:.1f}",
                     f"{d['cpu_decode']:.1f}", f"{d['cpu_recover_kmers']:.1f}", w, t,
                     f"run.log:{no}; counters.tsv:{cl} ({name})"])
    for no, d in recs["iblt_timed"]:
        p = d["pair"]
        name, w, t, cl = pc(f"iblt_timed {p}")
        a = d["a"]
        rows.append([p, f"km-peeler IBLT, n = 1.05 D, seed 1 (1 thread)", d["iblt_bytes"],
                     f"{d['iblt_bytes'] / d['D']:.2f}", 1, d["score"]["cancelled_exact"],
                     f"{a['user'] + a['sys']:.1f} (text alone {d['text']['user'] + d['text']['sys']:.1f})",
                     f"{d['diff']['user'] + d['list']['user'] + d['diff']['sys'] + d['list']['sys']:.1f}", "-", w, t,
                     f"run.log:{no}; counters.tsv:{cl} ({name})"])
    write(out, "timed", ["pair", "method", "bytes_one_side", "bytes_per_difference", "rounds", "exact",
                         "cpu_s_sketch_one_side (KP-10: hash+dedup+power sums; IBLT: text|build)",
                         "cpu_s_decode (KP-10: all rounds, both sides; IBLT: diff+list)", "cpu_s_back_to_kmers",
                         "pcount_wall_s_whole_run", "pcount_cpu_s_whole_run", "src"], rows)

    # 3. KP-10 over OS-random seeds at the true D
    rows = []
    for p in ("giga", "draft"):
        rs = [(no, d) for no, d in recs["kp10_seed"] if d["label"].endswith("_" + p)]
        if not rs:
            continue
        ds = [d for _, d in rs]
        n = len(ds)
        ex = sum(d["exact"] for d in ds)
        one = sum(d["rounds"] == 1 for d in ds)
        by = [d["sketch_bytes"] for d in ds]
        rows.append([p, n, ex, f"{cp_upper(n - ex, n):.3f}", one, f"{min(by)}-{max(by)}",
                     f"{statistics.mean(by) / D[p]:.2f}", ds[0]["bytes_round1"], f"{ds[0]['bytes_round1'] / D[p]:.2f}",
                     max(d["rounds"] for d in ds), sum(d["splits"] for d in ds) / n,
                     sum(d["missing"] for d in ds), sum(d["spurious"] for d in ds), lines(no for no, _ in rs)])
    write(out, "kp10_seeds", ["pair", "runs", "exact", "failure_upper95", "runs_finished_in_1_round",
                              "bytes_range", "mean_bytes_per_difference", "bytes_round1", "round1_bytes_per_difference",
                              "max_rounds", "mean_splits", "missing_total", "spurious_total", "src"], rows)

    # 4. KP-10 against the D estimate (fixed seeds)
    rows = []
    for no, d in recs["kp10_dest"]:
        p = d["label"].rsplit("_", 1)[1]
        rows.append([p, f"{d['dest'] / D[p]:g}", int(d["dest"]), d["bits0"], d["rounds"], d["splits"],
                     d["unresolved_buckets"], d["sketch_bytes"], f"{d['sketch_bytes'] / D[p]:.2f}", d["exact"],
                     f"{d['wall_reconcile']:.1f}", f"run.log:{no}"])
    write(out, "kp10_dest", ["pair", "dest_over_D", "dest", "bits0", "rounds", "splits", "unresolved", "bytes",
                             "bytes_per_difference", "exact", "wall_s_reconcile (8 threads)", "src"], rows)

    # 5. IBLT reliability by ratio
    rows = []
    groups = defaultdict(list)
    for no, d in recs["iblt_trial"]:
        groups[(d["pair"], d["ratio"])].append((no, d))
    for (p, r), g in sorted(groups.items()):
        ds = [d for _, d in g]
        n = len(ds)
        det = sum(d["list_exit"] == 0 for d in ds)
        raw = sum(d["score"]["exact"] for d in ds)
        canc = sum(d["score"]["cancelled_exact"] for d in ds)
        undetected = sum(d["list_exit"] == 0 and not d["score"]["cancelled_exact"] for d in ds)
        both = sum(d["score"]["listed_both_sides"] > 0 for d in ds)
        by = statistics.median(d["iblt_bytes"] for d in ds)
        rows.append([p, r, ds[0]["n"], n, det, raw, canc, f"{cp_upper(n - canc, n):.3f}", undetected, both, int(by),
                     f"{by / D[p]:.2f}", lines(no for no, _ in g)])
    write(out, "iblt_trials", ["pair", "ratio_n_over_D", "n", "trials", "list_exit_0", "exact_raw", "exact_cancelled",
                               "failure_upper95 (cancelled)", "undetected_failures", "trials_with_both_side_listings",
                               "median_bytes", "bytes_per_difference", "src"], rows)

    # 6. IBLT resend protocol from an underestimate
    rows = []
    groups = defaultdict(list)
    for no, d in recs["iblt_protocol"]:
        groups[(d["pair"], d["f"])].append((no, d))
    for (p, f), g in sorted(groups.items()):
        ds = [d for _, d in g]
        rows.append([p, f, len(ds), f"{statistics.mean(d['rounds'] for d in ds):.2f}", max(d["rounds"] for d in ds),
                     sum(d["last"]["score"]["cancelled_exact"] for d in ds),
                     f"{statistics.mean(d['total_bytes'] for d in ds):.0f}",
                     f"{statistics.mean(d['total_bytes'] for d in ds) / D[p]:.2f}", lines(no for no, _ in g)])
    write(out, "iblt_protocol", ["pair", "dest_over_D", "reps", "mean_rounds", "max_rounds", "exact_cancelled",
                                 "mean_total_bytes", "bytes_per_difference", "src"], rows)

    # 7. S collisions: distinct values vs distinct canonical k-mers, every KP-10 run
    rows = []
    kset = {d["genome"]: d["out"]["distinct_canonical"] for _, d in recs["kset"]}
    allk = recs["kp10_timed"] + recs["kp10_seed"] + recs["kp10_dest"]
    for p in ("giga", "draft"):
        rs = [(no, d) for no, d in allk if d["label"].endswith("_" + p)]
        if not rs:
            continue
        a, b = rs[0][1]["A"].rsplit("/", 1)[1][:-3], rs[0][1]["B"].rsplit("/", 1)[1][:-3]
        bad = [no for no, d in rs if d["values_A"] != kset[a] or d["values_B"] != kset[b]]
        rows.append([p, len(rs), kset[a], kset[b], len(bad), lines(no for no, _ in rs)])
    write(out, "collisions", ["pair", "kp10_runs", "canonical_A", "canonical_B",
                              "runs_where_distinct_S_values_differ_from_canonical_counts", "src"], rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
