"""Summarize round 3's owed runs (recon/round3_owed_job.sh) into TSV tables, each row citing its run.log line(s).

    python3 recon/round3_owed_summarize.py <run.log> [<run.log> ...] <out dir>

Reads the "@@ <kind> <json>" lines of either job (kp, iblt) and writes, for those present:
  kp_seeds.tsv     run 1 + 4: OS-seed runs per (pair, rule): rounds, exactness, bytes per difference per party,
                   two-way, one-way, stored and packed
  kp_protocols.tsv the fixed-seed one-way vs two-way pair (same buckets split?)
  kp_sweep.tsv     run 1: every fixed-seed sweep point (m, first-round buckets, load) and the exact-rule D_est sweep
  kp63.tsv         run 6: k = 63 at the 0.8 m rule
  kmc63.tsv        the k = 63 truth
  iblt_validate.tsv, kmp_trials.tsv (run 2 + 3), iblt_proto.tsv (run 5), ext.tsv + ext_trace.tsv (run 7)
Failure-rate bounds are exact one-sided 95% Clopper-Pearson. Packed IBLT bytes: 3 chunk cells of (62 + b) bits, b the
counter width that carries every difference counter seen at that (set, ratio) (max over its trials); packed sketch
bytes are kp10c's own (61-bit entries, 32-bit counts, see its header). Prints every table as Markdown.
"""

from __future__ import annotations

import json
import math
import pathlib
import statistics
import sys
from collections import defaultdict


def cp_upper(fail: int, n: int, alpha: float = 0.05) -> float:
    if n == 0 or fail == n:
        return 1.0

    def cdf(p: float) -> float:
        return sum(math.comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(fail + 1))

    lo, hi = fail / n, 1.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if cdf(mid) > alpha:
            lo = mid
        else:
            hi = mid
    return hi


def load(paths: list[pathlib.Path]) -> dict[str, list[tuple[str, dict]]]:
    recs: dict[str, list[tuple[str, dict]]] = defaultdict(list)
    for path in paths:
        tag = path.parent.name
        for no, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
            if line.startswith("@@ "):
                _, kind, js = line.split(" ", 2)
                recs[kind].append((f"{tag}/run.log:{no}", json.loads(js)))
    return recs


def span(srcs: list[str]) -> str:
    """'tag/run.log:a-b (n lines)' for a list of sources from one file, else a ';'-joined list."""
    if not srcs:
        return ""
    files = {s.rsplit(":", 1)[0] for s in srcs}
    if len(files) == 1 and len(srcs) > 2:
        nums = sorted(int(s.rsplit(":", 1)[1]) for s in srcs)
        return f"{files.pop()}:{nums[0]}-{nums[-1]} ({len(nums)} lines)"
    return "; ".join(srcs)


def write(out: pathlib.Path, name: str, rows: list[dict]) -> None:
    if not rows:
        return
    cols = list(rows[0])
    text = "\t".join(cols) + "\n" + "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows)
    (out / name).write_text(text)
    print(f"\n### {name}\n\n| " + " | ".join(cols) + " |\n|" + "---|" * len(cols))
    for r in rows:
        print("| " + " | ".join(str(r[c]) for c in cols) + " |")


def f3(x: float) -> str:
    return f"{x:.3f}"


def pair_of(label: str) -> str:
    return label.rsplit("_", 1)[1]


def kp_tables(recs, out):
    rows, groups = [], defaultdict(list)
    for src, j in recs.get("kp", []):
        lab = j["label"]
        if lab.startswith(("exact_os_", "pow2_os_")):
            groups[(pair_of(lab), lab.split("_")[0])].append((src, j))
    for (p, rule), v in sorted(groups.items()):
        js = [j for _, j in v]
        n = len(js)
        ex = sum(j["exact"] for j in js)
        rounds = defaultdict(int)
        for j in js:
            rounds[j["rounds"]] += 1
        m = lambda k: statistics.mean(j[k] for j in js)  # noqa: E731
        rows.append({"pair": p, "rule": rule, "runs": n, "exact": ex, "failure_upper95": f3(cp_upper(n - ex, n)),
                     "late_failures": sum(j["late_failures"] for j in js),
                     "buckets0": "/".join(sorted({str(j["buckets0"]) for j in js})),
                     "rounds": " ".join(f"{k}:{c}" for k, c in sorted(rounds.items())),
                     "mean_round2_buckets": f"{statistics.mean(j['buckets_sent_per_round'][1] if j['rounds'] > 1 else 0 for j in js):.1f}",
                     "decoder": "/".join(sorted({j["decoder"] for j in js})),
                     "per_party_B_diff": f3(m("bytes_per_true_difference")),
                     "per_party_min_max": f"{min(j['bytes_per_true_difference'] for j in js):.3f}-{max(j['bytes_per_true_difference'] for j in js):.3f}",
                     "per_party_packed": f3(m("packed_bytes_per_true_difference")),
                     "twoway_B_diff": f3(m("twoway_bytes_per_true_difference")),
                     "twoway_packed": f3(m("twoway_packed_bytes_per_true_difference")),
                     "oneway_B_diff": f3(m("oneway_bytes_per_true_difference")),
                     "oneway_packed": f3(m("oneway_packed_bytes_per_true_difference")),
                     "oneway_reply_kmers": "/".join(sorted({str(j["oneway_reply_kmers"]) for j in js})),
                     "mean_cpu_reconcile_s": f"{m('cpu_reconcile'):.1f}",
                     "src": span([s for s, _ in v])})
    write(out, "kp_seeds.tsv", rows)

    rows = []
    fx = {j["label"]: (s, j) for s, j in recs.get("kp", []) if j["label"].startswith("fixed_exact_")}
    for p in ("giga", "draft"):
        a, b = fx.get(f"fixed_exact_oneway_{p}"), fx.get(f"fixed_exact_twoway_{p}")
        if a and b:
            rows.append({"pair": p, "sent_oneway": a[1]["buckets_sent_per_round"], "sent_twoway": b[1]["buckets_sent_per_round"],
                         "failed_oneway": a[1]["buckets_failed_per_round"], "failed_twoway": b[1]["buckets_failed_per_round"],
                         "same_rounds": a[1]["buckets_failed_per_round"] == b[1]["buckets_failed_per_round"]
                         and a[1]["buckets_sent_per_round"] == b[1]["buckets_sent_per_round"],
                         "exact_both": a[1]["exact"] and b[1]["exact"], "src": f"{a[0]}; {b[0]}"})
    write(out, "kp_protocols.tsv", rows)

    rows = []
    for src, j in recs.get("kp", []):
        lab = j["label"]
        if not lab.startswith(("sweep_", "dest_exact_")):
            continue
        D = j["truth_A_minus_B"] + j["truth_B_minus_A"]
        rows.append({"pair": pair_of(lab), "kind": "dest" if lab.startswith("dest") else "buckets", "m": j["m"],
                     "dest_over_D": f3(j["dest"] / D), "buckets0": j["buckets0"],
                     "load_D_over_N0m": f3(D / (j["buckets0"] * j["m"])), "rounds": j["rounds"],
                     "sent": j["buckets_sent_per_round"], "exact": j["exact"],
                     "per_party_B_diff": f3(j["bytes_per_true_difference"]),
                     "per_party_packed": f3(j["packed_bytes_per_true_difference"]),
                     "twoway_B_diff": f3(j["twoway_bytes_per_true_difference"]),
                     "oneway_B_diff": f3(j["oneway_bytes_per_true_difference"]),
                     "oneway_packed": f3(j["oneway_packed_bytes_per_true_difference"]), "src": src})
    write(out, "kp_sweep.tsv", rows)

    rows = []
    v = [(s, j) for s, j in recs.get("kp63", []) if j["label"].startswith("k63_os_")]
    if v:
        js = [j for _, j in v]
        n, ex = len(js), sum(j["exact"] for j in js)
        rounds = defaultdict(int)
        for j in js:
            rounds[j["rounds"]] += 1
        rows.append({"run": "30 OS seeds, 0.8m rule", "runs": n, "exact": ex, "failure_upper95": f3(cp_upper(n - ex, n)),
                     "buckets0": "/".join(sorted({str(j["buckets0"]) for j in js})),
                     "rounds": " ".join(f"{k}:{c}" for k, c in sorted(rounds.items())),
                     "D": js[0]["truth_A_minus_B"] + js[0]["truth_B_minus_A"],
                     "per_party_B_diff": f3(statistics.mean(j["bytes_per_true_difference"] for j in js)),
                     "per_party_bytes_min_max": f"{min(j['sketch_bytes'] for j in js)}-{max(j['sketch_bytes'] for j in js)}",
                     "oneway_B_diff": f3(statistics.mean(j["oneway_bytes_per_true_difference"] for j in js)),
                     "oneway_packed": f3(statistics.mean(j["oneway_packed_bytes_per_true_difference"] for j in js)),
                     "missing": sum(j["missing"] for j in js), "spurious": sum(j["spurious"] for j in js),
                     "src": span([s for s, _ in v])})
    for s, j in recs.get("kp63", []):
        if j["label"] == "k63_fixed":
            rows.append({"run": f"fixed seed, {j['buckets0']} buckets", "runs": 1, "exact": int(j["exact"]),
                         "failure_upper95": "", "buckets0": j["buckets0"], "rounds": str(j["rounds"]),
                         "D": j["truth_A_minus_B"] + j["truth_B_minus_A"],
                         "per_party_B_diff": f3(j["bytes_per_true_difference"]), "per_party_bytes_min_max": str(j["sketch_bytes"]),
                         "oneway_B_diff": f3(j["oneway_bytes_per_true_difference"]),
                         "oneway_packed": f3(j["oneway_packed_bytes_per_true_difference"]),
                         "missing": j["missing"], "spurious": j["spurious"], "src": s})
    write(out, "kp63.tsv", rows)
    write(out, "kmc63.tsv", [{**{k: j[k] for k in ("tool", "A", "B", "k", "distinct_A", "distinct_B", "A_minus_B", "B_minus_A", "D")},
                              "src": s} for s, j in recs.get("kmc63", [])])


def iblt_tables(recs, out):
    rows = []
    for s, j in recs.get("iblt_validate", []):
        pg = j["iblt_grow_peel"]
        rows.append({"i": j["i"], "ratio": j["ratio"], "n": j["n"], "seed": j["seed"], "a_bytes": j["a_iblt_bytes"],
                     "d_full_eq_d_diff": j["d_full_equals_d_difference_only"],
                     "write_eq_kmpeeler_a": j["iblt_grow_write_equals_kmpeeler_a"], "kmp_list_exit": j["kmp_list_exit"],
                     "kmp_cancelled_exact": j["kmp_score"]["cancelled_exact"],
                     "kmp_listed": j["kmp_score"]["listed_A_minus_B"] + j["kmp_score"]["listed_B_minus_A"],
                     "peel_success": pg["success"], "peel_listed": pg["listed_i"] + pg["listed_j"],
                     "agree": (pg["success"] == 1) == (j["kmp_list_exit"] == 0)
                     and pg["listed_i"] == j["kmp_score"]["listed_A_minus_B"] and pg["listed_j"] == j["kmp_score"]["listed_B_minus_A"],
                     "max_abs_counter": j["stat"]["max_abs_counter"], "wall_s": j["wall_s"], "src": s})
    write(out, "iblt_validate.tsv", rows)

    g = defaultdict(list)
    for s, j in recs.get("kmp_trial", []):
        g[(j["set"], j["ratio"])].append((s, j))
    rows, seeds = [], []
    for (st, r), v in sorted(g.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        js = [j for _, j in v]
        seeds += [j["seed"] for j in js]
        n = len(js)
        ok = sum(j["list_exit"] == 0 and j["score"]["cancelled_exact"] for j in js)
        undet = sum(j["list_exit"] == 0 and not j["score"]["cancelled_exact"] for j in js)
        D = js[0]["D"]
        b = max(j["stat"]["counter_bits"] for j in js)
        cells = 3 * js[0]["stat"]["chunk"]
        rows.append({"set": st, "ratio": r, "D": D, "n": js[0]["n"], "trials": n, "exact": ok, "undetected": undet,
                     "both_side_listings": sum(j["score"]["listed_both_sides"] > 0 for j in js),
                     "failure_upper95": f3(cp_upper(n - ok, n)),
                     "listed_frac_when_failed": f3(statistics.mean((j["score"]["listed_A_minus_B"] + j["score"]["listed_B_minus_A"]) / D
                                                                   for j in js if j["list_exit"] != 0)) if n - ok else "",
                     "bytes": js[0]["iblt_bytes"], "B_diff": f3(js[0]["iblt_bytes"] / D),
                     "max_abs_counter": max(j["stat"]["max_abs_counter"] for j in js), "counter_bits": b,
                     "packed_B_diff": f3(cells * (2 * 31 + b) / 8 / D), "src": span([s for s, _ in v])})
    write(out, "kmp_trials.tsv", rows)
    if seeds:
        trip = [x + d for x in seeds for d in range(3)]
        print(f"\nkm-peeler trials: {len(seeds)} root seeds, hash functions shared between trials: {len(trip) - len(set(trip))}")

    g = defaultdict(list)
    for s, j in recs.get("iblt_proto", []):
        g[(j["label"].split("_")[0], j["proto"], j["f"])].append((s, j))
    rows = []
    for (p, pr, f), v in sorted(g.items()):
        js = [j for _, j in v]
        n = len(js)
        ex = sum(j["exact"] for j in js)
        rounds = defaultdict(int)
        for j in js:
            rounds[j["rounds"]] += 1
        b = max(j["counter_bits"] for j in js)
        rows.append({"pair": p, "proto": pr, "f": f, "reps": n, "exact": ex, "undetected": sum(j["undetected_wrong"] for j in js),
                     "failure_upper95": f3(cp_upper(n - ex, n)),
                     "rounds": " ".join(f"{k}:{c}" for k, c in sorted(rounds.items())),
                     "B_diff_mean": f3(statistics.mean(j["bytes_per_diff"] for j in js)),
                     "B_diff_min_max": f"{min(j['bytes_per_diff'] for j in js):.2f}-{max(j['bytes_per_diff'] for j in js):.2f}",
                     # repacked with the widest counter seen at this (pair, proto, f)
                     "packed_B_diff_mean": f3(statistics.mean(j["cells_sent"] * (62 + b) / 8 / j["D"] + (8 * (j["rounds"] - 1) if pr == "subtract" else 0) / j["D"]
                                                              for j in js)),
                     "counter_bits": b, "src": span([s for s, _ in v])})
    write(out, "iblt_proto.tsv", rows)
    all_seeds = [r["seed"] for _, j in recs.get("iblt_proto", []) for r in j["per_round"] if j["proto"] != "grow" or r is j["per_round"][0]]
    if all_seeds:
        trip = [x + d for x in all_seeds for d in range(3)]
        print(f"iblt_grow tables: {len(all_seeds)} root seeds, hash functions shared between tables: {len(trip) - len(set(trip))}")

    rows = []
    for s, j in recs.get("ext_set", []):
        rows.append({"what": "set", "pair": j["pair"], "trial": "", "D_ext": j["D_ext"], "detail": json.dumps(j), "src": s})
    for s, j in recs.get("ext_trial", []):
        sc = j["score"]
        rows.append({"what": "trial", "pair": j["pair"], "trial": j["trial"], "D_ext": j["D_ext"],
                     "detail": f"exit {j['list_exit']}; listed {j['listed']}; bytes {j['iblt_bytes']}; "
                               f"B per k-mer diff {j['iblt_bytes'] / sc['truth']:.3f}; missing {sc['missing']} "
                               f"({sc['missing_A_minus_B']} A-side); spurious {sc['spurious']} "
                               f"({sc['spurious'] / sc['truth']:.3f} of D); wall {j['wall_s']} s", "src": s})
    for kind in ("ext_miss_md5", "ext_trace"):
        for s, j in recs.get(kind, []):
            rows.append({"what": kind, "pair": j["pair"], "trial": "", "D_ext": "", "detail": json.dumps(j), "src": s})
    write(out, "ext.tsv", rows)
    write(out, "tally_meta.tsv", [{"detail": json.dumps(j), "src": s} for s, j in recs.get("tally", [])])
    write(out, "sets.tsv", [{"kind": k, "detail": json.dumps(j), "src": s} for k in ("kset", "diff", "subset") for s, j in recs.get(k, [])])


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        return 2
    out = pathlib.Path(sys.argv[-1])
    out.mkdir(parents=True, exist_ok=True)
    recs = load([pathlib.Path(p) for p in sys.argv[1:-1]])
    kp_tables(recs, out)
    iblt_tables(recs, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
