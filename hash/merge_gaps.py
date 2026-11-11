"""Merge the long gaps logged by `stats_seeds ... gaps=<file>` into distinct loci across all runs.

    python3 hash/merge_gaps.py <gaps.tsv> [--out loci.tsv] [--json summary.json]

Every hash is a function of the k-mer, so one repeat gaps in run after run, and gap EVENTS (rows) are not independent
(round-1 attack L2 F1, adjudication L2 F1; the pseudo_replication card). A locus is a maximal union of gap intervals
[start, end) on one record that overlap or touch, pooled over every method and seed (a global merge). Per locus:
gaps, runs (method, seed) and methods hitting it, and whether it lies in soft-masked repeats, meaning every gap in it
has at least half its bases lowercase (the same per-gap rule stats_seeds uses for long_outside_repeats). Reported:
gaps and loci in total and per method (global loci touched by that method), loci in repeats, loci hit by every
method, and the largest locus. Standard library only.
"""

from __future__ import annotations

import argparse
import collections
import json
import sys


def read_gaps(path: str) -> list[dict]:
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            v = dict(zip(header, line.rstrip("\n").split("\t")))
            rows.append({"method": v["method"], "seed": int(v["seed"]), "record": v["record"], "start": int(v["start"]),
                         "end": int(v["end"]), "gap_kmers": int(v["gap_kmers"]), "lower": int(v["lowercase_bases"])})
    return rows


def merge(rows: list[dict]) -> list[dict]:
    by_rec = collections.defaultdict(list)
    for r in rows:
        by_rec[r["record"]].append(r)
    loci = []
    for rec in sorted(by_rec):
        cur = None
        for r in sorted(by_rec[rec], key=lambda r: (r["start"], r["end"])):
            if cur is not None and r["start"] <= cur["end"]:
                cur["end"] = max(cur["end"], r["end"])
                cur["gaps"].append(r)
            else:
                cur = {"record": rec, "start": r["start"], "end": r["end"], "gaps": [r]}
                loci.append(cur)
    return loci


def summarize(rows: list[dict], loci: list[dict]) -> dict:
    methods = sorted({r["method"] for r in rows})
    runs = sorted({(r["method"], r["seed"]) for r in rows})
    per_method = {m: {"gaps": sum(r["method"] == m for r in rows),
                      "loci": sum(any(g["method"] == m for g in L["gaps"]) for L in loci)} for m in methods}
    in_rep = [all(2 * g["lower"] >= g["end"] - g["start"] for g in L["gaps"]) for L in loci]
    all_m = [len({g["method"] for g in L["gaps"]}) == len(methods) for L in loci]
    big = max(loci, key=lambda L: len(L["gaps"])) if loci else None
    return {"gaps": len(rows), "runs_with_gaps": len(runs), "methods": methods, "loci": len(loci),
            "records": len({L["record"] for L in loci}), "loci_in_repeats": sum(in_rep),
            "loci_outside_repeats": len(loci) - sum(in_rep), "loci_hit_by_all_methods": sum(all_m),
            "per_method": per_method,
            "largest_locus": None if big is None else {
                "where": f"{big['record']}:{big['start']}-{big['end']}", "bases": big["end"] - big["start"],
                "gaps": len(big["gaps"]), "runs": len({(g["method"], g["seed"]) for g in big["gaps"]})},
            "gaps_per_locus": round(len(rows) / max(len(loci), 1), 2)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("gaps")
    ap.add_argument("--out")
    ap.add_argument("--json")
    a = ap.parse_args()
    rows = read_gaps(a.gaps)
    loci = merge(rows)
    if a.out:
        with open(a.out, "w") as fh:
            fh.write("record\tstart\tend\tgaps\truns\tmethods\tin_repeats\n")
            for L in loci:
                ms = sorted({g["method"] for g in L["gaps"]})
                rep = all(2 * g["lower"] >= g["end"] - g["start"] for g in L["gaps"])
                fh.write(f"{L['record']}\t{L['start']}\t{L['end']}\t{len(L['gaps'])}\t"
                         f"{len({(g['method'], g['seed']) for g in L['gaps']})}\t{','.join(ms)}\t{'yes' if rep else 'no'}\n")
    s = summarize(rows, loci)
    print(json.dumps(s, indent=1))
    if a.json:
        with open(a.json, "w") as fh:
            fh.write(json.dumps(s, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
