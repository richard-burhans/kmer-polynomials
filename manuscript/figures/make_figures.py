#!/usr/bin/env python3
"""Draw the manuscript's figures, as SVG and PDF, from committed results files.

    make -C manuscript figures      (runs: uv run --no-project --with matplotlib==3.9.2 python figures/make_figures.py)

Every plotted number is read from a file under results/, never typed here:
  fig1_encoding     diagram only. Base values and the reverse-complement rule are the paper's definitions, and the
                    script asserts that F(rc(s)) = R(s) for its example before drawing it.
  fig2_collapse     results/kp7_symmetric_check_2026-11-10.txt (k = 7 block)
  fig3_speed        results/paper_bench_2026-11-12/ raw counters, through manuscript/paper_bench_table.py
  fig4_sketch       results/kp10_vs_shibuya_2026-11-11/laila_b69c5c76/tables/{kp10_dest,iblt_protocol}.tsv (giga)
Output is byte-for-byte deterministic: fixed SVG hash salt, no dates in either format.
"""
import csv
import pathlib
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
MANUSCRIPT = HERE.parent
RESULTS = MANUSCRIPT.parent / "results"
sys.path.insert(0, str(MANUSCRIPT))
import paper_bench_table  # noqa: E402

# Two categorical slots (blue, orange) from the dataviz reference palette, validated together, plus ink tokens.
OURS, OTHER = "#2a78d6", "#eb6834"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
MM = 1 / 25.4

plt.rcParams.update({
    "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans", "font.size": 7,
    "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "lines.linewidth": 1.4,
    "legend.frameon": False, "svg.hashsalt": "kmer-polynomials", "svg.fonttype": "path", "pdf.fonttype": 42,
})


def panel(ax, letter):
    ax.text(-0.02, 1.02, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", ha="right", va="bottom")


def save(fig, name):
    fig.savefig(HERE / f"{name}.svg", metadata={"Date": None})
    fig.savefig(HERE / f"{name}.pdf", metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)


# ---------- Figure 1: the encoding, and strand symmetry of the product ----------

# A base value is (sign, unit) with unit "1" or "w" (omega): A = 1, C = w, G = -w, T = -1, as paper.md defines them.
BASE = {"A": (1, "1"), "C": (1, "w"), "G": (-1, "w"), "T": (-1, "1")}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def rc(s):
    return "".join(COMP[b] for b in reversed(s))


def coeffs_F(s):
    return [BASE[b] for b in s]


def coeffs_R(s):
    """R = -x^(k-1) p(1/x): coefficient t of R is minus coefficient k-1-t of p."""
    return [(-sg, u) for sg, u in reversed(coeffs_F(s))]


def poly_tex(cs):
    out = []
    for t, (sg, u) in enumerate(cs):
        mono = "" if t == 0 else ("x" if t == 1 else f"x^{t}")
        term = (r"\omega " if u == "w" else ("" if mono else "1")) + mono
        out.append(("-" if sg < 0 else ("+" if out else "")) + term)
    return "$" + " ".join(out) + "$"


def fig1():
    s = "AACG"
    assert coeffs_F(rc(s)) == coeffs_R(s) and coeffs_R(rc(s)) == coeffs_F(s), "reverse complement must swap F and R"
    fig = plt.figure(figsize=(170 * MM, 62 * MM))
    ax = fig.add_axes([0.03, 0.06, 0.30, 0.86])
    for x, letter in ((0.01, "A"), (0.36, "B")):
        fig.text(x, 0.97, letter, fontsize=11, fontweight="bold", va="top")
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(Circle((0, 0), 1, fill=False, lw=0.6, ec=MUTED, ls=(0, (2, 2))))
    ax.axhline(0, color=GRID, lw=0.6, zorder=0)
    ax.axvline(0, color=GRID, lw=0.6, zorder=0)
    pos = {"A": (1, 0), "C": (-0.5, 3 ** 0.5 / 2), "G": (0.5, -3 ** 0.5 / 2), "T": (-1, 0)}
    lab = {"A": "A = 1", "C": r"C = $\omega$", "G": r"G = $-\omega$", "T": "T = $-1$"}
    off = {"A": (0.12, 0.10), "C": (-0.08, 0.14), "G": (0.10, -0.16), "T": (-0.12, 0.10)}
    ha = {"A": "left", "C": "right", "G": "left", "T": "right"}
    for a, b in (("A", "T"), ("C", "G")):  # complement: reflection through the origin
        (x0, y0), (x1, y1) = pos[a], pos[b]
        ax.plot([x0, x1], [y0, y1], color=INK2, lw=0.9, ls=(0, (4, 2)), zorder=1)
    for a, b in (("A", "G"), ("C", "T")):  # transitions, squared distance 1
        (x0, y0), (x1, y1) = pos[a], pos[b]
        ax.plot([x0, x1], [y0, y1], color=OURS, lw=1.4, zorder=1)
    ax.plot([pos["A"][0], pos["C"][0]], [pos["A"][1], pos["C"][1]], color=MUTED, lw=0.8, zorder=1)
    for b, (x, y) in pos.items():
        ax.plot(x, y, "o", ms=6, mfc="white" if b in "CT" else INK, mec=INK, mew=1.0, zorder=3)
        ax.text(x + off[b][0], y + off[b][1], lab[b], ha=ha[b], va="center", fontsize=8)
    ax.plot(0, 0, "+", color=INK2, ms=5, mew=0.8)
    ax.text(0.86, -0.50, "transition\n$|A-G|^2 = 1$", color=OURS, fontsize=6.5, ha="left", va="center")
    ax.text(0.30, 0.56, "$|A-C|^2 = 3$", color=INK2, fontsize=6.5, ha="left", va="bottom")
    ax.text(-0.05, -0.10, "complement\n= negation", color=INK2, fontsize=6.5, ha="right", va="top")
    ax.set_xlim(-1.55, 1.75)
    ax.set_ylim(-1.25, 1.25)

    bx = fig.add_axes([0.37, 0.06, 0.62, 0.86])
    bx.axis("off")
    bx.set_xlim(0, 1)
    bx.set_ylim(0, 1)
    rows = [(0.80, "$s$", s), (0.30, r"$\mathrm{rc}(s)$", rc(s))]
    for y, name, kmer in rows:
        bx.text(0.0, y + 0.07, name, fontsize=8, va="center")
        for i, b in enumerate(kmer):
            bx.text(0.09 + 0.045 * i, y + 0.07, b, fontsize=9, family="DejaVu Sans Mono", ha="center", va="center")
            sg, u = BASE[b]
            bx.text(0.09 + 0.045 * i, y - 0.02, ("$-$" if sg < 0 else "") + (r"$\omega$" if u == "w" else "1"),
                    fontsize=6.5, color=INK2, ha="center", va="center")
        bx.text(0.27, y + 0.10, "$F$ = " + poly_tex(coeffs_F(kmer)), fontsize=7.5, va="center")
        bx.text(0.27, y - 0.02, "$R$ = " + poly_tex(coeffs_R(kmer)), fontsize=7.5, va="center")
    # reverse complementation swaps F and R: an X between the two rows
    for y0, y1 in ((0.90, 0.28), (0.78, 0.40)):
        bx.add_patch(FancyArrowPatch((0.565, y0), (0.665, y1), arrowstyle="<->", mutation_scale=6, lw=0.6,
                                     color=MUTED, shrinkA=1, shrinkB=1))
    bx.text(0.615, 0.59, "rc swaps\n$F$ and $R$", fontsize=6.5, color=INK2, ha="center", va="center",
            bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none"))
    for y in (0.84, 0.34):
        bx.add_patch(FancyArrowPatch((0.67, y), (0.73, 0.59), arrowstyle="-|>", mutation_scale=7, lw=0.8,
                                     color=INK2))
    bx.text(0.865, 0.59, r"$S = (F+\lambda)(R+\lambda)$" + "\nthe same value\nfrom both strands", fontsize=7.5,
            ha="center", va="center", bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=OURS, lw=1.0))
    bx.text(0.0, 0.03, "Coefficients of $x^0, x^1, \\ldots, x^{k-1}$, here $k = 4$.", fontsize=6.5, color=INK2)
    save(fig, "fig1_encoding")


# ---------- Figure 2: distinct values over the 7-mer classes ----------

def fig2():
    src = RESULTS / "kp7_symmetric_check_2026-11-10.txt"
    block = src.read_text().split("k=8:")[0]
    classes = int(re.search(r"(\d+) double-stranded classes", block).group(1))
    vals = [int(re.search(pat, block).group(1)) for pat in (r"F\+R\s+distinct values: (\d+)",
                                                          r"F\*R\s+distinct values: (\d+)",
                                                          r"\(F\+lam\)\(R\+lam\)\s+distinct values: (\d+)")]
    names = ["$F + R$", r"$F \cdot R$", r"$(F+\lambda)(R+\lambda)$"]
    fig, ax = plt.subplots(figsize=(85 * MM, 58 * MM))
    fig.subplots_adjust(left=0.17, right=0.97, bottom=0.14, top=0.90)
    xs = range(len(vals))
    ax.bar(xs, vals, width=0.55, color=[MUTED, MUTED, OURS], edgecolor=INK, linewidth=0.6,
           hatch=["", "", ""])
    ax.axhline(classes, color=INK, lw=0.8, ls=(0, (4, 2)))
    ax.text(-0.42, classes + 150, f"{classes:,} double-stranded classes", fontsize=6.5, va="bottom", ha="left")
    for x, v in zip(xs, vals):
        ax.text(x, v + 150 if v < classes else v - 600, f"{v:,}", ha="center", va="bottom" if v < classes else "top",
                fontsize=7, color=INK if v < classes else "white", fontweight="bold" if v == classes else "normal")
    ax.set_xticks(list(xs), names)
    ax.set_ylabel("distinct values, all 7-mers")
    ax.set_ylim(0, classes * 1.12)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.tick_params(axis="x", length=0)
    save(fig, "fig2_collapse")
    return classes, vals


# ---------- Figure 3: cycles per input k-mer against sampling rate ----------

def rate(meth):
    num, den = meth.split()[-1].split("/")
    return int(num) / int(den)


def fig3():
    _, rows, _, _ = paper_bench_table.compute(RESULTS / "paper_bench_2026-11-12")

    def series(b, branchy, batched):
        pts = [(rate(m), rows[(b, m)][0]) for m in branchy] + [(rate(m), rows[(b, m)][0]) for m in batched]
        return sorted(pts)

    fig, (a, b) = plt.subplots(1, 2, figsize=(170 * MM, 68 * MM), sharey=True)
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.17, top=0.86, wspace=0.10)
    xlim = (6e-4, 1.6)

    # A: bench_e, where the paper times every-k-mer hashing and the sampler's break-even
    nt_e = rows[("bench_e", "ntHash batched (ntbatch)")][0]
    every_e = rows[("bench_e", "e sym_seq L2")][0]
    samp_e = series("bench_e", ["e branchy 1/1000", "e branchy 1/100"], ["e batched 1/16", "e batched 1/4"])
    slope, b0 = paper_bench_table.slope_fit(rows, "bench_e", "e batched")
    be = (nt_e - b0) / slope
    a.axhline(nt_e, color=INK2, lw=1.0, ls=(0, (1, 1.5)), label=f"ntHash, batched ({nt_e:.2f})")
    a.plot(*zip(*samp_e), "-o", color=OURS, ms=4.5, mec="white", mew=0.6, label="Eisenstein sampler")
    a.plot(1, every_e, "D", color=OURS, mfc="white", mew=1.2, ms=5, label=f"Eisenstein, every k-mer ({every_e:.2f})")
    a.plot(be, nt_e, "o", ms=7, mfc="none", mec=INK, mew=1.0, zorder=5)
    a.axvline(be, color=INK, lw=0.6, ls=(0, (3, 2)), ymax=nt_e / 30)
    a.annotate(f"break-even 1/{1 / be:.1f}", (be, nt_e), xytext=(be * 1.25, 5.0), fontsize=6.5, ha="left",
               arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2))
    a.set_title("program of Table 2", loc="left", fontsize=7, color=INK2)

    # B: bench_kept, where the sampler meets the universal-hash baselines in one program
    nt_k = rows[("bench_kept", "ntHash batched")][0]
    uh_all = rows[("bench_kept", "uh all")][0]
    samp_k = series("bench_kept", ["lib-e 8 br 1/1000", "lib-e 8 br 1/100"], ["lib-e 8 bat 1/16", "lib-e 8 bat 1/4"])
    uh_k = series("bench_kept", ["uh 0 br 1/1000", "uh 0 br 1/100"], ["uh 0 bat 1/16", "uh 0 bat 1/4"])
    b.axhline(nt_k, color=INK2, lw=1.0, ls=(0, (1, 1.5)), label=f"ntHash, batched ({nt_k:.2f})")
    b.axhline(uh_all, color=OTHER, lw=1.0, ls=(0, (5, 1.5, 1, 1.5)),
              label=f"canonical code + universal hash, every k-mer ({uh_all:.2f})")
    b.plot(*zip(*samp_k), "-o", color=OURS, ms=4.5, mec="white", mew=0.6, label="Eisenstein sampler")
    b.plot(*zip(*uh_k), "--s", color=OTHER, ms=4.5, mec="white", mew=0.6, label="same filter, then universal hash")
    b.set_title("program of Table 3", loc="left", fontsize=7, color=INK2)

    for ax, letter in ((a, "A"), (b, "B")):
        panel(ax, letter)
        ax.set_xscale("log")
        ax.set_xlim(*xlim)
        ax.set_ylim(0, 30)
        ax.set_xlabel("sampling rate (fraction of k-mers hashed)")
        ax.grid(axis="y", color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        ax.legend(loc="upper left", handlelength=2.6, borderaxespad=0.2)
    a.set_ylabel("cycles per input k-mer")
    save(fig, "fig3_speed")
    return dict(nt_e=nt_e, every_e=every_e, samp_e=samp_e, be=be, slope=slope, nt_k=nt_k, uh_all=uh_all,
                samp_k=samp_k, uh_k=uh_k)


# ---------- Figure 4: the two-round sketch, and cost when the difference is underestimated ----------

def read_tsv(name):
    with open(RESULTS / "kp10_vs_shibuya_2026-11-11/laila_b69c5c76/tables" / name) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def bucket(ax, x, y, w, h, hatch=None, fc="white", ec=INK, lw=0.8):
    ax.add_patch(Rectangle((x, y), w, h, fc=fc, ec=ec, lw=lw, hatch=hatch))
    for j in range(1, 4):  # count and power sums, stacked
        ax.plot([x, x + w], [y + j * h / 4, y + j * h / 4], color=ec, lw=0.3)


def fig4():
    fracs = (0.125, 0.25, 0.5, 1.0)
    kp = {float(r["dest_over_D"]): (float(r["bytes_per_difference"]), int(r["rounds"]))
          for r in read_tsv("kp10_dest.tsv") if r["pair"] == "giga" and r["exact"] == "True"}
    ib = {}
    for r in read_tsv("iblt_protocol.tsv"):
        if r["pair"] == "giga":
            assert float(r["mean_rounds"]) == int(r["max_rounds"]), "IBLT rounds must not vary across repetitions"
            assert r["exact_cancelled"] == r["reps"], "every IBLT repetition must end exact"
            ib[float(r["dest_over_D"])] = (float(r["bytes_per_difference"]), int(r["max_rounds"]))
    kp = [(f, *kp[f]) for f in fracs]
    ib = [(f, *ib[f]) for f in fracs]

    fig = plt.figure(figsize=(170 * MM, 66 * MM))
    ax = fig.add_axes([0.02, 0.04, 0.50, 0.86])
    fig.text(0.01, 0.97, "A", fontsize=11, fontweight="bold", va="top")
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6.2)
    n, w, h, gap, x0 = 8, 0.78, 1.25, 0.18, 1.6
    bad = 5
    ax.text(0.0, 4.55, "round 1", fontsize=7.5, fontweight="bold", va="center")
    ax.text(x0, 6.05, r"sketch$(A)$ $-$ sketch$(B)$: per bucket, the count and $\sum h^j$, $j = 1, \ldots, m + c$",
            fontsize=6.5, va="center")
    for i in range(n):
        x = x0 + i * (w + gap)
        is_bad = i == bad
        bucket(ax, x, 3.9, w, h, hatch="////" if is_bad else None, fc="white", ec=OTHER if is_bad else INK,
               lw=1.2 if is_bad else 0.8)
        ax.text(x + w / 2, 3.62, "> m" if is_bad else "\u2713", fontsize=7 if is_bad else 8, ha="center",
                va="top", color=OTHER if is_bad else INK2)
    ax.text(x0, 5.45, "buckets: ranges of $S$", fontsize=6.5, color=INK2, va="center")
    xb = x0 + bad * (w + gap) + w / 2
    ax.add_patch(FancyArrowPatch((xb, 3.15), (xb, 2.15), arrowstyle="-|>", mutation_scale=8, lw=0.9, color=INK))
    ax.text(xb - 0.12, 2.65, "check fails:\nsplit and resend\nthis bucket only", fontsize=6.3, ha="right",
            va="center", color=INK)
    ax.text(0.0, 1.05, "round 2", fontsize=7.5, fontweight="bold", va="center")
    hw = w * 0.55
    for x in (xb - hw - 0.05, xb + 0.05):  # the failed bucket's two halves
        bucket(ax, x, 0.45, hw, h, hatch="..", ec=INK)
        ax.text(x + hw / 2, 0.17, "\u2713", fontsize=8, ha="center", va="top", color=INK2)
    ax.text(x0 + n * (w + gap) - gap, 0.17, "\u2713 decoded", fontsize=6.3, ha="right", va="top", color=INK2)

    bx = fig.add_axes([0.62, 0.17, 0.36, 0.69])
    panel(bx, "B")
    bx.plot([p[0] for p in kp], [p[1] for p in kp], "o", color=OURS, ms=5, mec="white", mew=0.6,
            label="power-sum sketch, one party")
    bx.plot([p[0] for p in ib], [p[1] for p in ib], "s", color=OTHER, ms=5, mec="white", mew=0.6,
            label="km-peeler table, resent at double size")
    for (f, y, r), off in [(p, (0, -9)) for p in kp] + [(p, (8, 7)) for p in ib]:
        bx.annotate(f"{r} round{'s' if r > 1 else ''}", (f, y), xytext=off, textcoords="offset points",
                    ha="center", va="center", fontsize=6, color=INK2)
    bx.set_xscale("log", base=2)
    bx.set_xticks(fracs, ["1/8", "1/4", "1/2", "1"])
    bx.minorticks_off()
    bx.set_xlim(0.1, 1.25)
    bx.set_ylim(0, 45)
    bx.set_xlabel("difference estimate / true difference")
    bx.set_ylabel("bytes per difference")
    bx.grid(axis="y", color=GRID, lw=0.5)
    bx.set_axisbelow(True)
    bx.legend(loc="lower left", borderaxespad=0.1)
    save(fig, "fig4_sketch")
    return kp, ib


def main():
    c, v = fig2()
    fig1()
    s = fig3()
    kp, ib = fig4()
    print(f"fig2: classes {c}, F+R {v[0]}, F*R {v[1]}, S {v[2]}")
    print(f"fig3 A: ntHash {s['nt_e']:.2f}, every k-mer {s['every_e']:.2f}, sampler "
          + ", ".join(f"1/{1 / r:.0f} {y:.2f}" for r, y in s["samp_e"])
          + f", per kept {s['slope']:.1f}, break-even 1/{1 / s['be']:.1f}")
    print(f"fig3 B: ntHash {s['nt_k']:.2f}, universal hash every k-mer {s['uh_all']:.2f}, sampler "
          + ", ".join(f"1/{1 / r:.0f} {y:.2f}" for r, y in s["samp_k"]) + ", filter + universal hash "
          + ", ".join(f"1/{1 / r:.0f} {y:.2f}" for r, y in s["uh_k"]))
    print("fig4 B: sketch " + ", ".join(f"{f:g}: {y:.2f} ({r})" for f, y, r in kp)
          + "; table " + ", ".join(f"{f:g}: {y:.2f} ({r})" for f, y, r in ib))


if __name__ == "__main__":
    main()
