# Round 1 attack, slot L3, lens: statistics

Attacked text: `reports/rounds/round-1/ATTACKED-paper.md` (draft d049b77). Every number below was re-derived from
the raw artifact named beside it. Method cards: the brief has no appendix, so I looked for them myself and found them
at `../library/methods/`. The cards that govern findings here are `pseudo_replication.md` (F5) and `power.md` (F5, F6).
Prior rounds: this is round 1 of this repository and no closed round has ruled on any of these claims, so I have no
adjudication to disclose.

## F1 -- the quoted 2.7e-17 is 2k/q, not the 2(k-1)/q the sentence derives, and the F_{p^2} figure is wrong too (C1)

The paper says they "collide with probability at most $2(k-1)/q$ ... At $k = 31$ that is about $2.7 \times 10^{-17}$".
With q = 2^61 - 1: 2(k-1)/q = 60/q = **2.60e-17**. The quoted 2.7e-17 is 62/q = 2k/q, a looser bound the paragraph
never derives. The same slip is in the Implementation section: "its collision bound is about $1.4 \times 10^{-17}$" for
F_{p^2}, p = 2^31 - 1. There 2(k-1)/p^2 = **1.30e-17**, and even 2k/p^2 = 1.34e-17 rounds to 1.3, not 1.4. The source of
both is `results/kp10_vs_shibuya_2026-11-11.md` ("Gaussian about 2k/p² = 1.4 × 10⁻¹⁷ ... Eisenstein about 2k/q = 2.7 ×
10⁻¹⁷"), where 2k/p^2 is miscomputed as well.
**Correction:** "at most 2(k−1)/q, about 2.6 × 10⁻¹⁷ at k = 31" (Methods and Abstract) and "about 1.3 × 10⁻¹⁷" for
F_{p^2}. If the abstract keeps "about 2k/q", it should say that 2k/q is a rounded-up form of 2(k−1)/q.

## F2 -- a per-pair bound of 2.6e-17 is not "separates" or "injective" at genome scale: it allows ~2.5 colliding pairs per assembly (C1, C6)

The paper says the factorization "separates distinct double-stranded k-mers", and the claim list puts it as "injective".
What is proven is a per-pair probability. Applied to the inputs the paper itself uses (distinct canonical 31-mers from
`results/kp10_vs_shibuya_2026-11-11/stepA_meta_*.json`):

| set | distinct classes | C(n,2) · 2(k−1)/q | ideal 61-bit hash | ideal 64-bit hash |
|---|---|---|---|---|
| Cannbio-2 | 442,151,894 | **2.54** | 0.042 | 0.0053 |
| cs10 | 430,060,850 | **2.41** | 0.040 | 0.0050 |

So the proven bound permits an **expected ~2.5 colliding pairs within one cannabis assembly's k-mer set**. That is 60×
the expectation of an ideal 61-bit hash and 480× that of the 64-bit hash the paper calls "measured, not proven". The
bound is honest but weak: "provable" here means "at most a few collisions per genome in expectation", not "none".
The exhaustive k = 7/8 checks (8,192 and 32,896 classes, ≤ 5.4e8 pairs, so ≤ 1.4e-8 expected collisions) cannot
detect anything at this scale. They test only the algebraic non-identity of S(s) − S(u). That is a valid test of the
factorization step, but it is not evidence about collision rates.
**Correction:** replace "separates distinct double-stranded k-mers" with "distinct double-stranded k-mers give distinct
polynomials, so any two collide with probability at most 2(k−1)/q", and add one sentence that over a 4.4 × 10⁸-k-mer
set the bound allows about 2.5 colliding pairs in expectation. The claim list's "injective" must go.

## F3 -- "faster than batched ntHash at 1/20" is false at 1/20 in both repeats; the break-even is an interpolation from one evaluation mode (C11)

`results/paper_bench_2026-11-11.md` gives the break-even as **1/20.3**: slope 28.5 and intercept 9.96, from the
*collected* ("batched") 1/16 and 1/4 rows only. A break-even of 1/20.3 means that **at 1/20 exactly the sampler is
predicted slower**: 9.96 + 28.5/20 = 11.385 against 11.37. I recomputed per repeat from `counters.tsv`, subtracting
each repeat's own load run:

| repeat | ntHash batched | intercept | slope | break-even | predicted at 1/20 |
|---|---|---|---|---|---|
| 1 (runs 2, 16, 17) | 11.293 | 9.868 | 28.63 | 1/20.1 | **11.300 > 11.293** |
| 2 (runs 26, 40, 41) | 11.440 | 10.050 | 28.37 | 1/20.4 | **11.469 > 11.440** |

In both repeats, 1/20 sits on the slow side. The margin is about 0.01 to 0.03 cycles, far inside the 0.7% (0.08
cycle) spread of the ntHash row itself. So 1/20 is a tie at best, and it is never measured. The table measures 1/16
(0.97×) and 1/100 (1.14×).

Second, "a kept k-mer costs about 28.5 cycles" holds for collected evaluation only. The inline (branchy) rows give
(10.01 − 9.40)/(0.01 − 0.001) = **68 cycles per kept k-mer**, with break-even (11.37 − 9.33)/68 ≈ 1/33. The table's
sparse rows are inline and its dense rows are collected, so the 28.5 is not the cost of the rows printed at 1/100 and
1/1000. "And below" survives only if the faster mode is chosen at each rate: collected is monotone, and both measured
sparse points beat ntHash.

The project's own between-job variance is 3–6% ("All within the 3–6% between-job variance", same file). It moved this
break-even from 1/21.5 (corrected) to 1/20.3 between two jobs.
**Correction (Abstract and Results):** "faster than batched ntHash at sampling rates below about 1/20 (interpolated
between the measured 1/16, 0.97×, and 1/100, 1.14×)". Also state that 28.5 cycles per kept k-mer is for collected
evaluation, and that inline evaluation costs about 68.

## F4 -- every-timing provenance and precision: two identically ordered repeats cannot carry four significant figures, and one quoted timing is neither counted nor repeated (C10, C14)

1. **n = 2 and between-job variance.** With two repeats, (max − min)/mean is a range, not an interval. A t-interval
   on n = 2 has half-width 12.71 × |d|/2 ≈ 6.4 × the difference, so the 0.7% spread on batched ntHash is a 95%
   interval of about ±4.5%. The project's own file reports 3–6% between-job variance. The table quotes 11.37, 23.53,
   24.69, 25.86 and 26.43, and the ratios 0.46, 0.44 and 0.43 to two figures, but those three ratios cannot be ordered
   at this precision. "Matches the 16-lane AVX-512 code" (24.69 against 25.86, 4.5% apart) is right *because* the
   difference is unresolved. The text should say so rather than read the ordering as a result.
2. **"Interleaved" is a fixed order.** `run.log` shows that repeats 1 (runs 1–24) and 2 (runs 25–48) run the
   methods in the *same* sequence. The repeats alternate, but the order within each repeat is never varied, so any
   position-in-sequence drift is confounded with method. (Repeat 2's batched ntHash is 1.3% slower than repeat 1's.)
3. **The load subtraction is half the measurement.** The `__load_only__` runs cost 9.44–9.51 cycles per k-mer
   (runs 1, 25, 49 in `counters.tsv`), comparable to every hash cost reported. Subtracting it assumes load and hash
   costs add with no overlap, which is never tested. End to end, batched ntHash against the every-k-mer hash is
   (9.48 + 11.37)/(9.48 + 24.69) = **0.61×**, not 0.46×, and the 1/1000 sampler is 1.10×, not 1.21×. The paper should
   give the load cost's size next to "subtracted".
4. **The spaced-seed timing breaks "every timing ... counted ... twice".** The paper says timings come "from hardware
   performance counters" and that each method ran "twice, in interleaved order". The 0.76× spaced-seed figure is run_51:
   **one** run of `bench_spaced_e`, timed by its own wall clock, not by counters. Its denominator is also unequal.
   SeedNtHash emitted 42,667,269 hashes against 42,088,079 for ours. 0.76 is the time ratio 1.482/1.953; the per-hash
   rate ratio is 21.55/28.79 = **0.75**.
**Correction:** say that the spaced-seed rate is a single wall-clock run (0.75× per hash). Report cycle costs with an
honest uncertainty (± a few percent, between jobs), and give the load cost.

## F5 -- the long-gap evidence is pseudo-replicated, and the Eisenstein "none outside" is a zero with no power (C13)

1. **Pseudo-replication (`methods/pseudo_replication.md`: "effective sample size is the number of independent
   decisions, not the number of rows").** The 1,774 gaps are 5 methods × 5 seeds run over the *same* genome, and the
   same loci recur. In `results/kp7_stats_seeds_2026-11-10.eisenstein.tsv` the longest gap takes one of two values in
   every row: either about 14,290 k-mers (14,275 to 15,354, in 15 of 25 rows) or about 730 (713 to 788). That is one
   locus, present or absent depending on whether its repeat-unit k-mer is sampled. Per-seed counts are overdispersed:
   ntHash gave 8, 87, 56, 94 and 96, where Poisson with mean 68 has SD 8.3, and the seed-0 count of 8 sits 7 SD low.
   The independent unit is the repeat locus, and the paper counts each locus up to 25 times. "All but 2 of the 1,774"
   overstates the evidence by that factor. The conclusion that long gaps are repeats is plausible, but its support
   must be stated per locus.
2. **Power (`methods/power.md`: "State the expected count and P(zero) ... next to the zero").** "None of the
   Eisenstein hash's 367 lay outside them" reads as a contrast with ntHash's 2. At the pooled rate (2/1774) the
   expected count for 367 gaps is **0.41**, and **P(zero) = 0.66**. The zero is the modal outcome and says nothing
   about the Eisenstein hash. (The definition in `hash/stats_seeds.cpp` lines 16–17 and 86: a gap is "inside" when
   ≥ half its bases are soft-masked.)
3. **The bit and density tests are uncalibrated and cannot tell hashes apart.** Over 736,547,799 k-mers the i.i.d.
   SD of a bit fraction is 1.84e-5, so the expected maximum over 64 bits is about 4.4e-5. The observed 0.0009 is
   **about 20×** that, for *every* method, with overlapping seed ranges (ntHash 0.00065–0.00102, Eisenstein
   0.00078–0.00110). Minimizer-density deviations of 0.3–0.4% are similarly dozens of Poisson SEs. These deviations
   are floored by duplicated k-mers in the genome (an effective n of roughly 2 × 10⁶ independent values), not by the
   hash. No random-oracle baseline on the same input is reported, so "0.00090 against 0.00088" shows only that no hash
   defect larger than the genome's repeat floor exists.
**Correction:** "Sampling quality matches ntHash" should read "is indistinguishable from ntHash's at the resolution
the repetitive input allows". Drop "and none of the Eisenstein hash's 367 lay outside them", or give the expected
0.41 and P(zero) = 0.66 beside it. Report long gaps per distinct locus.

## F6 -- "of equal reliability" compares a two-round protocol with a one-round one, on a 30-trial subsample that the paper's own draft run contradicts (C15)

1. **Rounds.** `full_kp10_eis_giga.json`: `rounds_max 2`, `splits 1275`. The 11.27 B/diff *does* include the split
   resends: 4,561,762 entries = 65,536 × 67 + 1,275 × 2 × 67. I checked this premise, and the bytes accounting is
   consistent. But the paper's table caption does not say that the GigaDB sketch needed a second round, while the IBLT
   figure is one-shot. The IBLT could equally be sized at n = D (19.6 B/diff, 15/30) and resent on failure.
   "Equal reliability" therefore compares an adaptive protocol with a non-adaptive one.
2. **30/30 is not "reliability".** 30/30 successes leave a 95% Clopper–Pearson upper bound of **11.6%** on the failure
   rate. The trials cover partitions 0–3 only (202,557 of 3,238,101 differences, 6%), and the result is extrapolated
   to the whole genome ("66 MB vs 36.5 MB extrapolated").
3. **The same ratio fails on the paper's other pair.** At 1.05×, the draft pair failed **4 of 30** (13%;
   `full_kmp_kmers_draft_r1.05.json`: successes 26). The paper reports this itself, two sentences after asserting
   "equal reliability" at 1.05×.
4. **The sketch's "every run" is n = 1 per field per pair.** Deterministic decoding makes that adequate *if* F7 did
   not apply.
The direction (the sketch is smaller) survives: even the IBLT at n = D is 19.6 B/diff. The phrase does not.
**Correction:** "against 20.5 for an invertible Bloom lookup table sized to succeed in 30 of 30 trials on 4 of 64
partitions (one round; our sketch used two rounds, 1,275 buckets resent)". Drop "of equal reliability" in the
Abstract and the Discussion's "0.55 times the bytes ... of equal reliability".

## F7 -- overflow detection is probabilistic, not deterministic, and the cited source says so (C16)

The paper says "a bucket holding more than $m$ differences is detected rather than resolved wrongly", and the
Discussion says it "detected its own overflow deterministically". The decoder (`recon/recon61.py`, `decode`, lines
117–129) accepts a bucket when two conditions hold. First, the elements flagged on each side reproduce only the **c = 2
extra** power sums S_{m+1}, S_{m+2}. Second, the flagged counts differ by delta. It never checks the flagged elements
against S_1..S_m.

A bucket holding more than m differences is accepted wrongly whenever the Padé solution's roots among each side's own
elements happen to satisfy two field equations and the count. Nothing excludes that event. It is a coincidence over the
random hash values, heuristically of order q^{-2}, not an impossibility. By the Vandermonde argument, a signed set
whose first m + c power sums vanish needs only m + c + 1 = 67 distinct elements, and a bucket with m + 1 differences
plus an m-element wrong solution has up to 2m + 1 = 129, so the algebra allows it.

The paper's own citation states the opposite of "deterministic". `minsky_2003_set-reconciliation` (library summary,
"Unknown bound") says that confirmation at extra points is "tested probabilistically, with error per check at most
(D − 1)/|E|". The deterministic Minsky protocol is the *known-bound* one, which has no overflow check.
**Correction:** "the c extra sums check the result, so a bucket holding more than m differences is detected except
with probability about q^{-c} over the random hash" (Methods), and "detected its own overflow" without
"deterministically" (Discussion). `results/kp10_vs_shibuya_2026-11-11.md` (section C and table D) makes the same
overclaim.

## F8 -- the Ts/Tv range mixes a subset row with its parent and counts two haplotypes of one plant as separate pairs (C9)

`results/ideas_eisenstein_2026-11-10.md` lines 21–26. The upper end, 1.847 (64.9%), is "KOMPa vs KOMPb, uppercase only
(outside soft-masked repeats)". That is a filtered subset of the row above it (1.669, 62.5%), not a fourth alignment.
cs10 vs KOMPa (1.585) and cs10 vs KOMPb (1.584) compare cs10 with the two haplotypes of **one plant**, so they are
nearly one observation. Over whole alignments the range is **61.3–62.5% (Ts/Tv 1.58–1.67)**, from effectively two
independent comparisons on one chromosome. The 65% applies only to non-repeat columns of the haplotype pair. The
arithmetic of the conversion is right: 1.585/2.585 = 61.3%, 1.847/2.847 = 64.9%.
**Correction:** "61 to 63 percent of substitutions are transitions (Ts/Tv 1.58 to 1.67), 65 percent outside
soft-masked repeats".

## F9 -- the "bit for bit" baseline check is partial and did not run in the timed job (C12)

The timed job ran `bench_e` with `only=...`, and `bench_e.cpp` returns before its identity checks when `only` is set
(`if (only) return 0;`, line 162). Even unfiltered, its check pairs (lines 165–166) do **not** include
batched-vs-library ntHash. The only check is in `bench_pack.cpp` (lines 274–294). It compares the XOR of all hashes
plus the count, a linear checksum in which paired errors cancel, and it compares element by element only the
**first 2,000,000** k-mers of the first record (0.27% of 736.5M). No committed raw log under `results/` records that
check's output, and the "bit-identical" claim exists only in prose (`IDEAS.md:200`,
`kp7_next_steps_2026-11-10.md:18`). "About twice" holds: 23.53/11.37 = 2.07.
**Correction:** "matches the library's hashes element by element on the first 2 × 10⁶ k-mers and by count and XOR
checksum on all of cs10", with the log committed.

## Verdicts

| claim | verdict | finding |
|---|---|---|
| C1 | REFUTED | F1 (2.7e-17 is 2k/q; 2(k−1)/q = 2.6e-17), F2 (per-pair bound ≠ injective; ~2.5 expected pairs/assembly) |
| C2 | UNCHALLENGED | |
| C3 | UNCHALLENGED | |
| C4 | SURVIVED | 47.37/24.69 = 1.92 from `paper_bench_2026-11-11.md` counters |
| C5 | UNCHALLENGED | |
| C6 | NARROWED | F2 (the bound's genome-scale meaning is unstated) |
| C7 | UNCHALLENGED | |
| C8 | UNCHALLENGED | |
| C9 | NARROWED | F8 |
| C10 | NARROWED | F4 |
| C11 | REFUTED | F3 |
| C12 | NARROWED | F9 |
| C13 | NARROWED | F5 |
| C14 | NARROWED | F4.4 (single wall-clock run; 0.75× per hash) |
| C15 | NARROWED | F6 |
| C16 | REFUTED | F7 |
| C17 | SURVIVED | 137,216 B = 2^8 × 67 × 8, 1 round, 0 splits (`full_kp10_eis_draft.json`) |
| C18 | UNCHALLENGED | |
| C19 | UNCHALLENGED | |
| C20 | UNCHALLENGED | |
| C21 | UNCHALLENGED | |

C1 REFUTED · C2 UNCHALLENGED · C3 UNCHALLENGED · C4 SURVIVED · C5 UNCHALLENGED · C6 NARROWED · C7 UNCHALLENGED · C8 UNCHALLENGED · C9 NARROWED · C10 NARROWED · C11 REFUTED · C12 NARROWED · C13 NARROWED · C14 NARROWED · C15 NARROWED · C16 REFUTED · C17 SURVIVED · C18 UNCHALLENGED · C19 UNCHALLENGED · C20 UNCHALLENGED · C21 UNCHALLENGED
