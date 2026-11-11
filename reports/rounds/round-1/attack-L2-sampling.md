# Round 1 attack, slot L2, lens: sampling

Attacker L2, sampling lens. Attacked text: `reports/rounds/round-1/ATTACKED-paper.md`.

**Prior rounds.** `scripts/round history "sampling"` and `scripts/round history "gaps"` both report no prior
mention ("a searched absence"). Round 1 is the first round on this repository, so no finding below re-files an
adjudicated point.

**Method cards.** The brief has no appendix, so I looked for the cards myself. They are in
`/home/rcb112/project/claude/bench/repos/library/methods/`. Three apply here: `pseudo_replication.md` ("your
effective sample size is the number of independent decisions, not the number of rows"), `power.md` ("before
reporting that something is missing, compute how often it would be missing anyway"), and `substitution_scoring.md`.

**Re-run.** I rebuilt `hash/stats_seeds.cpp` unchanged except for one `fprintf` that logs each long gap's record and
coordinates. The copy is `/home/agent/.claude/jobs/0d0eac3e/tmp/L2/ss.cpp`, linked against the repo's
`kp7*.o` and the pinned ntHash. I ran it on `bench/data/kmer_poly_work/cs10.fa` with k = 31, all five methods and
5 seeds. Its 25 statistics rows are **byte-identical** to `results/kp7_stats_seeds_2026-11-10.eisenstein.tsv`
(`diff` is empty), and it logged exactly 1,774 long gaps. Gap coordinates are in `.../L2/gaps.<method>.tsv`.

## F1 -- the 1,774 long gaps are 154 loci counted up to 203 times each (C13)

Paper, Results, "Sampling quality matches ntHash": "All but 2 of the 1,774 such gaps lay in soft-masked repeats,
and none of the Eisenstein hash's 367 lay outside them."

The 1,774 is the sum over 25 runs (5 methods × 5 seeds) on **one genome**. Because every hash is a function of the
k-mer, the same repeat produces a gap in run after run. Merging overlapping gap intervals across all runs gives:

- **154 distinct loci** for the 1,774 gaps, a ratio of 11.5 rows per locus. The gaps fall on 12 records.
- **One locus**, NC_044374.1:21,995,814–22,046,193 (50 kb), holds **203** of the gaps. It contributes 10–14 gaps in
  each of 16 of the 25 runs.
- 59 loci are hit by all five methods.
- Per method, gaps → loci: ntHash 341 → 102, v2gate 363 → 95, v1gate 371 → 93, Gaussian 332 → 89,
  **Eisenstein 367 → 96**.

The `pseudo_replication.md` card covers exactly this: rows (gap events) inflate n (independent loci) 11.5-fold.
**The conclusion survives the correct n.** 152 of 154 loci are ≥ 50% soft-masked, against a 47.7% base rate. The
attack the claim list suggests (that 47.7% of cs10 is soft-masked, so a gap would sit in a repeat by chance) does
NOT refute the inference. I tested that premise and it fails. The flaw is the denominator, not the direction.

**Correction.** Replace "All but 2 of the 1,774 such gaps lay in soft-masked repeats" with: "The 1,774 such gaps
(25 runs) fall at 154 distinct loci, of which 152 lie in soft-masked repeats; one 50-kb repeat on NC_044374.1
accounts for 203 of the gaps."

## F2 -- "none of the Eisenstein hash's 367 lay outside" is a zero with no power (C13)

Both gaps outside repeats belong to ntHash: seed 3 at NC_044379.1:58,815,747, length 208 with 0 lowercase bases,
and seed 4 at NC_044371.1:66,268,409, length 202. Both barely clear the 200 threshold. Those are exactly what the
chance null predicts. The script's own null is 0.14 long gaps per run (I re-derived it: about 70.3 M gaps ×
(1 − 2/21)^200 = 70.3 M × 2.03 × 10⁻⁹ = 0.14). About half the genome is uppercase, so the null predicts about
0.07 per run outside repeats, or 1.8 over 25 runs. Observed: 2.

For the Eisenstein hash alone, the chance expectation over 5 seeds is about 0.35 outside-repeat gaps, so
P(zero) ≈ e^−0.35 ≈ 0.70. Under `power.md`, the zero says nothing about the hash. The sentence puts "ntHash: 2" next
to "Eisenstein: 0" (and the results file has a column for it). That reads as a contrast the data cannot support.

**Correction.** Drop "and none of the Eisenstein hash's 367 lay outside them", or replace it with "the 2 gaps
outside repeats (both ntHash, 202 and 208 k-mers) match the about 2 expected by chance over 25 runs."

## F3 -- the bit-balance comparison sits at the genome's floor, not the hash's (C13)

"The largest per-bit deviation from one half averaged 0.0009 to 0.0010, and 0.00090 for the Eisenstein hash
against 0.00088 for ntHash."

With n = 736,547,799 k-mers, i.i.d. sampling gives each bit fraction an SD of 0.5/√n = 1.84 × 10⁻⁵, and the expected
maximum |deviation| over 64 bits is about 4.8 × 10⁻⁵ (99th percentile 7.0 × 10⁻⁵; simulated). Even counting only
the 430,060,850 distinct canonical 31-mers (`results/kp10_vs_shibuya_2026-11-11/stepA_meta_cs10_NCBI_RefSeq_softmasked.json`),
the floor is about 6 × 10⁻⁵. Every method, ntHash included, measures 0.00065–0.00132 (TSV rows), which is 10–25×
that floor.

So a factor the methods share sets the deviation: high-copy repeats, where one hash value counts thousands of times.
Hash quality does not. Within each method, the seed-to-seed range (Eisenstein 0.00078–0.00110, ntHash
0.00065–0.00102) is 15× the 0.00002 gap between the two means the sentence sets side by side. The test can only
reject a bias larger than about 10⁻³. "Matches ntHash" is an equivalence statement from an instrument that would not
distinguish hashes below that level.

**Correction.** "At this genome's repeat-limited resolution (about 10⁻³), no method's bits deviated from balance
more than ntHash's."

## F4 -- "tens of gaps per seed" is false for ntHash seed 0, and the per-seed count is mostly one locus (C13)

"Every method left tens of threshold-sampling gaps longer than 200 k-mers per seed." ntHash seed 0 left **8**
(TSV row `ntHash 0 ... 8 0.14 0 729`). Across runs the per-seed count runs 8–96, and the main driver is whether the
NC_044374.1 locus (F1) is hit, which adds 10–14 gaps at once. The whole section also rests on **one** assembly
(cs10) at one k (31). The heading "Sampling quality matches ntHash" generalizes past that sample.

**Correction.** "Every method left 8 to 96 such gaps per seed on cs10 at k = 31." Scope the heading to cs10.

## F5 -- the 1/20 break-even is a single-job point estimate, and at 1/20 itself the model says slower (C11)

From `results/paper_bench_2026-11-11/counters.tsv` (load mean 6.983 × 10⁹ cycles subtracted, ÷ 736,547,799):

- Batched ntHash: 11.365 cycles.
- Batched-sampler model: intercept 9.96, slope 28.45, giving **11.382** at r = 1/20. That is 0.15% *slower* than
  ntHash, with break-even at 1/20.3. The abstract's "at sampling rates of 1/20 and below" therefore includes a rate
  where the paper's own model predicts a loss. Results says "about 1/20", which is right.
- Pairing the two within-job repeats of each of the three inputs in all 8 combinations gives break-evens from
  **1/19.0 to 1/21.7**.
- The same results file says earlier jobs differ by "3–6% between-job variance". A 3% relative shift in the
  baseline moves break-even to 1/16.3 or 1/26.7, and 6% moves it to 1/13.6 or 1/39.2. The previous job's
  (corrected) value was 1/21.5.

**Premise check on the brief's framing.** The brief says 1/20 is interpolated "between rows using DIFFERENT
evaluation modes". That is half true. The table's sparse rows are branchy (inline), but the results file also holds
*batched* sparse rows, and they fit the batched slope: the model gives 10.245 at 1/100 (measured 10.22) and 9.99 at
1/1000 (measured 9.95). So the linear model holds within one mode, and "and below" holds wherever it was measured.
The defect is the precision of the crossing point, not the modes.

**Correction (abstract).** "...is faster than batched ntHash at sampling rates below about 1/20 (1/19 to 1/22 across
repeats in one job)."

## F6 -- the Ts/Tv range unions two strata, and the sample is one plant on one chromosome (C9)

"In pairwise alignments of cannabis assemblies on chromosome NC_044378.1, 61 to 65 percent of substitutions are
transitions (Ts/Tv 1.58 to 1.85)." From `results/ideas_eisenstein/kp1_tstv_NC_044378.1.json`:

- **All columns:** cs10–KOMPa 1.585 (61.3%), cs10–KOMPb 1.584 (61.3%), KOMPa–KOMPb 1.669 (62.5%).
- **Uppercase only:** 1.780 / 1.769 / **1.847 (64.9%)**.

The quoted range takes its minimum from the all-columns stratum and its maximum from the repeat-excluded subset of
one pair. The memory rule "an extremum is not a spread" and `substitution_scoring.md` (estimate per stratum) both
reject that. Within either stratum the range is narrower: all columns 61.3–62.5% (1.58–1.67), uppercase 63.9–64.9%
(1.77–1.85).

The sample is small, and the phrase "pairwise alignments of cannabis assemblies" hides how small. There are three
comparisons, involving cs10 and the **two haplotypes of one plant**. The two cs10 comparisons share the reference and
a plant, so they are not independent. The third is not a pairwise alignment at all: it is "compared only at cs10
positions both align (projected through cs10 coordinates)" (`ideas/eisenstein/kp1_tstv.py` docstring). All three
are synNet-filtered.

**Correction.** "In synNet alignments of the two haplotypes of one plant (KOMP) against cs10 on NC_044378.1, 61 to 63
percent of substitutions are transitions (Ts/Tv 1.58 to 1.67), 64 to 65 percent outside soft-masked repeats."

## F7 -- "equal reliability" equates 30/30 trials on 6% of the data with 1 run per field (C15)

Abstract and Discussion: "20.5 for an invertible Bloom lookup table of equal reliability".

- **IBLT side.** 30/30 successes at 1.05·D on partitions 0–3 only (`full_kmp_kmers_giga_r1.05.json`: `"parts":
  [0,1,2,3]`, D = 202,557, 6.25% of the 64 partitions). The exact Clopper–Pearson 95% upper bound on the failure
  rate after 30/30 is **11.6%**, so the data allow up to about one failure in nine.
- **The same sizing on the other pair** failed 4 of 30 (`full_kmp_kmers_draft_r1.05.json`, failure_rate 0.133). The
  paper reports that, but "equal reliability" in the abstract is stated without the pair qualifier.
- **KP-10 side.** One run per field (`full_kp10_eis_giga.json`, `full_kp10_gau_giga.json`). The two runs share the
  same bucket assignment (identical `splits: 1275`), so they are not two independent draws.
- **The regimes differ.** KP-10's 11.27 bytes include a second, interactive round (4,561,762 entries = 2¹⁶ × 67 +
  1,275 × 2 × 67; the round-2 bytes ARE counted, which I checked). The IBLT figure is one-shot.
- **The configuration was tuned on the test pair.** KP-10's (2¹⁶, m = 64) is the best of four configurations
  calibrated on partition 0 of this same GigaDB pair (`kp10_*_b1{6..9}_m*.json`), and then reported on the whole
  pair, partition 0 included.

The direction survives. Even the IBLT's 50%-failure point at n = D costs 19.6 B/diff, above 11.3. "Equal
reliability" is not established by these sample sizes, though.

**Correction.** "...against 20.5 for an invertible Bloom lookup table sized to succeed in 30 of 30 trials on a
sixteenth of the pair."

## Claim-list steering note

The C13 entry steers toward the base-rate attack ("47.7% of cs10 is soft-masked"), which fails when tested (F1). It
does not mention the unit of count, which is the real defect: 1,774 rows from 154 loci. The C11 entry's
"different evaluation modes" premise is likewise half wrong (F5). Neither looks deliberate, but both point an
attacker at the weaker objection.

## Verdicts

| claim | verdict | finding |
|---|---|---|
| C1 | UNCHALLENGED | not reached under this lens |
| C2 | UNCHALLENGED | |
| C3 | UNCHALLENGED | |
| C4 | UNCHALLENGED | |
| C5 | UNCHALLENGED | |
| C6 | UNCHALLENGED | |
| C7 | UNCHALLENGED | |
| C8 | UNCHALLENGED | |
| C9 | NARROWED | F6 |
| C10 | SURVIVED | counters reproduce 11.37 and 24.69 (load-subtracted); 2-repeat spread ≤ 1.3% verified; ratios at 2 s.f. are inside the 3–6% between-job variance |
| C11 | NARROWED | F5 |
| C12 | UNCHALLENGED | |
| C13 | NARROWED | F1, F2, F3, F4 |
| C14 | UNCHALLENGED | |
| C15 | NARROWED | F7 |
| C16 | UNCHALLENGED | |
| C17 | UNCHALLENGED | |
| C18 | UNCHALLENGED | |
| C19 | UNCHALLENGED | |
| C20 | UNCHALLENGED | |
| C21 | UNCHALLENGED | |

C1 UNCHALLENGED · C2 UNCHALLENGED · C3 UNCHALLENGED · C4 UNCHALLENGED · C5 UNCHALLENGED · C6 UNCHALLENGED · C7 UNCHALLENGED · C8 UNCHALLENGED · C9 NARROWED · C10 SURVIVED · C11 NARROWED · C12 UNCHALLENGED · C13 NARROWED · C14 UNCHALLENGED · C15 NARROWED · C16 UNCHALLENGED · C17 UNCHALLENGED · C18 UNCHALLENGED · C19 UNCHALLENGED · C20 UNCHALLENGED · C21 UNCHALLENGED
