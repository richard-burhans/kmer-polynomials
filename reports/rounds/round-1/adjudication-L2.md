# Round 1 adjudication of L2 (sampling lens)

Adjudicator for L2. The attack file `attack-L2-sampling.md` has sha256 `7192e64acbe22e29...`. That matches the sealed
digest `7192e64acbe2` in `round.json`, so I adjudicated the sealed attack.

**Prior rulings.** `scripts/round history "sampling" --kind adjudication` and `... "gaps" ...` report no prior mention.
This is round 1 on this repository, so nothing here re-files a settled point.

**Method cards.** Three cards govern here: `library/methods/pseudo_replication.md` (F1), `power.md` (F2) and
`substitution_scoring.md` (F6). The memory rule "an extremum is not a spread" also applies to F6. F3, F5 and F7 have no
governing card. There I ruled on the arithmetic.

**What I re-ran.**
- The attacker's instrumented copy `L2/ss.cpp` differs from `hash/stats_seeds.cpp` only by one logging `fprintf` and
  its three globals (`diff` checked). Its 25 statistics rows, sorted, are identical to
  `results/kp7_stats_seeds_2026-11-10.eisenstein.tsv`. The attacker's "byte-identical" holds only after sorting,
  because the row order differs. The gap log holds 1,774 rows, and every per-(method, seed) count matches the TSV's
  `long_gaps`.
- I independently merged overlapping gap intervals across all 25 runs.
- From `counters.tsv` and `run.log` I refit the sampler model per repeat.
- I read `kp1_tstv_NC_044378.1.json` and the `kp10_vs_shibuya_2026-11-11/*.json` files directly.

---

## F1: 1,774 gaps are 154 loci (C13). NARROWED

**Premise verified.** My own merge gives **154 loci on 12 records**. The largest locus is
NC_044374.1:21,995,814–22,046,193 with **203 gaps**, hit in **16 runs** at 10–14 gaps each. **59 loci** are hit by all
five methods. At **152 loci** every gap is at least half lowercase. The two exceptions are the two ntHash gaps.

One detail does not reproduce. My per-method locus counts (ntHash 98, v2gate 91, v1gate 89, Gaussian 81,
Eisenstein 91) differ from the attacker's (102/95/93/89/96). This is probably a merge-convention difference
(per-method merge versus loci of the global merge). It does not bear on the finding.

**Scope.** The published sentence ("All but 2 of the 1,774 such gaps lay in soft-masked repeats") is literally true as
a count of events. It is not an error of fact. The defect is that 1,774 is presented as the weight of evidence when
the independent unit is the locus. That is exactly the `pseudo_replication.md` case: "the number of independent
decisions, not the number of rows". The attacker says correctly that the conclusion survives at n = 154. Its
dismissal of the base-rate objection, which it tested at 152/154 against 47.7%, is also correct.

One point strengthens the conclusion, and the attacker underplays it. The same loci recur across five unrelated hash
functions (59 of 154 hit by all five). That recurrence is direct evidence that the sequence, not the hash, makes the
gaps, and a corrected sentence should keep it.

**Repair.** The diagnosis holds and the repair is ENDORSED, with one addition: "(59 of them hit by all five methods)".
**The repair is not free.** The numbers exist only in the attacker's scratch copy (`/home/agent/.claude/jobs/.../L2/`).
Before the paper cites them, the gap logging and the locus merge must be committed as a versioned script and re-run
from the repo. Cost: a small code change plus one stats_seeds run (about the cost of the original run).

## F2: the Eisenstein zero has no power (C13). UPHELD

**Re-derived.** (1 − 2/21)^200 = 2.03 × 10⁻⁹. Times 70.3 M gaps that is 0.1425 per run. The uppercase share is 52.3%,
so for Eisenstein over 5 seeds λ = 0.373 and P(0) = e^(−0.373) = **0.69**. Over all 25 runs about 1.9 are expected
outside repeats, and 2 were observed (both ntHash, lengths 202 and 208, with 21 and 0 lowercase bases). Under
`power.md` a zero that is the modal outcome under the null "says nothing about your subject". The clause "none of the
Eisenstein hash's 367 lay outside them" reads as a contrast with ntHash's 2 that the data cannot support.

**Repair.** ENDORSED, either variant. Deleting the clause is **free**. The replacement ("the 2 gaps outside repeats,
both ntHash, match the about 2 expected by chance over 25 runs") is also free, because it uses only TSV values and
the script's own null.

## F3: bit-balance resolution (C13). NARROWED on the diagnosis; OVERTURNED on the repair

**Diagnosis.** The i.i.d. floor reproduces: σ = 1.84 × 10⁻⁵. The simulated expected maximum over 64 bits is
4.8 × 10⁻⁵, with a 99th percentile of 7.0 × 10⁻⁵. At 430 M distinct k-mers it is 6.3 × 10⁻⁵. Observed values are
0.00065–0.00132. The direct evidence that the instrument resolves only about 10⁻³ is the seed-to-seed range within
each method: Eisenstein 0.00078–0.00110, ntHash 0.00065–0.00102. That range is about 15 times the 0.00002 difference
between the means the sentence sets side by side. This part holds.

The attacker's causal leg is inferred, not measured. "High-copy repeats, where one hash value counts thousands of
times" is plausible, since identical k-mers share a hash and multiplicity inflates the variance. But no
multiplicity-weighted variance was computed. So the defect is narrowed to "the comparison sits at a resolution of
about 10⁻³". The mechanism stays unproven.

**Repair OVERTURNED: it installs a false claim.** "No method's bits deviated from balance more than ntHash's" is false.
From the TSV, the means are ntHash **0.00088**, Eisenstein 0.00090, v1gate 0.00101, fourth-root provable 0.00103 and
v2gate 0.00103. All four other methods have a higher mean than ntHash. The single largest value, 0.00132, belongs to
the fourth-root provable hash. The attacker should have checked its replacement against the same TSV rows it cited.

A correct and **free** repair: "The largest per-bit deviation from one half averaged 0.0009 to 0.0010 for every method
(range over seeds 0.00065 to 0.00132), so this test resolves only deviations near 10⁻³, and at that resolution the
Eisenstein hash and ntHash are indistinguishable."

## F4: "tens of gaps per seed" (C13). NARROWED

**The count leg is UPHELD.** ntHash seed 0 left **8** long gaps (TSV), and 8 is not "tens". The range is 8–96. The
repair "8 to 96 such gaps per seed" is **free** and correct.

**The driver leg is overstated.** The NC_044374.1 locus contributes 10–14 gaps in 16 runs (verified). But runs without
it have 8–59 gaps and runs with it have 66–96. The separation (about 25–30) exceeds the locus's own contribution, so
co-occurring large repeats (the 14 kb `max_gap` rows) also drive the count. "Mostly one locus" is not established.

**The scope leg is weak.** The paragraph already opens "Across the whole cs10 reference", so the body is scoped. Only
the heading generalizes. Scoping the heading to cs10 is harmless and free, but it is optional.

## F5: the 1/20 break-even (C11). NARROWED

**Premise verified.** With load (6.98 × 10⁹ cycles) subtracted, batched ntHash costs 11.33 and 11.40 cycles per k-mer
in its two repeats (mean 11.37). The table's own line reads "beats batched ntHash (11.37) below 1/20.3", which is a
two-point fit on 1/16 and 1/4. I refit it three ways:

| fit points | break-even | 8-combination range |
|---|---|---|
| 1/16, 1/4 | 1/20.33 | 1/19.1 to 1/21.7 |
| 1/100, 1/16, 1/4 | 1/20.13 | 1/19.2 to 1/21.1 |
| all four | 1/20.02 | 1/19.1 to 1/21.0 |

Exactly at r = 1/20 the model predicts 11.37–11.39 against ntHash's 11.37. That is a tie, not a win. The attacker is
right that "1/20 and below" in the abstract claims a point that is a coin flip, and that the Results' "about 1/20" is
the defensible form.

The attacker's premise check of the brief also holds. The batched 1/100 and 1/1000 rows (10.22 and 9.95) sit on the
batched line, so "and below" holds wherever it was measured. The defect is cosmetic precision at the crossing point.
This is a narrowing of the brief's own steer, which was "different evaluation modes".

**Repair.** ENDORSED in the short form "below about 1/20", which is **free**. The parenthetical "(1/19 to 1/22 across
repeats in one job)" is accurate but optional for an abstract. The extrapolation to 3–6% between-job shifts (1/13.6 to
1/39.2) is arithmetic on a variance quoted from earlier jobs. It is not a measurement, so it should not enter the
text.

## F6: Ts/Tv range unions two strata (C9). UPHELD

**Verified from the JSON.** All columns: 1.585, 1.584 and 1.669 (61.3%, 61.3%, 62.5%). Uppercase only: 1.780, 1.769
and 1.847 (64.0%, 63.9%, 64.9%). The published "61 to 65 percent (Ts/Tv 1.58 to 1.85)" takes its minimum from one
stratum and its maximum from another. That is the extremum-not-a-spread error, and `substitution_scoring.md` (per
stratum) governs it.

The description is also verified. The `kp1_tstv.py` docstring confirms synNet MAFs, cs10 against KOMPa and KOMPb ("the
two haplotypes of one plant"). It also confirms that the third pair is "compared only at cs10 positions where both
align (projected through cs10 coordinates)", so it is not a pairwise alignment.

What survives is the qualitative statement "a transition, the commoner substitution", which holds in every stratum
and every pair.

**Repair.** ENDORSED and **free**, with one refinement: give Ts/Tv for the uppercase stratum too. Suggested text: "61 to
63 percent of substitutions are transitions (Ts/Tv 1.58 to 1.67), 64 to 65 percent outside soft-masked repeats (1.77
to 1.85)".

## F7: "equal reliability" (C15). UPHELD on the defect; the repair is NARROWED (incomplete)

**Verified.**
- `full_kmp_kmers_giga_r1.05.json` has `parts [0,1,2,3]`, D = 202,557 and 30/30 successes.
- The Clopper–Pearson 95% upper bound on the failure rate is **11.6%** (1 − 0.025^(1/30)).
- The draft pair failed 4 of 30 at the same 1.05 sizing.
- KP-10 ran once per field. `full_kp10_eis_giga.json` has `rounds_max: 2` and `sketch_elements` 4,561,762 =
  2¹⁶·67 + 1,275·2·67, so the second-round bytes are counted.
- The configuration (b16, m64) is the cheapest of the four `kp10_eisenstein_b1*_m*` calibrations, at 11.32 B/diff.
  Those calibrations were run on part of this same pair.

The two sides differ in kind. KP-10's reliability is exactness by detect-and-resend, which costs an extra round and
depends on C16, attacked under another lens. The IBLT's reliability is a one-shot success probability bounded only to
11.6% or less. "Equal reliability" equates the two, and the data do not establish it. The direction survives: even
the IBLT at n = D (15/30) costs more than 11.3 B/diff.

**Repair narrowed.** The attacker's replacement fixes the abstract only. Rule 3 applies: the same equation appears at
three more sites. These are "the same reliability" (Results, line 128), "0.55 times the bytes of an invertible Bloom
lookup table of equal reliability" (Discussion, line 136) and the abstract's "20.5 ... of equal reliability". All
must change together. The repair should also say that the sketch used a second round on this pair, so that "one-shot
versus two rounds" is not hidden. With those additions the repair is **free**, because it is text only and every
number is already in the JSON.

---

## Minority and steering notes

- F1's and F5's challenges to the claim-list steering are correct and recorded. The brief's suggested attacks (the
  47.7% base rate for C13, and "different evaluation modes" for C11) both fail when tested. This lone attacker found
  the real defects behind them.
- F3 is the round's instance of rule 2: the diagnosis is right and the prescribed repair is a stronger false claim.

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | 154 loci / 152 in repeats / 203-gap locus reproduce; the sentence is true as an event count, and the defect is the unit; repair endorsed (add "59 hit by all five methods") but not free: commit the gap-logging and merge script and re-run |
| F2 | UPHELD | P(zero Eisenstein outside-repeat gaps) ≈ 0.69 under the script's own null; dropping the clause is free |
| F3 | NARROWED | resolution of about 10⁻³ holds (the seed range is 15 times the mean gap), the repeat mechanism is unmeasured, and the repair "no method deviated more than ntHash" is FALSE (four methods have higher means); a free correct wording is given |
| F4 | NARROWED | "tens per seed" is false (ntHash seed 0 = 8) and "8 to 96" is a free fix; "mostly one locus" is overstated and the body is already scoped to cs10 |
| F5 | NARROWED | break-even is 1/20.0 to 1/20.3 (repeats 1/19 to 1/22), so "1/20 and below" claims a tie; "below about 1/20" is free; the between-job extrapolation should not enter the text |
| F6 | UPHELD | the range unions the all-columns minimum with the uppercase maximum; the per-stratum wording is free |
| F7 | UPHELD | 30/30 on 4 of 64 parts bounds IBLT failure only to 11.6% or less against a two-round deterministic sketch; the repair must cover all three "equal/same reliability" sites and the second round, and is then free |
