# Round 1 adjudication of L3 (statistics)

Attack read: `reports/rounds/round-1/attack-L3-statistics.md`, sha256 prefix **df828a763b2c**, which matches the sealed
digest. I read the attack first and the claims second, as the brief requires. Scratch work is in
`/home/agent/.claude/jobs/0d0eac3e/tmp/ADJ-L3/` (`pb.py` recomputes the per-repeat benchmark numbers from
`counters.tsv` and `run.log`).

Prior rulings: `scripts/round history "collision bound"` and `"overflow"` with `--kind adjudication` return a searched
absence. Round 1 is the first round in this repository, so no finding re-files a settled point.

Method cards: these live at `../library/methods/` (there is no `refs/methods/` in this repo). `power.md` governs F5.2
and `pseudo_replication.md` governs F5.1. No card governs benchmark design (F3, F4) or the algebraic arguments
(F1, F2, F7), and I say so where it matters.

---

## F1: 2.7e-17 is 2k/q, not 2(k−1)/q; the F_{p²} figure is wrong (C1)

**Premise re-run.** With q = 2^61−1, 60/q = 2.602e-17 and 62/q = 2.689e-17. With p = 2^31−1, 60/p² = 1.301e-17 and
62/p² = 1.344e-17. The total degree of S(s)−S(u) is at most 2(k−1): FR has degree 2(k−1), λ(F+R) has degree k, and
λ² has degree 2. The Methods sentence derives 2(k−1)/q and then quotes 62/q. The Implementation section's 1.4e-17 is
not what either form gives. The Abstract's "about 2k/q, 2.7 × 10⁻¹⁷" is internally consistent, and 2k/q is a valid
(looser) upper bound, so the Abstract is not wrong.

**Fix.** Endorsed and free. Use "about 2.6 × 10⁻¹⁷" where the paper says 2(k−1)/q, and "about 1.3 × 10⁻¹⁷" for
F_{p²}. The Abstract can stay as written. **Whole set (rule 3):** the same figures also appear in
`manuscript/ledger.md` lines 16–17 and in `results/kp10_vs_shibuya_2026-11-11.md` lines 43 and 129. Correct the
quantity on every surface, not just the paper.

**UPHELD.**

## F2: a per-pair bound is not "separates"/"injective" at genome scale (C1, C6)

**Premise re-run.** The distinct counts are 442,151,894 (Cannbio-2) and 430,060,850 (cs10), from `stepA_meta_*.json`.
C(n,2)·60/q = 2.54 and 2.41. An ideal 61-bit hash gives 0.042 and 0.040, and an ideal 64-bit hash gives 0.0053 and
0.0050. The arithmetic is correct.

**Diagnosis, split by leg:**
- **"Injective"** appears only in the claim list (C1's wording), not in the paper. C1's own refutation condition names
  this exact risk ("if 'injective' overstates a probabilistic statement"). A random evaluation into F_q is not
  injective. What is injective is the map from classes to *formal polynomials* S. This leg is UPHELD. The fix is
  free: the ledger and claim wording should say "injective as polynomials; per-pair collision ≤ 2(k−1)/q".
- **"Separates" in the paper.** Both the Abstract and the Methods follow "separates distinct double-stranded k-mers"
  immediately with the probabilistic bound (Abstract: "...k-mers: any two collide with probability at most...").
  Read in context, "separates" refers to the polynomial and the probability is stated in the same sentence. It is
  loose, not false. The attacker's replacement wording is accurate and free, and I endorse it as an improvement, not
  as a correction of a falsehood.
- **"~2.5 expected colliding pairs per assembly".** This is correct only as what the *bound* permits. It is an upper
  bound on the expectation, not an estimate of it, and the attacker words it that way in the body ("the proven bound
  permits"). The comparison with ideal hashes compares a worst-case proven bound with an idealised expectation. It
  shows the bound is loose, not that the hash is worse. Adding one sentence on the bound's genome-scale meaning is
  endorsed and free, **provided** it says "the bound allows at most about 2.5 colliding pairs in expectation", not
  "the hash produces about 2.5 collisions".
- **The k = 7/8 checks.** The attacker concedes they validly test the non-identity step, which is all the paper uses
  them for ("Exhaustive enumeration agrees"). No defect.

**NARROWED.** The defect is the claim list's "injective" plus loose wording. The paper states the probabilistic bound
next to "separates", so nothing in the paper is false.

## F3: "faster than batched ntHash at 1/20" (C11)

**Premise re-run** (`pb.py`). With each repeat's own load subtracted, I reproduce the attacker's numbers exactly.
Repeat 1: ntHash 11.293, collected slope 28.62, intercept 9.87, predicted cost at 1/20 = 11.300. Repeat 2: ntHash
11.440, slope 28.37, intercept 10.05, predicted 11.469. With the file's mean-load subtraction, repeat 1 still gives
11.333 against 11.327. So 1/20 does sit marginally on the slow side, in every variant.

**But the scale of the "falsity" is 0.006–0.03 cycles (≤ 0.3%).** That is below the resolution the attacker argues
for in F4, and the ledger and Results already say "about 1/20 (1/20.3)". "Refuted" overstates a rounding issue. The
real defects are:
- The **Abstract** drops the "about" that the Results carries. The fix is free: restore it, or write "below about
  1/20".
- **"A kept k-mer costs about 28.5 cycles"** is the *collected-evaluation* slope (1/16 and 1/4 rows only). The
  inline rows give 64 cycles from the 1/16–1/4 rows ((25.39−13.32)/0.1875) and 68 from the attacker's 1/100–1/1000
  rows, with a break-even near 1/31–1/33. The caption already says which rows are inline, but the 28.5 sentence does
  not. The fix is free: "about 28.5 cycles per kept k-mer with collected evaluation".
- **"And below" survives.** Collected evaluation is linear in rate, and both measured sparse points (inline) beat
  ntHash (1.14×, 1.21×). The best mode at each rate is faster below the break-even.

**Repair: OVERTURNED as written.** The attacker's prescribed parenthetical, "interpolated between the measured 1/16,
0.97×, and 1/100, 1.14×", misdescribes the derivation. The attacker's own table shows the break-even comes from the
collected 1/4 and 1/16 rows and is a short *extrapolation* beyond 1/16, not an interpolation to 1/100. Interpolating
between 1/16 (collected) and 1/100 (inline) is exactly the cross-mode slope that C11's refutation condition warns
against, so the fix would install the very error the claim flagged. Endorsed wording, also free: "below about 1/20,
extrapolated from the collected-evaluation slope (28.5 cycles per kept k-mer, 1/4 and 1/16 rows); measured 0.97× at
1/16 and 1.14× at 1/100".

**NARROWED.** The diagnosis holds as a precision and labelling defect, not as a refutation. The prescribed repair is
overturned.

## F4: timing provenance and precision (C10, C14)

1. **n = 2.** The arithmetic is right (6.35·|d| half-width gives about ±4.5% from a 0.7% spread). But a t-interval on
   two points is a weak instrument. The relevant fact is that every ratio in the paper is *within* one job, and the
   in-process wall-clock loop timings agree far better than the counters do (ntHash 2.743 s and 2.742 s across the
   two repeats, 0.04%; `run.log`). The attacker concedes that "matches AVX-512" is stated correctly. Four-figure
   table entries with a stated spread are conventional. Endorsed (free): one clause noting that between-job variance
   is 3–6%, so ratios are good to about two figures. **Narrowed.**
2. **Fixed order within repeats.** True (`run.log`), and it is a design weakness. However, the in-loop wall-clock times
   reproduce across positions to under 1% for every method I checked, so there is no evidence of position drift.
   This is a note, not a defect. **Narrowed to a disclosure.**
3. **Load subtraction "assumes no overlap, never tested".** **OVERTURNED.** `bench_e.cpp` loads the whole FASTA (line
   122) before any timed loop, so load and hashing are sequential phases of one process and additivity holds by
   construction. Independent check: the in-process loop timers, which never include the load, give
   ntHash/e-L2 = 2.743/5.941 = **0.462** and ntHash/branchy-1/1000 = 2.743/2.260 = **1.214**. These match the
   counter-derived 0.46× and 1.21×. The attacker's 0.61× and 1.10× "end-to-end" ratios fold FASTA parsing into a
   hash comparison, which is not what the paper compares. The attacker should have checked the per-round wall-clock
   lines that sit beside every counter run in `run.log`. Stating the load size (about 9.5 cycles per k-mer) is a
   harmless optional disclosure, not owed.
4. **Spaced-seed timing.** **UPHELD.** run_51 is a single process timed by its own wall clock (3 internal rounds), not
   counters, and it was not repeated. That contradicts Methods' "every timing ... in its own counted process ...
   twice". Per-hash rate 21.55/28.79 = 0.749 against time ratio 0.759, so 0.75× against 0.76×. Which denominator is
   fairer depends on why SeedNtHash emits 1.4% more windows (N handling). Either way the difference is 0.01, so I do
   not require the number to change. The fix is free: say in Methods or Results that the spaced-seed rate is one
   wall-clock run in the same job, and optionally give 0.75× per hash.

**NARROWED.** Leg 4 is upheld, legs 1–2 are reduced to disclosures, and leg 3 is overturned.

## F5: long-gap evidence, power, uncalibrated tests (C13)

1. **Pseudo-replication.** UPHELD on diagnosis under `pseudo_replication.md`. In
   `kp7_stats_seeds_2026-11-10.eisenstein.tsv` the max gap is about 14,290 in exactly 15 of 25 rows and 713–788
   otherwise. Per-seed counts are bimodal and overdispersed (ntHash 8, 87, 56, 94, 96). The 1,774 gaps are not 1,774
   independent observations. **Scope narrowed:** the *conclusion* ("long gaps reflect repeats") survives any
   plausible effective n. Against a 52.3% unmasked base rate, 2 outside out of even a few dozen distinct loci is
   decisive, and every method/seed exceeds the 0.14 chance expectation by more than 50×. Repair: a free wording fix
   ("the same repeat loci recur across seeds and methods") is endorsed. Reporting per distinct locus would cost a
   re-run, because gap coordinates are not in the TSV. That is desirable but not owed for this sentence.
2. **Power of "none of the Eisenstein's 367".** UPHELD under `power.md`, which requires the expected count and P(zero)
   beside any reported zero. 367 × 2/1774 = 0.41 and P(0) = 0.66. These are recomputed and correct, with the caveat
   that the pooled rate is itself pseudo-replicated, so they are illustrative. The fix is free: drop the clause or add
   "(0.41 expected at the pooled rate; P(zero) = 0.66)".
3. **Uncalibrated bit and density tests.** The i.i.d. SD of 1.84e-5 is recomputed. The observed 0.0009 is about
   15–20× the i.i.d. expectation for the maximum, for every method, so the input's repeat structure sets the floor.
   The "effective n ≈ 2 × 10⁶" is back-solved, not measured, so I give it no weight. The paper's heading, "matches
   ntHash", is not false: "indistinguishable at the input's resolution" is what it means. Endorsed (free): add that
   resolution caveat. A random-oracle baseline would cost a run and is not required.

**NARROWED.** Leg 2 is upheld, leg 1 is upheld on diagnosis but its conclusion survives, and leg 3 is a free caveat.

## F6: "of equal reliability" (C15)

**Premise re-run** (JSONs). Sketch: `rounds_max 2`, `splits 1275`, sketch_elements 4,561,762 = 65,536·67 + 1,275·2·67,
and 36,494,096 B / 3,238,101 = 11.27 B/diff, so the bytes include the resends (the attacker agrees). IBLT at 1.05×:
30/30 on parts 0–3 (D = 202,557, 6.3%), 20.53 B/diff. At 1.0×: 15/30, 19.55 B/diff. Draft pair at 1.05×: 26/30.
Clopper–Pearson upper bound for 0/30 = 1 − 0.025^(1/30) = 11.6%. All of this is correct.

- **Rounds.** Methods *does* disclose the mechanism ("It is then split and resent, at the cost of one more round").
  What is missing is that *this run* used it (1,275 buckets) while the IBLT row is one-shot. This leg is narrowed to
  a disclosure.
- **"Equal reliability."** UPHELD. 30/30 on 6% of the difference bounds the failure rate only below 11.6%. The paper's
  own draft-pair run failed 4/30 at the same ratio, two sentences later. The sketch's "every run" is deterministic
  decoding plus probabilistic overflow detection (see F7). Calling the two "equal" asserts a measured parity that
  was never measured. The direction survives: even the IBLT at n = D (19.6 B/diff, which fails half the time) is
  larger than 11.3.
- **Repair.** Endorsed and free: drop "of equal reliability" from the Abstract and Discussion, and drop "for the same
  reliability" in Results. The attacker missed this third surface, the Results sentence "The k-mer table needed 1.8
  times as many bytes for the same reliability"; per rule 3 it must change too. Replace with the attacker's
  descriptive wording ("sized to succeed in 30 of 30 trials on 4 of 64 partitions; our sketch used two rounds, 1,275
  buckets resent").

**UPHELD.** This holds on the phrase, with the round-count leg reduced to a disclosure because Methods already
describes the mechanism.

## F7: overflow detection is probabilistic, not deterministic (C16)

**Premise re-run.** I read `recon/recon61.py` `decode` (lines 100–129). Acceptance requires (a) Padé solvability,
(b) the flagged own elements reproducing only S_{m+1}..S_{m+c}, and (c) flagged counts differing by δ. S_1..S_m are
matched implicitly when P and Q split completely over flagged elements, so the attacker's "never checks S_1..S_m" is
not an independent gap. The core point stands, though. A bucket with > m differences is wrongly accepted iff a
wrong signed set of ≤ m own elements matches the true one in m + c power sums. The difference of the two sets can
have more than m + c + 1 distinct elements, so Vandermonde does not exclude the coincidence. Excluding it rests on
the randomness of x and λ. That makes detection probabilistic over the hash, which is exactly C16's refutation
condition. The `minsky_2003` summary (lines 23–25) confirms that unknown-bound verification is probabilistic in the
source the paper cites.

**Repair: NARROWED.** Drop "deterministically" in the Discussion and "rather than resolved wrongly" in Methods; this
is free. **I do not endorse the attacker's "except with probability about q^{-c}".** That figure is a heuristic. No
derivation is given, and the event also needs the Padé numerator to split over the own set, so the true bound is
neither stated nor proven. Writing q^{-c} would install an unproven number where an unproven adverb used to be.
Endorsed wording: "detected except with negligible probability over the random choice of x and λ". A stated bound
can be added later only if one is derived. Same correction applies to `results/kp10_vs_shibuya_2026-11-11.md`
section C / table D.

**UPHELD** on the defect, **OVERTURNED** on the q^{-c} repair. The governing verdict is UPHELD.

## F8: Ts/Tv range mixes a subset row (C9)

**Premise re-run** (`results/ideas_eisenstein_2026-11-10.md` lines 21–24). 1.847 / 64.9% is "same, uppercase only",
a subset of KOMPa vs KOMPb (1.669 / 62.5%). The paper's "(Ts/Tv 1.58 to 1.85)" presents it as one of the alignments.
The conversions are correct. The "two haplotypes of one plant = nearly one observation" leg does not matter for a
descriptive range, which is not an inference, so I give it no weight. **Fix endorsed and free:** "61 to 63 percent
(Ts/Tv 1.58 to 1.67), 65 percent outside soft-masked repeats".

**UPHELD** (the independence leg is immaterial).

## F9: the "bit for bit" baseline check is partial and not in the timed job (C12)

**Premise re-run.** `bench_e.cpp` line 162 `if (only) return 0;` precedes its identity pairs, and those pairs omit
NT_BATCH. `bench_pack.cpp` lines 274–294 hold the only check: XOR + count over all k-mers, plus element-by-element on
the first 2,000,000 k-mers of record 0. No committed log under `results/` records that output. Every hit for
"bit-identical" is prose. **I re-ran the check** with the local `hash/bench_pack` on `cs10_head6m.fa`: "ntHash batched
(ntbatch) identical to ntHash canonical (fwd+rev): yes" and "positions and hashes identical on the first 2000000
k-mers of NC_044371.1: yes", exit 0. The claim is very likely true. The local binary is not verified to be the timed
job's build (f56af99), and this was 6 Mb, not all of cs10. "About twice": 23.53/11.37 = 2.07.

**Fix endorsed.** The wording change is free. Committing a log of the check on whole cs10 costs one run of
`bench_pack` (minutes) and is owed before "bit for bit" can stand unqualified.

**UPHELD.** The evidence is partial and uncommitted, though my re-run supports the claim's direction.

---

## Minority and corroboration

I can only speak to this lens. Two of my rulings depend on a lone observation:
- **F7's repair.** The attacker's q^{-c} figure is overturned. The diagnosis is upheld.
- **F4.3.** This leg is overturned by evidence the attacker did not consult: the in-loop wall-clock timers in `run.log`.

If other lenses filed F7 or "equal reliability", agreement among same-family agents is not corroboration. The
corroboration here is the decoder code and the JSONs I re-read.

## Verdict table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | 2(k−1)/q = 2.6e-17 and 2(k−1)/p² = 1.3e-17; Methods and Implementation quote 2.7 and 1.4. The fix is free; also correct ledger.md:16–17 and kp10_vs_shibuya.md:43,129 |
| F2 | NARROWED | "injective" (claim list only) must go and the 2.5-pairs arithmetic is right as an upper bound; the paper already states the bound next to "separates"; the wording fixes are free |
| F3 | NARROWED | 1/20 is on the slow side by ≤ 0.3%, a precision and labelling defect; the 28.5 needs "collected"; the prescribed "interpolated 1/16–1/100" repair is OVERTURNED as a cross-mode misdescription |
| F4 | NARROWED | the spaced-seed single wall-clock run is upheld (free disclosure); the load-subtraction "overlap" leg is OVERTURNED (sequential phases; loop timers reproduce 0.46× and 1.21×); n = 2 and fixed order reduce to disclosures |
| F5 | NARROWED | the zero-with-no-power leg is upheld (0.41 expected, P(0) = 0.66, per power.md, free); pseudo-replication is real but the repeat conclusion survives; the bit-test calibration needs a free caveat |
| F6 | UPHELD | "equal reliability" is unmeasured (30/30 on 6%, CP ≤ 11.6%; 26/30 on the draft pair); drop it on all three surfaces, including Results' "for the same reliability"; the round mechanism is already disclosed in Methods |
| F7 | UPHELD | wrong acceptance is excluded only by the hash's randomness, so drop "deterministically" (free); the q^{-c} figure in the repair is OVERTURNED as unproven |
| F8 | UPHELD | 1.85 is a soft-mask-filtered subset row; "61–63% (1.58–1.67), 65% outside repeats" is free |
| F9 | UPHELD | the identity check is partial (XOR + first 2M) and unlogged, and it did not run in the timed job; my 6 Mb re-run passes; committing a whole-cs10 log costs one run |
