# Round 2 adjudication of L3 (statistics)

Adjudicator: skeptic for L3. Attack read: `attack-L3-statistics.md`, sha256 `c4f6e91c0e76...`, which matches the
seal in `round.json:173`. Read in the harness order: attack first, then claims C1-C15 and the frozen text
(`ATTACKED-paper.md` lines 14, 63, 67, 135-156).

**What I re-ran.** All checks were pure-Python arithmetic on committed numbers, under `nice -n 19 ionice -c3`. No
genome was read and nothing ran on laila.

- **Tables:** `kp10_seeds.tsv`, `kp10_dest.tsv`, `iblt_trials.tsv`, `iblt_protocol.tsv`, `truth.tsv`.
- **`run.log` greps:** the seed fields, and `n_per_round` for every protocol repetition.
- **Code:** `recon/kp10c.c`, in particular the sizing loop at :541, `reconcile` at :319-397 and `decode` at :257-266.
  Also `recon/aswritten.py:136-142` and `recon/kp10_vs_iblt_job.sh:114-204`.
- **Earlier results:** `results/round1_owed_2026-11-11.md:119-122` and `results/kp10_vs_shibuya_2026-11-10_prior.md:84-96`.

**Prior rounds.** I ran `scripts/round history --kind adjudication` for "power of two", "bits_for", "phase", "second
round" and "1001".

- Round 1 L2 F7 (UPHELD) required disclosing that the sketch used a second round. The text now does, and nothing here
  re-files that.
- No earlier adjudication covers the power-of-two quantization, the protocol lattice, the m/2 sizing rule or the IBLT
  seed scheme. All five findings are new.

**Premise of the overflow model (F1, F2, F5), tested.**

- **Capacity.** `decode` accepts a bucket iff its symmetric difference is at most `mbar`, which is m or m-1 with the
  parity of delta. Because d has the parity of delta, a bucket overflows exactly when it holds more than 64
  differences.
- **Cost of an overfull bucket.** Each one costs two buckets in the next round (`kp10c.c:388-392`). That is the
  attacker's `1 + 2P` form.
- **The giga pair.** I reproduce every number the attacker gives: 1,253.4 expected splits against 1,253.97 measured,
  11.263 B/diff modelled against 11.26 measured, and 1.6 x 10^-6 expected overfull buckets at 2^17.
- **A second load the attacker did not use.** At the draft pair with D_est = D/2 the load is 62.0 per bucket
  (bits 7). The model expects 46.9 splits and 14.99 B/diff; the run measured 46 splits and 14.87 (`kp10_dest.tsv:10`).
- **Three measured loads.** The model is now checked against measurements at loads of 31, 49.4 and 62 (17.30, 11.26
  and 14.87 B/diff). A one-round sketch costs exactly 536 B per bucket, so the low end of the sawtooth is also pinned
  by measurement: 21.70 at a load of 24.7.

The model is sound. Its use in F1 and F2 is mostly sound, but the prescribed fixes overreach, as set out below.

**Method cards.** Two cards govern, both from `../library/methods/`.

- `pseudo_replication.md`: n is the number of independent decisions. This decides F1.
- `power.md`: state P(the outcome) before reporting it. This decides F5.

No card governs F3's protocol lattice. That ruling rests on arithmetic from the code alone.

---

## F1: the 11.3 vs 20.5 sits at the sketch's cheapest phase (C1, C13)

**Diagnosis: UPHELD.**

- **Where D falls.** It is 3,238,101 = 0.965 x 51.2 x 2^16, so the first round runs at a load of 49.4 of the 51.2 the
  0.8m rule allows.
- **The boundary.** Any D_est above 3,355,443 (1.036 D) selects 2^17 buckets, and the bytes are then exactly
  131,072 x 536 = 70,254,592, the 21.70 already in Table 4. This leg is measured, not modelled: a one-round sketch's
  bytes are fixed by its bucket count, and the expected number of overfull buckets there is 1.6 x 10^-6.
- **The unit of analysis.** Under `pseudo_replication.md` the decision that sets the size ratio is the pair and its D.
  The 30 seeds replicate the hash, with a byte spread of 0.36 percent. "In 30 of 30 runs" is a statement about
  recovery, and the seeds do replicate recovery. So the abstract does not misstate what was scored. It does present one
  favourable phase as the comparison.

**Scope: narrower than filed.**

- **What the attacker calls owed is largely measured already.** The "sawtooth outside D" is called a model projection
  owing a run. Its ends are measured (11.26 at a load of 49.4, 21.70 at 24.7), and its interior is checked at loads of
  31 and 62 (17.30 and 14.87).
- **The only free choice is the average.** The one unmeasured quantity is the "log-uniform average 15.2", and that
  depends on a prior over D that nothing in the paper supplies.

**Fix: NARROWED.**

- **Endorsed, and the fix is free.** Say in the abstract or in Results, beside the 11.3, that with power-of-two bucket
  counts this difference fell just under a boundary. An estimate 4 percent higher would have doubled the first round to
  2^17 buckets, the 21.7 the paper already reports.
- **Endorsed in Results only, also free.** Add the measured range across loads: 17.3 at 31 per bucket, 14.9 at 62 and
  11.3 at 49, with one round costing 21.7 at 24.7.
- **OVERTURNED:** "about 15 on average" and the "phase-averaged 1.35x as the headline". Both are prior-dependent model
  outputs. Putting either in the abstract would install an unmeasured number in place of a measured one.

## F2: "its smaller size comes from the second round" confounds two causes (C4, C13)

**Diagnosis: UPHELD on the defect.**

- **The comparator.** The single-round comparator is the next power of two, at an expected 1.6 x 10^-6 overfull
  buckets. It is not the smallest one-round sketch.
- **The paper's own counter-example.** The paper reports 17.30 B/diff in one round on the draft pair, at a load of 31
  per bucket, with 0 splits in 30 of 30 runs (`kp10_seeds.tsv:3`; also C6). That is below the table's 20.5 with no
  second round.
- **What follows.** "Comes from the second round" therefore attributes to the second round a saving that the
  power-of-two lattice shares. The same holds for the Discussion's "Because an overfull bucket is resent on its own,
  the sketch spent 11.3 ..." (`ATTACKED-paper.md:156`).

**OVERTURNED on the quantitative decomposition.**

- **The reliability match.** The "matched reliability" rows match the sketch's P(any overflow) to the table's 9.5
  percent. That is the table's one-sided 95 percent UPPER BOUND, not its failure rate. The table went 30 of 30 at
  1.05 D after 18 of 30 at 1.0 D, so its true rate could be far below 9.5 percent.
- **The selected row.** The 1.26x / 1.44x split rests on the 0.047 row, and the "10^5 times lower" comparison has the
  same flaw. Matched to a 1 percent failure rate, the model gives 17.1, and nothing measured fixes which row is right.
- **What the attacker should have done.** Anchor the claim on the measured 17.30 at a load of 31, and bracket the
  table's failure rate instead of equating it with its bound.

**Fix: endorsed in wording only, and the fix is free.**

- "Sized for a single round at the next power of two (2^17 buckets), the sketch needed 21.7. At the draft pair's load
  of 31 per bucket it used 17.3 in one round, so its advantage here comes from the second round and the power-of-two
  bucket count together."
- Narrow line 156 to match.
- Do not publish the 1.26x / 1.44x split or the 16 to 17 B/diff one-round figures, which are unrun.
- A measured non-power-of-two one-round sketch would cost a code change and a run. It is not required for the
  narrowed wording.

## F3: the eightfold row is the table's best case (C2, C9)

**Diagnosis: UPHELD, and stronger than filed.**

- **Every lattice point.** Every D_est in the sweep is D/2^k. Doubling from 1.05 D/2^k therefore lands the table's last
  resend on exactly 1.05 D, the very ratio already tuned to succeed in 30 of 30. That holds at all four points in both
  pairs, not only at the eightfold row. `run.log` shows `n_per_round` ending at 3,400,006 or 3,400,008 for giga, and at
  8,328 for draft, in every repetition.
- **Which way it cuts.** As the attacker says, this favours the table. The abstract's 38.5 is the table's most
  favourable case, and the conclusion is understated, not reversed.

**Error in the attack's table.** The sketch column is wrong at D_est = 0.15 D and 0.20 D. Those estimates select
bits0 = 14, not 13 (485,715 / 2^13 = 59.3, which exceeds 51.2). The sketch cost there is the measured 19.40
(`kp10_dest.tsv:3`), not 20.76. The direction is unchanged.

**The 72 figure.** It is anchored, not invented. At 0.1131 D the fourth table is 0.95 D, where the measured success is
0 of 30, and a fifth table at 1.9 D certainly lists. Still, no protocol run sampled that point.

**Fix: NARROWED.**

- **Endorsed, and the fix is free.** Disclose in the Results text and the Figure 4B caption that each estimate is D/2^k,
  which puts the table's last resend at its tuned 1.05 D, its most favourable case.
- **Endorsed in Results only, as a labelled projection from the measured success rates, also free.** At an estimate
  of 0.113 D the table would need five rounds and about 72.
- **Not endorsed in the abstract:** the range "38.5 to 72". The abstract may say "38.5, the table's most favourable
  case", which is free.
- **Costs a run:** sampling off the lattice.

**Figure 4B leg.** I verified that the table's ten repetitions are byte-identical. `mean_rounds` equals `max_rounds`
and is an integer at every point, and the table sizes are deterministic given the rounds. Endorsed, and free: the
caption should say "identical in all 10 repetitions", not "whose mean is plotted". C9 otherwise survives.

## F4: the table's "independent random seeds" are fixed integers (C7; C5 survives)

**Diagnosis: UPHELD on the seed leg.**

- **The seeds.** `kp10_vs_iblt_job.sh:176` passes `1000 + t`. `run.log` contains 300 such seeds, each of 1001 to 1030
  exactly 10 times: five giga ratios and five draft ratios.
- **The sketch's seeds.** The sketch's 30 runs are OS-seeded: 60 `"seed_x": "os"` lines.
- **The verdict.** "Independent random seeds" is false as written for the table. Whether its trials are independent
  rests on km-peeler's seeded hash family. The source is not in this repo, so that check is owed. It costs reading
  code, not a run.

**Sibling surface the attacker missed.** Implementation (`ATTACKED-paper.md:67`) says "Our experiments used fixed seeds
for reproducibility". That is true of the D_est sweep (14 `0x5eed` lines) but not of the 30 sketch runs. Any rewording
of C7 must make line 67 agree.

**Tuning leg: NARROWED.**

- **Verified.** 1.05 was chosen in the 2026-11-10 run as the smallest grid ratio with 30 of 30, on partitions 0-3 of
  the same GigaDB pair (`prior.md:91`). The current 30 of 30 is therefore not a pre-fixed test.
- **The caption already says** "sized to succeed in 30 of 30", which discloses the outcome selection.
- **Which way it cuts.** Choosing the smallest passing ratio biases the table's size DOWN, toward the table, so it is
  not unfair to the table.
- **Unverified.** The contrast "the sketch's 0.8m rule had no such tuning step" was asserted without checking.
  `aswritten.py`'s docstring says earlier KP-10 runs used m/2, so the 0.8m rule was adopted at some point, and the
  attack did not establish when or why.

**Fix: endorsed, and the fix is free.**

- "each method ran 30 times on one host, the sketch with independent OS-random seeds and km-peeler with seeds 1001 to
  1030, the same at every table size; 1.05 was the smallest ratio that succeeded in 30 of 30 on four of 64 parts of
  the GigaDB pair in an earlier run"
- Reconcile line 67 with that sentence.

**C5.** I recomputed 1 - 0.05^(1/30) = 0.09503. It matches the claim, which survives.

## F5: the k = 63 run used the m/2 rule, not the stated 0.8m (C3)

**Diagnosis: UPHELD.**

- **The rule used.** `aswritten.py:136-142` sizes for an expected difference per bucket of at most m/2, and
  `round1_owed_2026-11-11.md:119` confirms the k = 63 run used it. Methods (`ATTACKED-paper.md:63`) states 0.8m for "the
  first round" and attaches no exception to the k = 63 sentence.
- **The numbers.** I recomputed 2^10 at a load of 20.9, 548,864 B and 25.7 B/diff. The expected number of overfull
  buckets is 9.4 x 10^-12. Under `power.md`, "one round" was certain before the run.
- **What survives.** The deterministic content of C3 survives, as the attacker says. The S-set difference equals
  KMC's canonical difference element by element, 21,382 and 7,931, against 398,217,248 forward-only.

**One overstatement in the attack.** "Under the stated rule it would be 274,432 B, 12.8 B/diff" omits that at 0.27
expected overfull buckets about 24 percent of runs would take a second round. The 12.8 is the round-1 floor, not the
expected cost.

**Fix: NARROWED.**

- **Endorsed, and the fix is free.** Add the clause "buckets sized from the true difference at m/2 = 32 differences
  each, twice the count the 0.8m rule gives", and keep the 549 KB out of comparison with Table 4.
- **The alternative costs a run:** rerunning at 2^9 is a whole-genome k = 63 job on laila. It is not needed for the
  claim.

---

**Correct minority.** I read no other round-2 attack, as the brief requires, so I cannot count agreement. The
power-of-two phase dependence (F1, F2) and the protocol lattice (F3) appear in no earlier adjudication. This record
upholds both diagnoses on the strength of one attacker's model, which I re-ran and checked against a load the
attacker did not use.

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | Diagnosis upheld: D sits at 0.965 of a power-of-two boundary, and a 3.6% overestimate gives the measured 21.7. Free disclosure endorsed; the "about 15 average" and "1.35x headline" fixes are overturned as prior-dependent model numbers |
| F2 | NARROWED | Attribution to the second round alone is upheld as a defect (measured 17.3 in one round at a load of 31). The 1.26x/1.44x decomposition is overturned: it matches the sketch to the table's 9.5% upper bound, not its rate. Free wording fix endorsed |
| F3 | UPHELD | Every D/2^k estimate lands the table's last resend on its tuned 1.05 D, so 38.5 is its best case and the claim is understated. The attack's sketch column is wrong at 0.15 and 0.20 D (19.40, bits 14). Free caption, Results and "identical in all 10" fixes endorsed; abstract range not endorsed; the off-lattice sweep costs a run |
| F4 | NARROWED | The table's seeds are fixed 1001-1030 reused at every ratio (verified), and 1.05 was chosen on parts 0-3 of the same pair. The selection favours the table, and "0.8m had no tuning" is unverified. Free C7 wording and line-67 reconciliation endorsed; the km-peeler seed-independence check is owed (code read, no run); C5 survives |
| F5 | UPHELD | The k = 63 run used the m/2 rule (9e-12 expected overfull buckets), not the stated 0.8m, so its "one round" had no power. The free disclosure is endorsed; the rerun at 2^9 costs a run and is not required |
