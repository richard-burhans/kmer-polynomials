# Round 2 adjudication of L2 (sampling)

Adjudicator: skeptic for L2. The attack file's sha256 is `3f8cc8e39824e031...`. It matches the seal in
`round.json:172`, so I proceeded.

**How I worked.** I read the attack first and the claims second. I re-ran every premise I could from the committed raw
outputs (`results/kp10_vs_shibuya_2026-11-11/laila_b69c5c76/`, abbreviated `L/`), and did one light local scan. All
local work ran single-process under `nice -n 19 ionice -c3`:
- I recomputed the runs test from `run.log` `@@ iblt_trial`.
- I recomputed the Poisson overflow model.
- I read the km-peeler source at the pinned commit b8ec39b (`run.log:2` confirms the commit).
- I computed per-record lengths of the two GigaDB-pair FASTAs in `bench/data/kmer_poly_work/` with awk. This is a
  linear read with no hashing.

**Prior rounds.** I ran `scripts/round history` (adjudications) on "sawtooth", "power of two", "independent random
seeds" and "comes from the second round". None has a prior ruling.
- Two prior rulings bear on this attack:
  - **Round 1 adjudication-L2 F7** ruled that KP-10's (2^16, m = 64) configuration was calibrated on this pair.
  - **Round 1 adjudication-completeness Gap 5** ruled on the 61%-against-3% point. It NARROWED the point and endorsed
    reporting both figures.
- F1 discloses and extends the first ruling. F6 partly re-files the second without citing it; see F6.

**Method cards.** The cards are at `../library/methods/`, library commit 2601816.
- `pseudo_replication.md:8` says "effective sample size is the number of independent decisions". It governs F1 and
  F3, and I apply it below.
- No card covers Clopper–Pearson under dependent trials. I say so rather than invent one.

---

## F1: 11.3 against 20.5 is one point on a sawtooth. UPHELD on the defect; fix NARROWED

**Premise re-run.**
- **The sizing rule.** `kp10c.c:541` reads `while (... dest / 2^bits > 0.8 * m) bits++`.
- **Bytes per bucket.** `kp10c.c:574` gives (m + c + 1) × 8 = 536 B.
- **Round 1 alone is closed-form.**
  - GigaDB: 536 / (3,238,101 / 65,536) = **10.85**. `L/tables/kp10_seeds.tsv` round1_bytes_per_difference is
    **10.85**.
  - Draft: 536 / (7,931 / 256) = **17.29**. The measured value is **17.30**, one round, 0 splits.
- **Round 2.** My Poisson recomputation at λ = 49.41 gives 1,253.4 expected overflowing buckets. The measured mean is
  1,254 (`kp10_seeds.tsv`).
- **What the bytes depend on.** They are a deterministic function of D's position on the power-of-two grid, plus
  ±0.18% seed noise. I also confirmed the cliff: `kp10_dest.tsv:7` gives 21.70 at 2^17.
- **The draft pair shows the sawtooth too.** At D_est = D/2 it used 2 rounds and cost **14.87**, less than its 17.30
  at the true D (`kp10_dest.tsv:10-11`). That is a direct measurement of the effect F1 describes, and the attacker did
  not cite it.
- **The admissible rule.** The rule constant 0.8 must lie in [0.772, 0.968) to reproduce both 2^16 (GigaDB) and 2^8
  (draft) from the true D. 0.8 sits near the lower edge. The attacker's 0.75 counterexample is correct: it gives 2^17
  and 21.7. The `kp10c.c:30-31` comment does say the rule was fixed to reproduce the earlier runs' bucket counts.

**Scope.**
- The defect is real, and under the sampling lens it is the right one. The paper's unit of evidence for "how big
  against the table" is D, and n = 2.
- The paper already says "These results come from one genome and two pairs". That covers genome breadth, not D's
  position in the octave.
- The Discussion's causal sentence ("Because an overfull bucket is resent on its own, the sketch spent 11.3 ... where a
  one-round table spent 20.5") is where the point-estimate does the most harm.

**The prescribed fix.**
- **Endorsed, and free:** a Results sentence disclosing three things:
  - the sizing rule's history (it reproduces bucket counts calibrated on this pair; round 1 F7);
  - that a 4% overestimate moves GigaDB to 2^17 at the measured 21.70;
  - that at the same overestimate the table is about 21.35 (modelled), so it is a near tie, not a reversal.
- **NARROWED: the abstract wording.** "Because buckets come in powers of two, the sketch's cost ... ranges from about
  11 to 21" has three problems:
  - **The range is simulated.** It must be labelled modelled. Only the 11.3 and 17.3 endpoints are measured.
  - **It is conditional on this rule.** The range belongs to this rule (m = 64, 0.8 m, power-of-two counts), not to
    the sketch. Buckets are ranges of S and could be cut at any count, so "because buckets come in powers of two"
    states an implementation choice as if it were intrinsic.
  - Endorsed instead, and free: the abstract gets a short qualifier, e.g. "at a difference near the top of its bucket
    range". The modelled 11–21 range goes in Results, labelled as from the sizing rule.
- **Verified for the whole set.** The draft row (C6) must get the same treatment. Its 17.3 is load 31.0, which the
  attacker correctly noted.

## F2: "its smaller size comes from the second round" confounds the round with the bucket grid. UPHELD on the defect; fix NARROWED

**Premise re-run.** The Poisson numbers reproduce exactly:

| m | expected overflowing buckets | P(any) | B/diff |
|---|---|---|---|
| 84 | 0.176 | 0.161 | 14.09 |
| 88 | 0.017 | 0.017 | 14.73 |
| 90 | 0.005 | 0.005 | 15.06 |

- `kp10c` accepts `--m` up to MAXM (`kp10c.c:515`).
- With m = 88 the paper's own 0.8 m rule keeps 2^16 buckets (49.41 < 70.4). So "sized for one round" by raising m is
  the natural reading of the paper's own design, not a contrivance.
- The "sized for one round" row (`kp10_dest.tsv:7`) is in fact D_est = 2D. It is a doubled estimate, not a one-round
  design.
- The attribution "comes from the second round" is therefore unsupported. The model says most of the gap to the table
  is available in one round.

**The fix.**
- **Endorsed, and free:**
  - Delete "so its smaller size comes from the second round".
  - Relabel the 21.7 row and ledger rows 153 and 177 as "sized for one round by doubling the buckets (m = 64)".
  - Apply the same deletion to the Discussion's "Because an overfull bucket is resent on its own" attribution for the
    true-D comparison. It stays valid for the D/8 point (C2).
- **OVERTURNED: publishing the modelled split.** The attacker's replacement sentence would install "14.7" and "3.4 of
  9.2" as paper numbers on a model alone. Three reasons:
  - The model is calibrated (1,253 against 1,254), but it is unrun at m = 88.
  - The one-round m = 88 design succeeds with P ≈ 0.983. That is not the 30-of-30 standard the table was held to.
  - Its own split sits at the same favourable point of the octave (F1).
- **The owed run is endorsed.** It costs one laila job: `kp10c --m 88 --dest 3238101` on giga, 30 OS seeds. The
  number enters the paper only after that run.

## F3: km-peeler's 30 seeds give overlapping hash triples. UPHELD on the fact; NARROWED in consequence

**Premise re-run.**
- **The seeding.** `km-peeler/ibflib.c:27-29` (b8ec39b) has `j = root_seed + i`. `ibflib.c:331` places partition j at
  `hash % chunk_size + j * chunk_size`, and the partitions are equal-sized.
  - So trial t + 1's hypergraph reuses two of trial t's three hash functions, with identical within-partition
    collision patterns.
  - `kp10_vs_iblt_job.sh:176` passes `1000 + t`.
- **The clustering.** My runs-test recomputation on draft n = D gives 6 runs against 12.7 expected, P = 0.002. That
  survives a correction for the attacker's two tests. Giga at n = D gives P = 0.66.
- **The KP-10 side.** Its seeds are OS-random.

**Scope.**
- **C7.** "Independent random seeds" is false for km-peeler: its seeds are consecutive and not random. **REFUTED as
  worded.**
- **C5's conclusion survives, and is strengthened.** "The runs do not compare the two methods' reliability" still
  holds. Positive dependence widens the table's true bound beyond 9.5%; it does not narrow it. What fails is the 9.5%
  figure attached to the table, in ledger rows 152 and 179 and in "Each 30 of 30". **The attacker's "C5 REFUTED" is too
  strong.** C5 is NARROWED: restrict the 9.5% to the sketch.
- **What the attacker did not credit.** Selection is symmetric in kind:
  - the table's 1.05 was chosen on these seeds;
  - KP-10's m = 64 and its rule constant were chosen on this pair (round 1 F7; F1 above).
  - The prose should not imply that only the table was tuned on the data.

**The fixes.**
- **Endorsed, and free:** the C7 and C5 wording and ledger rows 152 and 179, as the attacker gives them.
- **Endorsed, and costs one laila run:** rerun `iblt_sweep` with root seeds spaced at least r = 3 apart, or drawn at
  random. One amendment: confirm the 1.05 sizing on the **fresh** seeds as a held-out check, since 1.05 was selected
  on 1001–1030. Spacing alone does not undo that selection.

## F4: Figure 4B's grid hits the table's cheapest phase, and lines join step functions. UPHELD on the defect; one leg of the fix NARROWED

**Premise re-run.**
- **The grid.** `make_figures.py:250` plots `fracs = (0.125, 0.25, 0.5, 1.0)` with joined lines (`"-o"` and `"--s"`)
  and `set_xlim(0.1, 1.5)`.
- **The table's resends.** `iblt_protocol.tsv` totals equal 1.05 D × (1 + 1/2 + ...). For example, at D/8:
  19.55 × 1.05 × 1.875 = 38.49, and the measured value is 38.49. So the final resend lands at exactly 1.05 D, the
  cheapest passing size.
- **The off-grid costs.**
  - At 0.9 D: (0.945 + 1.89) × 19.55 = 55.4.
  - At D/6: 2.625 × 19.55 = 51.3. The sketch at D/6 gets 2^14 buckets, the same as at D/4.
  - Both are arithmetic from measured inputs, and both are modelled.

**The fix.**
- **Endorsed, and free:**
  - the caption sentence that the estimates are powers of two, where the table's last resend is exactly 1.05 D;
  - markers only, or steps. This is a figure re-render with no re-measurement.
- **NARROWED: "show the sketch's 2D point".** Under rule 3, adding the sketch's 2D point alone shows the overestimate
  side for one method only. `iblt_protocol.tsv` has no giga point above D. At D_est = 2D the table would be about
  2.1 × 19.55 = 41 B/diff (modelled), against the sketch's measured 21.70. So the omitted region is mostly *favourable*
  to the sketch. Only the narrow window 1.036–1.06 D favours the table.
  - **Endorsed, free:** state in the caption or Results that a 4% overestimate moves the sketch to 21.7, about the
    table's cost at that estimate.
  - **Costs a run:** plotting a 2D point needs both methods. Only the sketch's exists.
- **Also,** the x-limit of 1.5 with no data above 1 should be trimmed to 1. That is free.

## F5: the 20,000-k-mer check vouches for neither k = 63 nor its hash path. UPHELD

**Premise re-run.**
- **The 20,000-k-mer check.** `recon/kp7e_hash.py:85-97` draws uniform random codes at `K = 31`
  (`prototype/kmers.py:9`) and compares them with `kp7e_direct`.
- **The k = 63 values.** They come from the C streamer (`aswritten.py:251`, `_stream_into`).
- **The truth.** It is KMC strings hashed by the C reference (`aswritten.py:224-245`).
- **What the attacker missed.** There is an additional k = 63-specific check: `aswritten.py:286-297` (c5) compares the
  streamer's k = 63 values with the C reference on genome substrings. It strengthens the attacker's point that the
  real guarantee lies elsewhere.

**The fix.** Endorsed as written, and free. The Methods sentence becomes "hashed with the C stream hasher, decoded in
a numpy prototype ... recovered values equal the C-reference hashes of KMC's canonical 63-mer difference, element by
element". The Results sentence (line 150) is already correct.

## F6: "61 percent against about 3 percent" compares different populations. NARROWED

**Disclosure.** This partly re-files round 1's Gap 5 without citing it. Gap 5 ruled that the discrepancy is real, and
endorsed reporting the regime explanation that the results file gives (`prior.md:106-107`). The new content is the
identity of the authors' datasets: SARS-CoV-2, plus random sequences at pm 0.001 and 0.01. I verified it at
`shibuya_2022_genomic-reconciliation.md:747-782`. The ratios are 3.4%, 3.5% and 2.3%. That part is a valid, free
addition.

**Two legs of the correction are OVERTURNED, because each installs a stronger claim:**
- **"Their method is exact only when k-mers are unique."** The authors (lines 382-383) state a *sufficient* condition:
  "In the extreme case where each k-mer occurs once, our method computes exact k-mer set differences". "Only when"
  turns it into a necessary one, which they do not claim.
- **"This difference consists of scaffolds that share sequence with the chromosomes".** It is offered as the
  explanation ("the population, not the harness, explains it"). The sharing is unmeasured. No file counts how many
  scaffold k-mers occur in chromosome extended syncmers. The direction is predicted by the authors' mechanism; the 61%
  magnitude is not. So the harness-error question stays open, as the paper already says ("which we did not trace").

**Endorsed, free:** "against about 3 percent in the authors' tests on viral and random sequences, which are nearly
free of repeats. Their method is exact when every k-mer occurs once."

## F7: "almost all from unplaced scaffolds" has no measurement behind it. NARROWED

**What the attacker should have checked.** It is correct that no file counts the 3,237,326 A∖B k-mers by record. But
"no measurement behind it" overstates the case. Two quick measurements bound it tightly:
- **The other direction is tiny.** `L/tables/truth.tsv:2` gives B∖A = **775**. Any chromosome edit that removes A
  k-mers by substitution also creates B k-mers. So such edits are bounded at roughly the 10^3–10^4 scale, under 1% of
  3.24 M.
- **The chromosomes match exactly in length.** My scan of the two FASTAs shows the ten chromosomes have **identical
  lengths** in both assemblies (CM028014–23 against Cs_Cb2_01–10, matched one for one). Cannbio-2's 137 other records
  sum to **10,317,051** bases, exactly the difference in `sets.tsv`. So there is no net indel in the chromosomes.

**The residual loophole.** Base-to-N masking in the release would remove A k-mers without creating B ones. Only a
per-record count closes it.

**The fix.**
- **Endorsed, and free:** cite the evidence in one clause: "the ten chromosomes have the same lengths in both, and
  only 775 k-mers are new in the release".
- **Owed, and cheap on laila:** the exact per-record count. One kmtools pass over 137 records.
- **Overturned as an alternative:** the attacker's "consistent with" wording. It is weaker than the evidence warrants.

## On the steering note and correct minorities

- **The steering note is fair on both counts.** C1's refutation conditions omit D's position, and C7 pre-asserted
  independence.
- **F3's seed overlap was found only by this attacker.** I record it as a correct finding from a single critic, upheld
  on the code.
- **F1's sawtooth** is likewise this lens's alone among the files I was permitted to read.

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | bytes are a closed-form function of D's octave position (10.85 and 17.29 reproduce exactly; draft at D/2 is cheaper than at D); fix NARROWED: the abstract range is modelled and rule-specific, not "because powers of two" |
| F2 | UPHELD | "comes from the second round" rests on a doubled-bucket row; deleting it is free; the modelled 14.7 / 3.4-of-9.2 must not enter the paper until the owed m = 88 run (one laila job) |
| F3 | UPHELD | km-peeler seeds are root + i with roots 1001–1030 (ibflib.c:27-29), runs test P = 0.002 reproduces; C7 wording false, but C5's conclusion is strengthened, not refuted; rerun with spaced or random seeds plus a held-out check of 1.05 costs one laila run |
| F4 | NARROWED | grid-phase and joined-lines defects are real and the caption and marker fixes are free; showing only the sketch's 2D point is one-sided (the table at 2D is about 41 B/diff, modelled) |
| F5 | UPHELD | the 20,000 check is k = 31 numpy against C on random codes; the k = 63 guarantee is the C-reference check against KMC (plus c5); fix free |
| F6 | NARROWED | the dataset identity is a valid free addition; "exact only when" and the "population, not harness" explanation install stronger, unmeasured claims and are overturned |
| F7 | NARROWED | no per-record count, but B∖A = 775 and identical chromosome lengths support "almost all"; cite them (free) and owe the exact count; the "consistent with" wording is overturned |
