# Round 2 attack, slot L5, lens: INFERENCE

Attacker brief: `reports/rounds/round-2/BRIEF-inference.md`. Text under attack: `ATTACKED-paper.md`, `ATTACKED-ledger.md`.
Raw outputs read: `results/kp10_vs_shibuya_2026-11-11/laila_b69c5c76/tables/{kp10_dest,kp10_seeds,iblt_trials,iblt_protocol}.tsv`,
`run.log:44,48`; `results/round1_owed_2026-11-11/run_*.json`; `recon/kp10c.c:1-31,541`;
`results/kp10_vs_shibuya_2026-11-10_prior.md:6-8,15-30,99-116`; `results/paper_bench_2026-11-12.md`.
Literature read from the library at its working tree (`../library/papers_md/`): `eppstein_2011_set-difference.md`
(full text, lines 600-620, 1210-1275), `summaries/shibuya_2022_genomic-reconciliation.md` and full text lines 770-790,
`summaries/goodrich_2011_iblt.md`. I found no held copy of a rateless IBLT paper (grep "rateless" over the library:
no hit).

**Method cards.** `../library/methods/` holds 32 cards (power, pseudo_replication, causal_designs, ...). None governs
sketch sizing or protocol comparison. `causal_designs.md` is the closest: a causal reading ("because") needs a contrast
that nets out what both arms share, which is the shape of F3. No card is violated outright.

**Prior rounds.** `scripts/round history` for "second round", "partition", "one round", "21.7", "power of two",
"underestimate", "resent whole", "must be resent", "bucketed": round 1 adjudicated the two-round vs one-round
*reliability* framing (round 1's L5 F5, L3 F6, L7 F6: NARROWED/UPHELD, now repaired in the text) and the two-party
decoder (round 1's L5 F5, UPHELD; the attacker's replacement sentence OVERTURNED). No prior round raised bucket-count
granularity, the D/8 decomposition, set-size estimation of D, or partitioned IBLTs. Where a finding touches a settled
point I say so.

**Steering.** The claim list asks whether the table's resend protocol is a straw man (C2) but frames the sketch's
sizing rule as given. The weakest inference in the section is that the sketch's bucket count is restricted to powers
of two while the table's size is continuous (F1, F2), and the list does not point at it.

Local checks: one `nice -n 19` Python process, Poisson tail sums only, no genome data touched.

## F1 -- "its smaller size comes from the second round" is mostly power-of-two rounding (C4, C6)

**Claim.** "Sized for a single round, the sketch needed 21.7, so its smaller size comes from the second round."

**Flaw.** 21.7 is not the cost of a one-round sketch. It is the cost of the next power of two. `kp10c.c:541` takes
the smallest `2^bits` with `dest/2^bits <= 0.8m`, and buckets are top-bit ranges (`kp10c.c:12`). So the only
one-round size tried was 2^17 buckets, a load of 24.7 = 0.39m. That run is the `dest = 2D` row of the sweep
(`kp10_dest.tsv:7`, bits0 17, 70,254,592 B). It was not sized for one round. Buckets are ranges of S, so any bucket
count works. The table, by contrast, was sized continuously (1.05 D).

**Demonstration.** A Poisson model of bucket loads predicts **1,253.4** overfull buckets at 2^16 (mean 49.4 per
bucket, overflow when >64). Observed: **1,253.97** (`kp10_seeds.tsv:2`, `mean_splits`). The model is therefore
calibrated on this very run. The same model gives these one-round costs at 536 B per bucket (67 entries x 8 B):

| buckets | load | expected overfull buckets | B/diff |
|---|---|---|---|
| 90,000 | 36.0 | 0.78 | 14.90 |
| 95,000 | 34.1 | 0.15 | 15.73 |
| 100,000 | 32.4 | 0.03 | **16.55** |
| 131,072 (2^17) | 24.7 | 1.6e-6 | 21.70 |

So at 100,000 buckets the chance of any overflow is about 3%, comparable to what 30/30 can show. That sketch would
finish in one round at about 16.5 B/diff against the table's 20.5. Of the 9.2 B/diff gap, roughly 4 comes from the
sketch's per-difference efficiency in one round and roughly 5 from the second round. That is not "comes from the
second round".

**Same defect in C6.** The draft pair's 17.3 B/diff is 2^8 buckets at a load of 31.0 = 0.48m. The model gives one
round at 200 buckets for 13.5 B/diff (expected failures 0.03). The paper's own sweep beats 17.3 at D_est = D/2:
**14.87** B/diff in 2 rounds (`kp10_dest.tsv:10`). C6's "17.3 against 20.5" is where 7,931 happens to fall against
a power of two, not a property of either method.

**Correction.** "Sized for a single round with power-of-two buckets, the sketch needed 21.7. A one-round sketch with
about 100,000 buckets would need about 16.5 by a Poisson model that predicts the observed 1,254 overfull buckets, so
the second round accounts for roughly half the gap." Alternatively, rerun kp10c with a non-power-of-two first round
(this is owed, and cheap on laila). Table 4's row label "sized for one round" should read "2^17 buckets (next power of
two)".

**Verdict: C4 REFUTED (attribution). C6 NARROWED.**

## F2 -- 11.3 vs 20.5 is the best point of a sawtooth, reported as the comparison (C1)

**Flaw.** Under the paper's own sizing rule, the sketch's cost depends on where D falls between two powers of two.
At 0.8m the load per bucket runs from 25.6 to 51.2 across an octave of D. The table's cost does not vary (1.05 D x
19.55 B). With the calibrated model from F1 (round 1 plus split halves for overfull buckets) across one octave:

- minimum 11.20 B/diff, at loads just under 0.8m, which is where the giga pair sits (49.4 = 0.77m);
- maximum **20.92**, at loads just over 0.4m, no better than the table's 20.5;
- mean **15.16** for D log-uniform within the octave, a ratio of 1.35x rather than 1.82x;
- the ratio is 1.82x or better in only **5.8%** of the octave.

The abstract reports the 1.82x as what the sketch does against the table. It is what the sketch does at one D, near
its luckiest load. The draft pair shows the other side: 0.48m load and 17.3 B/diff. The text never says that the
ratio varies with D.

**Also in C1 (disclosed: round 1's L5 F5 settled the decoder point).** The abstract's "it recovered every difference"
counts the union of two parties' outputs. Party A gets its 3,237,326 k-mers and never gets B's 775. The Results says
so ("Each party ... recovers only its own"); the abstract does not, even though it compares bytes against a table
whose single listing gives both sides. This is a scope issue, not a new measurement.

**Correction.** In the abstract: "... with 11.3 bytes per difference over two rounds, against 20.5 ... ; with
power-of-two buckets the sketch's cost depends on where the difference falls, and over that range it averages about 15
by a model fitted to this run." Or report both pairs side by side (11.3 and 17.3) as the range seen. Add "each party
recovers its own differing k-mers".

**Verdict: C1 NARROWED.**

## F3 -- "Because an overflowing bucket is resent alone" does not explain the D/8 result (C2)

**Claim.** "Because an overflowing bucket is resent alone, an eightfold underestimate of the difference raised this
to 20.8, against 38.5 for the table."

**Flaw.** At D/8, the sketch's first three rounds split **every** bucket. That is exactly resending the whole sketch
at double size. `kp10_dest.tsv:2` (giga, 1/8): bits0 13, rounds 5, splits **58,602 = 8,192 + 16,384 + 32,768 +
1,258**. Bytes 67,212,256 = (2^13 + 2^14 + 2^15 + 2^16) x 536 + 2 x 1,258 x 536, exactly. So the sketch paid the same
geometric series as the table, 1 + 1/2 + 1/4 + 1/8 = 1.875 times its one-round cost:

- sketch: 10.85 x 1.875 = 20.34, plus 0.42 for the final splits, gives 20.76;
- table: 20.53 x 1.875 = **38.49** (`iblt_protocol.tsv:6`, exact).

Resending "alone" mattered only in the last round, and that saving was already counted at the true D. The degradation
from the true D is the same for both: 20.76/11.27 = **1.84x** for the sketch against 38.49/20.53 = **1.87x** for the
table. The ratio of sketch to table is 0.549 at D and 0.539 at D/8, and likewise at D/4 (1.72x vs 1.75x) and D/2 (1.48x
vs 1.50x). So the 20.8 vs 38.5 gap is the true-D gap of F2 carried through an identical doubling, not something an
underestimate exposes. The results file says so itself in §4 ("1.84x its own true-D cost"; "1.87x at D/8").

On the draft pair splitting alone did help (1.23x against 1.88x, `kp10_dest.tsv:5`). It helped because the doubling
reached 128 buckets at 0.97m load, below the 0.48m sizing that the rule gives at the true D: the granularity of F1
again, on a pair the abstract does not cite.

**Correction.** "An eightfold underestimate raised this to 20.8, against 38.5 for the table: both resend at double
size until the buckets fit, so both cost about 1.85 times their own true-D size." Drop "Because an overflowing bucket
is resent alone" from the abstract, and from the caption of Figure 4 ("Resending only overfull buckets keeps the
sketch small when the difference is underestimated").

**Verdict: C2 REFUTED (the causal clause; the numbers reproduce).**

## F4 -- an eightfold underestimate cannot happen on these pairs, and the table's protocol ignores free information (C2, C13)

**Premise checked.** |A| and |B| are each party's own counts (`truth.tsv`). |A| - |B| is a lower bound on D, and it
is tight when the difference is one-sided:

- giga: 442,151,894 - 438,915,343 = **3,236,551** = 0.99952 D;
- draft: 442,159,825 - 442,151,894 = **7,931** = D exactly.

Exchanging one 8-byte set size therefore estimates D to 0.05% on both pairs the paper uses. km-peeler's cells carry
an `int64 counter` (results file §2), so the counters of a failed, subtracted table also give r(|A∖B| - |B∖A|). The
"resend at double" protocol throws that away. A table protocol that resends at 1.05(|A| - |B|) after one failure would
cost about 20.53 x (1/8 + 1) ≈ **23** B/diff at D/8, not 38.5. Caveat: 1.05 x 0.99952 = 1.0495 D, between the
measured 18/30 at 1.0 D and 30/30 at 1.05 D, so its success rate is owed, not shown. The sketch could also read its
per-bucket counts and does not, so both protocols are naive. Only the table's naivety is priced into the abstract.

This is an inference defect, not a data defect. The D/8 sweep is a synthetic misestimate for a regime these pairs, being
"each dominated by one structural difference", do not have. Two-sided differences (where |A| - |B| says little) were not
tested.

**Correction.** Either drop claim 5 from the abstract, or add: "on these pairs the set sizes alone estimate the
difference to 0.05 percent, so the underestimate is imposed, and a table resized from its own failed counters would
cost about 23." The Results should state the resend protocol as the one km-peeler's interface supports, not as the
table's only option.

**Verdict: supports C2 REFUTED; C13 NARROWED (see F5).**

## F5 -- "an IBLT that fails must be resent whole" is a protocol choice; the advantage is not the power sums' (C13)

**Claim (Introduction).** "a bucket that holds too many differences is split and resent on its own in a second round,
whereas an invertible Bloom lookup table that fails must be resent whole". Discussion: "Because an overfull bucket is
resent on its own, the sketch spent 11.3 ... where a one-round table spent 20.5."

**Flaw.**
1. Bucketing by hash range and resending only the failed buckets is a protocol over any per-bucket summary. Nothing
   about it needs power sums. This project's own earlier comparator ran km-peeler as independent IBLTs on 4 of 64
   key-space partitions (`kp10_vs_shibuya_2026-11-10_prior.md`; results §7, "partitions 0-3"). Eppstein et al.'s
   strata estimator partitions the universe into independent IBFs and decodes each separately (eppstein_2011 full
   text, lines 600-620). A partitioned IBLT that resends only failing partitions is the obvious stronger baseline, and
   it was not run.
2. The premise has a cost on the IBLT side that cuts the other way. Eppstein's Fig. 7 text (lines 1216-1222) says
   "Deltas below 200 all require at least 50% overhead". At about 50 differences per partition, a partitioned IBLT would
   pay heavily, so I do not claim it would win. The defect is the universal "must be resent whole" and the attribution of
   the 1.82x to a mechanism the competitor could share. The results file concedes this ("A smarter IBLT protocol ... is
   not what km-peeler implements", §4), but the paper turned "km-peeler does not" into "an IBLT must".
3. Rateless IBLTs extend a failed table incrementally rather than resending it (Yang, Gilad and Alizadeh, SIGCOMM
   2024, "Practical Rateless Set Reconciliation"). **I cite this from memory: it is not held in the library, and the
   prior-art lens should verify it before anyone relies on it.**

The Discussion sentence "this shows that the hash supplies the symmetry, not that it beats canonicalization" is
correctly hedged and SURVIVES.

**Correction.** Introduction: "whereas km-peeler's table, when it fails, is resent whole". Discussion: "Because our
protocol resends only overfull buckets ... ; a table partitioned the same way was not tested."

**Verdict: C13 NARROWED.**

## F6 -- "needs no canonicalization" and the niche "values that must be added" are not what the controls test (C3, C14)

1. **Linearity is the sketch's property, not S's.** Power sums are linear over any field-valued element IDs. A
   canonical 2-bit code with a universal hash into F_q gives IDs that are just as addable, and the paper measures that
   hash as 2.4x faster (Discussion). The only control run in the same sketch was a forward-only hash
   (`run_..._k31_fwd_b8_r1.json`: 256 of 256 buckets unresolved). Its failure was certain in advance because it is not
   strand-symmetric, so it cannot separate "S" from "any strand-symmetric ID". The Discussion admits as much ("not that
   it beats canonicalization") and then still lists "values that must be added" as S's niche. At k ≤ 32 that niche
   belongs equally to the canonical code. Nothing measured shows S is needed for it.
2. **Wording.** "a forward-only polynomial hash in the same sketch saw 398,217,248 differences". The sketch saw none: it
   recovered 0 and left 256 buckets unresolved (`D_recovered: 0`, `unresolved_buckets: 256`). The 398,217,248 is the
   exact set difference, computed outside the sketch (`D_true_values`).
3. **k = 63 is sized differently from the stated rule and has no comparator.** The numpy run used `bits_rule: "auto:
   bits_for(exact value D, m)"` and chose bits 10 (`run_..._k63_sym_bauto_r12.json`), a load of 20.9 = 0.33m. The C
   rule of the Methods ("smallest power of two ... at most 0.8m") gives 2^9 = 274,432 B. So 549 KB is 25.7 B/diff
   under a different, unstated rule. The Methods sentence about the first round does not describe this run. No table
   or canonical-code baseline ran at k = 63, so "where that recipe is costly" (Discussion) is unmeasured at the only
   k > 32 point.

**Correction.** Abstract: "Power sums of these hashes form a linear set sketch, and the hash's strand symmetry lets it
run on k-mers as written." Results: "a forward-only hash ... faced 398,217,248 differences and overflowed every
bucket". Methods: state that the k = 63 run sized buckets by a different rule (2^10). Discussion: drop "values that
must be added" from S's niche, or qualify it as "values that must be added, for k > 32".

**Verdict: C3 NARROWED.**

## F7 -- the 2.4x crosses two programs the paper says not to compare, and the long-k niche is unmeasured (C14)

The ledger's own row (`ATTACKED-ledger.md:57`) says 2.4x = 24.67/10.24 "across the two programs". The caption of
Figure 3 says "The panels are not compared with each other, because the two programs give the same sampler code
28.7 and 22.5 cycles per kept k-mer". The batched-ntHash baselines agree to 0.5% (11.37 vs 11.43), so the cross-program
ratio of every-k-mer paths is probably close. Still, the paper uses a comparison it forbids two paragraphs earlier,
and does not say why the forbidding reason does not apply. "Earns its place for long k-mers, where canonical codes span
several words" has no timing behind it: no two-word canonical code was benchmarked. "With a proof" for multiply-shift
on distinct canonical codes is correct (2-universal), and that part SURVIVES.

**Correction.** "about 2.4 times as fast (across two benchmark programs whose ntHash baselines agree to 0.5 percent)".
Mark the long-k niche as expected, not measured.

**Verdict: C14 NARROWED.**

## F8 -- "66 multiplications per element" is the build cost, cited as the recovery cost (C12)

The paper says: "Recovery costs more computation than a table listing. The sketch needs about m + c = 66 field
multiplications per element". The ledger sources the 66 (`ATTACKED-ledger.md:113`) to `prior.md:8`: "Its power sums
cost about m multiplications per element". That is sketch construction. km-peeler's "3 cell updates" is also
construction. The measurements point the other way on each half:

- build: kp10c 89.4 CPU-s per side against km-peeler's 240.3 (results §5, implementation-dominated);
- decode: 110.3 s against 3.5 s, a factor of **31**.

Decode costs that much because a degree-~50 polynomial is evaluated at every element of every differing bucket, and
every bucket differed: `buckets_with_difference_round1: 65536` (run.log:44). So the recovery cost per element is a
second ~m multiplications on top of the 66. The sentence states a build count, calls it recovery, and leaves out both
measured figures.

**Correction.** "Building the sketch takes m + c = 66 field multiplications per element against 3 cell updates, and
recovery evaluates a degree-m polynomial at every element of each differing bucket. Single-threaded decode took 110 s
against 3.5 s for the table's listing."

**Verdict: C12 NARROWED.**

## F9 -- "the GigaDB sketch ran before x and λ were separated" describes a superseded run (C11)

The paper's GigaDB numbers come from job b69c5c76 (kp10c). `run.log:44` (the fixed-seed 11.27 run) records
`"seed_x": "0x5eed", "seed_lam": "0x1a3d"`, two separate seeds. `run.log:48` and its 29 siblings (the 30/30) record
`"seed_x": "os", "seed_lam": "os"`. The same holds for the k = 31 draft run (seed_x 24301 = 0x5eed, seed_lam 6717 =
0x1a3d). The sentence describes the 2026-11-10 numpy run, which the paper no longer reports. It wrongly casts doubt
on the independence that the bound needs (Methods: "as long as x and λ are drawn independently").

The first sentence, "hashes every 31-mer as written ... a bucket is a range of S", matches `kp10c.c:8-12` and
SURVIVES.

**Correction.** Delete "and the GigaDB sketch ran before x and λ were separated", or say it of the earlier run that the
paper does not report.

**Verdict: C11 REFUTED (second sentence).**

## F10 -- the extended-syncmer row: wrong unit, unlike datasets, and a superseded run beside current ones (C8)

1. **Unit.** The source says "**1,986,807 spurious** k-mers" (`prior.md:105`). The paper says "spurious syncmers
   numbering 61 percent". The 507,221 extended syncmers listed are a different, smaller count.
2. **61% vs 3%, unlike regimes.** Shibuya et al.'s 3% comes from covid and random sequences with point mutations
   (pm = 0.01; full text lines 770-783, "within about 3% of the true set difference size"). Here the difference is
   whole unplaced scaffolds that, per `prior.md:106`, "share much sequence with the chromosomes". That sharing is
   exactly what produces spurious k-mers. The source names this cause. The paper sets 61 against 3 without it, which
   invites the reader to infer a defect in the method or the harness.
3. **Not like for like in Table 4.** The row is 5 trials from the superseded 2026-11-10 run (`prior.md:100`), with n
   sized from the true extended-syncmer difference. The other rows are 30 trials from b69c5c76. "2.7 times smaller"
   compares an inexact output (it missed 7) with an exact one. Round 1 (critic gap 5) left this open, and the ledger
   still marks it OPEN (`ATTACKED-ledger.md:199`).

**Correction.** "returned 1,986,807 spurious k-mers, 61 percent of the true difference, against about 3 percent in the
authors' point-mutation tests; here the dropped scaffolds share sequence with the chromosomes (5 trials, earlier
run)". Mark the row in Table 4 as from the earlier run.

**Verdict: C8 NARROWED.**

## F11 -- "almost all from unplaced scaffolds" is inferred, not counted (C10)

The only source is `prior.md:26-27`: "GigaDB keeps only the 10 chromosomes, so the difference is the unplaced
scaffolds". That is a deduction from assembly structure. No artifact counts how many of the 3,237,326 A-side k-mers lie
on Cannbio-2's unplaced scaffolds, and the 775 B-side k-mers show that the chromosomes are not identical either. The
ledger has no row for the phrase (grep "unplaced" in `ATTACKED-ledger.md`: no hit). The counts 3,238,101, 0.73% and
7,931 reproduce (`truth.tsv`), and KMC independence is the provenance lens's question.

**Correction.** Owed check: locate the recovered A-side k-mers (kp10c already rescans and emits them) by sequence name
and report the fraction on unplaced scaffolds. Until then write "consistent with the release dropping the unplaced
scaffolds".

**Verdict: C10 NARROWED.**

## F12 -- "each party checks its buckets with its own data only", but acceptance is joint (C7)

`kp10c.c:21`: "A bucket is accepted when both sides pass and deg P_rev - deg Q_rev = delta". Each check uses only one
side's data, as the paper says. The decision to accept or split, though, needs both sides' outcomes, so a real
deployment sends a pass/fail bit per bucket per round. That traffic is not counted (about 8 KB for 65,536 buckets,
0.003 B/diff, so negligible). More important, each "round" is a round trip, while the table's "1 round" is a one-way
exchange. Separately, the table's ratio of 1.05 was chosen as the smallest grid value that succeeded 30/30 on the very
trials then reported. That favours the table in bytes, so it does not hurt the sketch's claim, but it means the 9.5%
bound is a post-selection figure for the table (the statistics lens's ground). km-peeler unmodified, r = 3 and the
canonical input all match the results file §2. I did not rebuild it.

**Correction.** "Each party checks its buckets with its own data, and the two exchange one bit per bucket to agree which
to split."

**Verdict: C7 NARROWED (minor).**

## Claims not refuted

- **C5.** The Clopper–Pearson bound is 1 - 0.05^(1/30) = 0.0950, which is correct. "Neither method returned a wrong
  difference": `kp10_seeds.tsv` has missing_total 0 and spurious_total 0, and `iblt_trials.tsv` has
  undetected_failures 0 throughout. The hedge "do not compare reliability" is the right inference. SURVIVED.
- **C9.** Plotting one sketch run against the table's 10-run mean does not mislead here. Every row of
  `iblt_protocol.tsv` has mean_rounds = max_rounds, so all 10 repetitions sent identical bytes and the mean is each
  repetition's value. The sketch's seed spread at the true D is 36.42-36.55 MB (0.36%, `kp10_seeds.tsv:2`). SURVIVED on
  the inference point. The plotted values themselves are the correctness lens's ground. The caption title inherits
  F3's causal defect.
- **C15.** NARROWED. Row 113 sources the 66 to a build-cost sentence that the paper uses as a recovery cost (F8). Row
  57 records the cross-program 2.4x that Figure 3's caption forbids (F7). No row covers "almost all from unplaced
  scaffolds" (F11), and none covers the causal "because" of claim 5 (F3). Row 180 (20.8 vs 38.5) cites a line that does
  say what the row says.

| C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 | C9 | C10 | C11 | C12 | C13 | C14 | C15 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| NARROWED | REFUTED | NARROWED | REFUTED | SURVIVED | NARROWED | NARROWED | NARROWED | SURVIVED | NARROWED | REFUTED | NARROWED | NARROWED | NARROWED | NARROWED |
