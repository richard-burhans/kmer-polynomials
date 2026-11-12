# Round 2 adjudication of L5 (inference)

Attack read: `reports/rounds/round-2/attack-L5-inference.md`. Its sha256 begins `509e54e1891e`, matching the seal
(`round.json`: `509e54e1891e4f60...`). Read in the harness order: the attack first, then the claims and their evidence.

**What I re-ran.** One `nice -n 19 ionice -c3` Python process, on committed tables and Poisson arithmetic only. No
genome data was touched. Checks:

- the raw tables `laila_b69c5c76/tables/{kp10_dest,kp10_seeds,iblt_protocol,truth}.tsv`;
- `run.log:44,48,53` (seeds);
- the `results/round1_owed_2026-11-11` JSONs for k63 and fwd;
- `recon/kp10c.c:1-31,540-545`;
- `results/kp10_vs_shibuya_2026-11-11.md` §2, §4, §5 and §6;
- `results/kp10_vs_shibuya_2026-11-10_prior.md:1-30,95-116`;
- `ATTACKED-paper.md:14,23,63,67,123,135-156`;
- `ATTACKED-ledger.md` rows 50-60, 110-115, 178-200.

Every number the attacker computed reproduces:

| quantity | reproduced value |
|---|---|
| Poisson overfull buckets at 2^16 | 1,253.36, against 1,253.97 observed |
| one round at 100,000 buckets | 16.55 B/diff, P(any overflow) 0.030 |
| sawtooth over one octave | min 11.20, max 20.94, log-uniform mean 15.16; ratio ≥ 1.82 in 5.75% of it |
| D/8 bytes | (2^13+2^14+2^15+2^16)·536 + 2·1,258·536 = 67,212,256, exactly |
| \|A\|−\|B\| | 3,236,551 (giga) and 7,931 (draft) |
| k63 run | bits 10 under `bits_rule: "auto: bits_for(exact value D, m)"` |
| fwd run | D_recovered 0, unresolved 256 |
| seeds | `seed_x`/`seed_lam` are 0x5eed/0x1a3d on line 44 and "os"/"os" on line 48 |

**Earlier rounds.** `scripts/round history "own differing"` returns round 1's L5 F5. That ruling UPHELD the
oracle-decoder point and OVERTURNED an attacker sentence that described a two-party decoder that had not been built. kp10c
has since built that decoder (`kp10c.c:15-23`), so this round's F2 tail and F12 are not re-filings. The attacker
disclosed the round-1 ruling. `scripts/round history "resent whole"` finds no prior mention.

**Method cards.** No card in `../library/methods/` governs sketch sizing or comparing protocols. `causal_designs.md`
governs the "because" in F3 (a causal reading needs a contrast that nets out what the arms share), and I apply it
there. `power.md` already settled C5.

---

## F1: "its smaller size comes from the second round" is power-of-two rounding (C4, C6). NARROWED

**What is right.** The numbers in the diagnosis are correct:

- `kp10c.c:541` allows only power-of-two bucket counts.
- The "sized for one round" row is the dest = 2D run (`kp10_dest.tsv:7`, load 24.7 = 0.39m).
- The Poisson model is calibrated on this run to 0.05%.
- Within that model, about 100,000 buckets would finish in one round (P(any overflow) ≈ 3%) at about 16.5 B/diff.

So the general inference that the advantage "comes from the second round" is confounded with granularity.

**Why not REFUTED.**

- As a statement about the sketch as implemented (power-of-two buckets), the sentence is a measured fact. The smallest
  bucket count that finished in one round, 2^17, cost 21.7. The next one down, 2^16, failed in round 1 every time
  (`runs_finished_in_1_round` 0 of 30).
- The 16.5 counterfactual and the "roughly half" split are model predictions. No run measures them.

**C6 leg.** 17.3 at a load of 0.48m is real and reproduces. The attacker is right that it marks where 7,931 falls
against a power of two. NARROWED, not refuted.

**Repairs.**

- **Endorsed, free:** "Sized for a single round **with power-of-two buckets**, the sketch needed 21.7".
- **Endorsed, free:** relabel the Table 4 row "2^17 buckets (next power of two)".
- **Endorsed, free:** note that the draft pair's one-round sketch sits at a load of 0.48m.
- **Not endorsed for the text:** "about 16.5" and "roughly half the gap". These stay out until a run measures them. A
  kp10c run with a non-power-of-two first round is **owed and costs a laila run**. It needs a code change to the
  bucket key (top bits → multiply-high on the value range).

## F2: 11.3 vs 20.5 is the best point of a sawtooth (C1). NARROWED

**What is right.**

- Under the implemented rule the sketch's cost depends on where D falls in the octave, and the table's does not. My
  recomputation gives the same figures: min 11.20, max 20.94, mean 15.16, 5.75% of the octave at ≥ 1.82×.
- The abstract never says that the ratio varies with D.
- The abstract also omits "each party recovers only its own k-mers". The Results discloses it (line 146), and it bears
  on an abstract that sets bytes against a table whose single listing gives both sides.

**Limits.** 11.3 is a measured value at a stated D, and the abstract does not claim it is typical. The defect is an
omission of scope, not a false number.

**Repairs.**

- **Endorsed, free:** give both pairs in the abstract (11.3 on the GigaDB pair, 17.3 on the draft pair) as the range
  seen.
- **Endorsed, free:** add "each party recovers its own differing k-mers".
- **Not endorsed:** "averages about 15 by a model fitted to this run" in the abstract. It puts a model-derived average
  of an unmeasured octave into an abstract-level claim.

## F3: "Because an overflowing bucket is resent alone" does not explain the D/8 result (C2). NARROWED

**What is right.** On the GigaDB pair the first three rounds split every bucket. The sketch and the table degrade by
the same factor from their own true-D cost: 1.84× and 1.87×. "Resent alone" therefore does **not** make the sketch less
sensitive to an underestimate than the table. A reader would take that from "an eightfold underestimate raised this
[only] to 20.8".

**Where the attacker overreaches.** The conclusion that the "because" is false does not follow. The correct contrast
for "because X" holds the sketch fixed and removes X (`causal_designs.md`), and the sketch's own measured sizes give
that contrast:

| sketch that resends whole at double size | B/diff | table |
|---|---|---|
| at the true D (2^16 fails, then 2^17) | 32.5 | 20.5 |
| at D/8 (2^13 … 2^17) | 42.0 | 38.5 |
| draft pair at D/8 | 32.4 | 38.55 |

Whole-resend sketch costs are my arithmetic from the measured bucket sizes and fail patterns in `kp10_dest.tsv`.

- Without per-bucket resend, the sketch loses to the table at both D and D/8 on the GigaDB pair.
- With per-bucket resend it costs 11.3 and 20.8.
- So, within the implemented design, resending alone does cause the level gap at D/8.
- What it does not cause is a smaller penalty for underestimating.
- The attacker's own observation that the D/8 gap "is the true-D gap carried through" concedes this: the true-D gap is
  itself the final-round partial resend (see F1).
- On the draft pair, partial resend reduced even the degradation: 1.23× for the sketch against 1.88× for the table.

**Repairs.**

- **Endorsed, free:** append to the abstract sentence, or substitute for its framing, "both resend at double size until
  they fit, so each costs about 1.85 times its own true-D size on this pair; the gap is the one at the true D".
- **OVERTURNED:** dropping "Because an overflowing bucket is resent alone" from the abstract.
- **OVERTURNED:** deleting the Figure 4 caption title "Resending only overfull buckets keeps the sketch small when the
  difference is underestimated". It is true against the sketch's own whole-resend counterfactual (20.8 vs 42.0), and
  more so on the draft pair. Deleting it would install the false converse.

## F4: an eightfold underestimate cannot happen on these pairs; the table's protocol ignores free information (C2, C13). NARROWED

**What is right.** |A| − |B| is 3,236,551 = 0.99952 D on the GigaDB pair and exactly 7,931 on the draft pair. On these
one-sided pairs the D/8 underestimate is imposed rather than realistic, and the paper should say so.

**Wrong premise.** "Only the table's naivety is priced into the abstract" is false:

- Both protocols are run equally naively, and both costs are priced: 20.8 against 11.3 for the sketch, 38.5 against
  20.5 for the table.
- The sketch carries per-bucket sizes (`kp10c.c:12-13`), so it could use the same information.
- A sketch that restarts at 2^16 after a failed round at 2^13 would cost about 10.85·1.125 + 0.42 ≈ 12.6. That keeps the
  ratio near 0.55.

**Repairs.**

- **OVERTURNED:** "a table resized from its own failed counters would cost about 23". It gives the information to one
  arm only. It also sets an untested success rate (1.0495 D sits between 18/30 and 30/30) against a measured one.
- **OVERTURNED:** "drop claim 5". The measurement is sound. Its scope is what needs saying.
- **Endorsed, free:** "On these pairs the two set sizes alone estimate the difference to within 0.05 percent, so the
  underestimate is imposed to test the protocols."
- **Endorsed, free:** state the doubling protocol as the one km-peeler's interface supports (results §4 already says
  its partial listing cannot be reused).

## F5: "an IBLT that fails must be resent whole" is a protocol choice (C13). UPHELD

The Introduction's universal "must" (line 23) generalizes from km-peeler to all IBLTs:

- Partitioned or strata-style IBFs (eppstein_2011) can resend only the failing parts.
- The results file itself concedes a smarter protocol was "not tested" (§4).

The attacker's own caveat is fair: Eppstein's small-delta overhead means a partitioned IBLT need not win. The finding
is about the universal wording and about attributing the gain to the power sums, not about the outcome.

**Repairs.**

- **Endorsed, free:** "whereas km-peeler's table, when it fails, is resent whole".
- **Endorsed, free:** in the Discussion, "a table partitioned the same way was not tested".
- The rateless-IBLT citation (Yang, Gilad and Alizadeh, SIGCOMM 2024) was cited from memory and is not held in the
  library. **I did not verify it.** It is owed to the prior-art lens before anyone cites it.

## F6: "needs no canonicalization", "values that must be added", and the k = 63 sizing (C3). NARROWED

1. **"Values that must be added".** Linearity belongs to the power-sum sketch over any field-valued IDs. A canonical
   code with a 2-universal hash into F_q is just as addable. The Discussion lists "values that must be added" as a
   place where "that recipe is costly or unavailable", and at k ≤ 32 nothing measured supports that.
   - **Endorsed, free:** qualify it as "for k > 32" or drop it.
   - **Not required:** rewording the abstract's "needs no canonicalization". As written it is true of what was run:
     the as-written runs, with canonicalization only for scoring against the truth (`kp10c.c:24-25`).
2. **"Saw 398,217,248 differences".** The fwd JSON shows the sketch recovered 0 and left 256 buckets unresolved. The
   398,217,248 is `D_true_values`, computed outside it. UPHELD.
   - **Endorsed, free:** "a forward-only hash faced 398,217,248 differences and overflowed every bucket" (line 150), and
     the same in the Introduction (line 23, "saw 398 million").
3. **k = 63 sizing.** `bits 10` (load 20.9 = 0.33m) is not the 0.8m rule, which would give 2^9 and 274,432 B. The paper
   says "buckets sized from the true difference" but not by which rule. UPHELD.
   - **Endorsed, free:** state the k = 63 run's rule and that it chose 2^10. 549 KB is 25.7 B/diff.
   - The point that the long-k niche has no comparator is the same as F7's tail.

## F7: 2.4× crosses two programs that Figure 3's caption says not to compare; the long-k niche is unmeasured (C14). NARROWED

1.12 = 11.43/10.24 (same program) and 2.4 = 24.67/10.24 (across programs) both reproduce from the ledger rows.

- The two programs' ntHash baselines agree to 0.5%, so the cross-program ratio is probably sound.
- The paper does not say why the caption's prohibition, which rests on a sampler-specific 28.7 vs 22.5 gap, does not
  apply here. That is the defect.
- "With a proof" survives.

**Repairs.**

- **Endorsed, free:** the attacker's parenthetical on the two programs.
- **Endorsed, free:** the long-k niche worded as expected, not measured.

## F8: "66 multiplications per element" is a build cost cited as a recovery cost (C12). UPHELD

**The defect.**

- The source (`prior.md:8`) gives 66 as the cost of the power sums, which is construction.
- 3 cell updates is IBLT construction too.
- The paper's preceding sentence says "Recovery", and the measured recovery cost (110.3 s vs 3.5 s, results §5) is
  omitted.

**The repair, amended.** The 110.3 s is "both rounds, **both sides'** root flagging" in one process, single-threaded.
km-peeler's 3.5 s is one listing. The attacker's sentence "decode took 110 s against 3.5 s" is therefore not per party.

- **Endorsed, free:** "Building the sketch takes m + c = 66 field multiplications per element, against 3 cell updates.
  Decoding, which evaluates the recovered polynomial at the elements of every differing bucket, took 110 CPU-s for
  both parties together on one thread, against 3.5 s for the table's listing."

## F9: "the GigaDB sketch ran before x and λ were separated" describes a superseded run (C11). UPHELD

- `run.log:44` carries separate fixed seeds (0x5eed/0x1a3d), and `run.log:48-105` carry `os`/`os`.
- The sentence was added in round 1 (`c3f3a11`), before job b69c5c76.
- The first sentence of C11 matches `kp10c.c:8-12` and survives.

**Whole-set check (rule 3).** Deleting only the clause, as prescribed, leaves the sentence before it at line 67: "Our
experiments used fixed seeds for reproducibility, so they test the construction rather than the bound". That is false
for the 30-run comparisons, which drew OS seeds. It also contradicts C7's "independent random seeds".

- **Endorsed, free, with that amendment:** delete the clause and write "The sweeps used fixed seeds; the 30-run
  comparisons drew both seeds from the operating system."
- The attacker missed the companion sentence, and nobody should apply this fix without it.

## F10: the extended-syncmer row has the wrong unit, unlike regimes and a superseded run (C8). NARROWED

**What holds.**

- The unit is wrong: 1,986,807 are spurious **k-mers** (`prior.md:105`; 61.4% of D), not syncmers.
- The 3% comes from point-mutation tests, and the source's own explanation (scaffolds sharing sequence with the
  chromosomes) is missing from the paper.
- The row comes from the 5-trial 2026-11-10 run, with n sized from the true extended-syncmer difference.

**What does not.** "2.7 times smaller" is not hidden as exact-vs-inexact: the same sentence reports the spurious and
missed counts.

- **Endorsed, free:** the attacker's replacement sentence.
- **Endorsed, free:** marking the Table 4 row as coming from the earlier 5-trial run.

## F11: "almost all from unplaced scaffolds" is inferred, not counted (C10). NARROWED

The phrase rests on a deduction (`prior.md:26-27`), and no artifact locates the k-mers.

The deduction is not airtight. The 775 B-side k-mers bound substitutions on the chromosomes, but not deletions or
trimming on GigaDB's side, which add A-side k-mers without adding B-side ones. All counts reproduce (`truth.tsv`).

- **Endorsed, free:** "consistent with the release dropping the unplaced scaffolds".
- The check that locates the k-mers by sequence name is **owed and costs a laila run**: a rescan that emits sequence
  names.

## F12: "each party checks its buckets with its own data only", but acceptance is joint (C7). NARROWED

**What holds.**

- `kp10c.c:21` makes acceptance joint.
- The 1.05 ratio was chosen from the trials it reports. That post-selection favours the table in bytes.

**What the attacker understates (whole-set check).**

- Each side's check needs the other polynomial's coefficients (results §2: "power sums of the other polynomial's roots
  computed from that polynomial's coefficients"). The party that did not receive the sketch must therefore get either
  the other sketch or the decoded polynomial for its side, not just one bit per bucket.
- On these pairs this is cheap only if the receiver is the side with the large difference. Cannbio-2 would receive and
  return Q, about 775 coefficients plus per-bucket headers.
- In the other direction it would cost up to 8 B per difference.
- The bytes stay "one side's message", as the results file says (§6), but that rests on the right choice of direction.

**Repair.**

- **Endorsed, free, amended:** "Each party checks its buckets with its own data. The receiver returns, per bucket, a
  pass/fail bit and the decoded polynomial for the other side, which are not counted (about 775 coefficients on the
  GigaDB pair when Cannbio-2 receives)."
- The exact uncounted bytes are free arithmetic from existing outputs. This is for the correctness lens.

## Claims the attacker left standing

- **C5: SURVIVED.** I concur. 1 − 0.05^(1/30) = 0.0950. `missing_total` and `spurious_total` are both 0.
- **C9: SURVIVED** on the inference point. I concur: `mean_rounds = max_rounds` in every `iblt_protocol` row. Per F3,
  the caption title is **not** defective.
- **C15: NARROWED.** I concur on the missing and mis-sourced rows (66, 2.4×, "unplaced"). I do not concur that the
  "because" needs a row as a defect (F3).

**Correct minority (rule 8).**

- F9's companion "fixed seeds" sentence and F12's direction-dependent uncounted traffic are additions of this
  adjudication. No attacker raised them here.
- F3's overreach is the attacker's alone.

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | true of power-of-two buckets as implemented and measured; the 16.5 counterfactual is a model, so a non-power-of-two run is owed; the "with power-of-two buckets" wording and the Table 4 label are free |
| F2 | NARROWED | numbers reproduce; an omission of scope, not a false number; quoting both pairs and "each party recovers its own" is free; the model mean in the abstract is not endorsed |
| F3 | NARROWED | the underestimate penalty is the same for both (1.84× vs 1.87×), but resending alone does cause the level gap (a whole-resend sketch costs 42.0 at D/8); add the equal-degradation clause for free; dropping "because" and the caption title is OVERTURNED |
| F4 | NARROWED | the underestimate is imposed on these pairs (disclosure free); the one-sided "table resized from counters ≈ 23" and "drop claim 5" are OVERTURNED |
| F5 | UPHELD | "must be resent whole" holds only for km-peeler; both rewordings free; the rateless citation is unverified and owed to the prior-art lens |
| F6 | NARROWED | "saw 398 M" and the unstated k = 63 rule are defects (fixes free); "values that must be added" needs "k > 32" (free); "needs no canonicalization" stands |
| F7 | NARROWED | cross-program 2.4× needs its justification stated; the long-k niche is expected, not measured; both fixes free |
| F8 | UPHELD | 66 is the build cost, not recovery; fix free, amended to say the 110 s covers both parties |
| F9 | UPHELD | the GigaDB runs had separate (and OS) seeds; the deletion is free but must also fix the companion "fixed seeds" sentence |
| F10 | NARROWED | spurious k-mers, not syncmers; earlier 5-trial run; unlike regimes; fixes free |
| F11 | NARROWED | a deduction, not a count; interim wording free; locating the k-mers costs a laila run |
| F12 | NARROWED | joint acceptance is real; the uncounted return is the other side's polynomial, not just a bit; direction-dependent; amended fix free |
| C5 | UPHELD | (attacker's SURVIVED stands) bound and zero-error counts reproduce |
| C9 | UPHELD | (attacker's SURVIVED stands) identical bytes across repetitions; the caption title is not defective |
| C15 | NARROWED | mis-sourced and missing rows real; no row is owed for the "because" |
