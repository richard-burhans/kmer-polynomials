# Round 2 attack, slot L3, lens: statistics

Attacker: L3 (statistics). Frozen text: `ATTACKED-paper.md`, `ATTACKED-ledger.md`. Raw outputs read:
`results/kp10_vs_shibuya_2026-11-11/laila_b69c5c76/` (`run.log`, `tables/*.tsv`), `recon/kp10c.c`,
`recon/kp10_vs_iblt_job.sh`, `recon/aswritten.py`, `results/round1_owed_2026-11-11.md`.

**Appendix.** The brief declared none. I found the method cards myself at `../library/methods/`. The ones
that govern findings here are `pseudo_replication.md` (F1: "effective sample size is the number of independent
decisions, not the number of rows") and `power.md` (F5: state P(the outcome) before reporting it). No
literature claim is made below, so no paper was opened.

**Prior rounds.** `scripts/round history` for "power of two", "phase", "1.05", "second round", "9.5" and
"bucket". Round 1's L2 F7 and L3 F6 (UPHELD) said "equal reliability" was wrong and the sketch's second round was
undisclosed. The draft now discloses both. Nothing below re-files them. No prior round mentions the power-of-two
bucket quantization (searched absence for "power of two"), the IBLT protocol's lattice phase, or the k = 63 sizing
rule.

**Steering.** The claim list asks whether the numbers reproduce and whether the comparison is fair to the table.
It never asks how the sketch's own cost depends on where D falls relative to its power-of-two bucket boundary. That
dependence is the weakest point of claims 4 and the C4 attribution (F1, F2).

All local checks were pure-Python arithmetic on numbers from the tables (no genome was read), run under `nice -n 19`.

## F1 -- the 11.3 against 20.5 is the sketch's best phase on a period-2 sawtooth, and one pair fixes the phase (C1, C13)

**Defect.** The sketch's first-round bucket count is a power of two (`recon/kp10c.c:541`:
`while (dest / 2^bits > 0.8 * m) bits++`). Its bytes per difference therefore depend on the load
L = D / 2^b, which runs from 25.6 to 51.2 as D moves through one octave. The bytes are 536 B per bucket
(67 entries of 8 B; `35,127,296 / 65,536 = 536`).

- On the GigaDB pair, L = 3,238,101 / 65,536 = 49.41. That is 96.5 percent of the 51.2 boundary, so the sketch
  sits at almost the cheapest point of its cycle.
- **The overflow model reproduces the run.** In a Poisson model of bucket loads, P(load > 64 | mean 49.41) × 65,536
  gives **1,253.4** expected overfull buckets. The 30 seeds measured a mean of **1,253.97** (`kp10_seeds.tsv:2`,
  `mean_splits`). The model's bytes per difference, 536/L × (1 + 2 P(overflow)), is **11.263**, against the measured
  **11.26**.
- **Over one octave of D, the same model gives a sawtooth.** The cost runs from **11.2 to 20.9** B/diff, with a
  log-uniform average of **15.2**. The table is sized continuously (n = 1.05 D) and costs about 20.5 at every D
  (20.53 on giga, 20.54 on draft, `iblt_trials.tsv:3,9`).
- **A 3.6 percent overestimate erases the advantage.** Any D_est above 51.2 × 2^16 = 3,355,443 (1.036 D) selects
  2^17 buckets. That is the measured 21.70 B/diff (`kp10_dest.tsv:7`), at parity with the table's 20.5.
  `kp10_dest.tsv` deliberately samples D_est only at exact powers of two times D, so this is not visible in the
  sweep or in Figure 4B.

**Why the statistics card applies.** The 30 seeds replicate the hash, not the comparison. By
`pseudo_replication.md`, the unit that decides the size ratio is the pair and its D. There is n = 1 such unit at
the 3.2M scale, and its phase happens to sit at 0.965 of the boundary. "In 30 of 30 runs" in the abstract
sentence presents seed replication beside a ratio that seeds cannot replicate. Seed spread at fixed D is 0.36
percent (36,415,840 to 36,547,696 B).

**Correction.**
- Abstract: "... with 11.3 bytes per difference over two rounds, against 20.5 ...". Add that this D fell just
  under a power-of-two bucket boundary. Over D the sketch's cost runs from 11 to 21 B/diff (about 15 on average),
  and an overestimate of more than 4 percent costs 21.7.
- Alternatively, report the phase-averaged ratio, about 1.35×, as the headline rather than 1.82×.
- C13's "the sketch's advantage over the table is due to resending only overfull buckets" must be narrowed with
  it; see F2.

**Premise check.** The model is fitted to nothing: it has no free parameter and reproduces both the split count and
the bytes to 0.05 percent. The sawtooth outside D is a model projection. A measured version is owed: run `kp10c`
with D_est at 1.04 D and at 0.6 D, or on a second pair.

## F2 -- "its smaller size comes from the second round" confounds one-round sizing with the power-of-two doubling (C4, C13)

**Defect.** The comparator, 21.7 B/diff, is a one-round sketch at L = 24.7. The expected number of overfull buckets
there is 131,072 × P(Pois(24.7) > 64) ≈ **1.6 × 10⁻⁶**, so it is sized to a failure rate about 10⁵ times lower
than the table's (at most 9.5 percent). It is the next power of two, not the smallest one-round sketch. Buckets are
ranges of S (paper, Methods), so nothing in the design forces the bucket count to a power of two. Only
`v >> (61 - bits)` in `kp10c.c:12` does.

At the reliability the table was held to, the same validated model (F1) gives the smallest one-round sketch on
the GigaDB pair:

| buckets | load | P(any bucket overflows) | B/diff, one round |
|---|---|---|---|
| 96,536 | 33.5 | 0.088 | 15.98 |
| 98,536 | 32.9 | 0.047 | 16.31 |
| 103,536 | 31.3 | 0.0095 | 17.14 |
| 131,072 (paper) | 24.7 | 1.6 × 10⁻⁶ | 21.70 |

The paper's own data already contain a one-round sketch at this load. On the draft pair, L = 7,931 / 256 = 31.0
gives **17.30 B/diff in one round, in 30 of 30 runs, with 0 splits** (`kp10_seeds.tsv:3`). That is below the
table's 20.54 with no second round.

So of the 1.82× gap, about 20.5 / 16.3 = 1.26× is available in one round at matched reliability. About
16.3 / 11.3 = 1.44× comes from the second round.

**Correction.**
- C4/Results: "Sized for a single round with a power-of-two bucket count, the sketch needed 21.7, but at the draft
  pair's load (31 per bucket) it used 17.3 in one round, so part of its advantage over the table comes from the
  second round and part from the bucket count."
- Discussion: "Because an overfull bucket is resent on its own, the sketch spent 11.3 ..." should be narrowed the
  same way.

**Premise check.** I did not run a non-power-of-two sketch, because `kp10c` cannot do one. The 16 to 17 B/diff
figures are model projections, anchored by the measured 17.30 at L = 31 and by the model's exact reproduction of
the 49.4-load run.

## F3 -- the eightfold row samples the table's cheapest phase; at other near-eightfold estimates the table costs up to about 72 (C2, C9)

**Defect.** The table's protocol starts at 1.05 D_est and doubles. At D_est = D/8 exactly, its fourth table is
n = 3,400,008 = 1.05 D (`run.log` `iblt_protocol` lines, `n_per_round` "425001, 850002, 1700004, 3400008" in all 10
reps). That is precisely the size already known to succeed in 30 of 30, so the protocol lands on its best case.

Using the measured success rates (0/30 at 0.95 D, 18/30 at 1.0 D, 30/30 at 1.05 D; `iblt_trials.tsv:7-9`) and
19.55 B per unit of n:

| D_est / D | table, resend at double (expected B/diff) | sketch |
|---|---|---|
| 0.1131 | 72.0 | 20.76 (bits0 = 13) |
| 0.120 | 50.2 | 20.76 |
| 0.125 (paper) | 38.5 | 20.76 |
| 0.150 | 46.2 | 20.76 |
| 0.200 | 61.6 | 20.76 |

- **The sketch's cost depends only on bits0.** bits0 = 13 for every D_est in (0.0648 D, 0.1295 D], and
  `kp10_dest.tsv` sampled exactly one D_est per bits0 value.
- **The table's cost is continuous in D_est** and ranges about 38.5 to 72 over that span. The paper reports its
  minimum as "the" cost of an eightfold underestimate.

**Direction.** This error favours the table. The abstract's "20.8, against 38.5" understates the sketch's
advantage under misestimation; the ratio is 1.85× here and up to about 3.5× nearby. I file it because the number is
phase-selected, not because the conclusion reverses.

**Correction.** Report the eightfold comparison as "20.8 against 38.5 to 72, depending on where the estimate falls
relative to the table's doublings". Or sample D_est off the power-of-two lattice (for example 0.12 D and 0.9 D) in
the sweep.

**On Figure 4B's "mean".** All 10 repetitions at every point had identical round counts and identical bytes
(`iblt_protocol` lines; `mean_rounds` = `max_rounds` in `iblt_protocol.tsv:2-9`). The plotted "mean" is therefore
each repetition's value, and the repetitions carry no information about bytes, only about success. Plotting one
fixed-seed sketch run against it is not misleading, because sketch bytes vary 0.36 percent over seeds at the true
D. The caption should say the table's bytes were identical in all 10. That leg of C9 survives.

## F4 -- the table's "independent random seeds" are the fixed integers 1001 to 1030, reused at every ratio and in both pairs (C7; C5 survives)

**Defect.** `recon/kp10_vs_iblt_job.sh:176` passes seed `1000 + t`, t = 1..30, to `ibltseq build -s`, and the same
30 seeds run at every ratio and in both pairs. `run.log` shows each of seeds 1001 to 1030 exactly 10 times: 5 giga
ratios plus 5 draft ratios. The protocol seeds are `100000k + 1000r + 1000f`. These seeds are not random. Their
independence as trials rests entirely on km-peeler's seeded hash family, which this round did not check. Checking
that `ibltseq`'s hash with seeds s and s+1 gives unrelated cell assignments is owed.

The sketch's 30 seeds are genuinely OS-random: distinct x and λ per run, `"seed_x": "os"` in `run.log` `kp10_seed`
lines.

The table's 1.05 ratio was also chosen in the 2026-11-10 run on parts 0 to 3 of the same giga pair, as the smallest
grid ratio with 30 of 30. The caption itself says "sized to succeed in 30 of 30". Its present 30 of 30 is
therefore a correlated confirmation on a superset of the selection data, not a test fixed before any outcome was
seen. The sketch's 0.8m rule had no such tuning step.

**C5.** The arithmetic is right: 1 − 0.05^(1/30) = **0.09503**, one-sided 95 percent Clopper-Pearson. "Neither
method returned a wrong difference" matches `kp10_seeds.tsv` (0 missing, 0 spurious) and `iblt_trials.tsv`
(`undetected_failures` 0). The conclusion "do not compare reliability" is the conservative one. C5 SURVIVES. The
seed caveat above weakens only the table's half of the bound, and in the direction of less information.

**Correction (C7).** Write "each method ran 30 times on one host, the sketch with independent OS-random seeds and
km-peeler with seeds 1001 to 1030, the same at every table size". Also say that the 1.05 sizing was chosen from an
earlier run on part of the same pair.

## F5 -- the k = 63 run used a different, doubled sizing rule from the one Methods states, so its "one round" had no chance to fail (C3)

**Defect.** Methods: "The first round uses the smallest power of two of buckets that holds at most 0.8m estimated
differences each." The k = 63 run instead used `bits_for` in `recon/aswritten.py:136-142`, which allows m/2 = 32 per
bucket ("expected difference per bucket <= m/2"). `results/round1_owed_2026-11-11.md:119` confirms that it followed
the bits_for rule.

- At D = 21,382 that gives 2^10 buckets at load 20.9: 548,864 B, 25.7 B/diff.
- Under the paper's stated 0.8m rule it would be 2^9 at load 41.8: **274,432 B, 12.8 B/diff**.

**Power (`power.md`).** At load 20.9 with m = 64, the expected number of overfull buckets is
1,024 × P(Pois(20.9) > 64) ≈ **10⁻¹¹**. "One round" was therefore certain before the run, and so was exact decoding
given correct hashing. Under the stated rule the expectation would be 0.27 overfull buckets, which is an
informative test.

The same holds at k = 31 as written: 256 buckets at load 31.0 give an expectation of 1.7 × 10⁻⁵. There the two
rules agree, so the paper's 137 KB is unaffected.

The informative part of C3 is deterministic and survives. The S-set difference equals the canonical difference
(7,931 and 21,382, matched element by element), against 398,217,248 for the forward-only hash. Only "549 KB and
one round" misdescribes what was tested.

**Correction.** Either state "buckets sized from the true difference at m/2 = 32 differences each, twice the
count the 0.8m rule gives", or rerun at 2^9 and report that size. The 549 KB must not be read against the 0.8m
numbers in Table 4.

## Verdicts

| claim | verdict | finding |
|---|---|---|
| C1 | NARROWED | F1: the ratio is set by D's phase against a power-of-two boundary (n = 1 pair) |
| C2 | NARROWED | F3: 38.5 is the table's best phase; the conclusion holds and is understated |
| C3 | NARROWED | F5: the k = 63 run used the m/2 rule, not the stated 0.8m; its "one round" had no power |
| C4 | NARROWED | F2: the one-round comparator is the next power of two; at matched reliability about 16 to 17 |
| C5 | SURVIVED | F4: the bound is 0.09503, and the no-wrong-difference claim matches the tables |
| C6 | SURVIVED | 17.30 and 20.54 reproduce (`kp10_seeds.tsv:3`, `iblt_trials.tsv:3`), one round each, 30/30 each |
| C7 | NARROWED | F4: the table's seeds are fixed 1001 to 1030, reused; 1.05 was tuned on part of the same pair |
| C8 | UNCHALLENGED | |
| C9 | SURVIVED | F3: the plotted values match `kp10_dest.tsv` and `iblt_protocol.tsv`; the table's 10 reps are byte-identical |
| C10 | UNCHALLENGED | |
| C11 | UNCHALLENGED | |
| C12 | UNCHALLENGED | |
| C13 | NARROWED | F2: part of the advantage needs no second round (17.3 in one round at load 31) |
| C14 | UNCHALLENGED | |
| C15 | SURVIVED | Ledger rows 102, 152, 153, 155 and 175 to 181 cite `results/kp10_vs_shibuya_2026-11-11.md` lines 104 to 107, 111 to 118 and 141, which say what the rows say |

C1 NARROWED · C2 NARROWED · C3 NARROWED · C4 NARROWED · C5 SURVIVED · C6 SURVIVED · C7 NARROWED · C8 UNCHALLENGED · C9 SURVIVED · C10 UNCHALLENGED · C11 UNCHALLENGED · C12 UNCHALLENGED · C13 NARROWED · C14 UNCHALLENGED · C15 SURVIVED
