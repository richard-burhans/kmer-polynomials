# Round 2 adjudication, lens L7: reproduction

Adjudicator: skeptic of `attack-L7-reproduction.md`. **Digest checked:** sha256 `2692cd5685d85b5d...` matches the seal
(`round.json:177`). I read the attack first and the claims second, as the harness fixes. I did not read other round-2
attacks or adjudications.

**What I re-ran (local, `nice -n 19 ionice -c3`, raw outputs only).** I parsed the 496-line `run.log` of
`laila_b69c5c76`: 300 `iblt_trial`, 80 `iblt_protocol`, 60 `kp10_seed`, 12 `kp10_dest` and 2 `kp10_timed` records. I
read `kp10c.c` (header, `decode`, `flag_roots`), `kp10_vs_iblt_job.sh` (seeds, `iblt_once`), `aswritten.py`
(`bits_for`, `s_of_strings`), `recon61.py:118-130`, the k = 63 run JSON, `prior.md:28-36,97-113`,
`paper_bench_2026-11-12.md:28-55`, and round 1's `adjudication-L5.md` F5. Nothing heavy was run.

**Method cards.** `power.md` governs F9 only, and the attacker applied it correctly. No card governs protocol
accounting (F1, F8), and I say so.

---

## F1: "every difference with 11.3 bytes" prices one of two messages. NARROWED

**What reproduces.**
- **The arithmetic.** 2 × 36,481,232 / 3,238,101 = 22.53, 2 × 67,212,256 / 3,238,101 = 41.51, and
  2 × 137,216 / 7,931 = 34.60.
- **The decoder.** `decode()` starts from `SA - SB` (`kp10c.c:259`). Acceptance runs both sides' checks
  (`:300-304`).
- **What each side recovers.** `run.log:28` records both sides' k-mers.
- **The table's bytes.** km-peeler's bytes are `a.iblt` alone (`iblt_once`). One `ibltseq list` gave both sides
  (`run.log:30`).

So, **for the GigaDB pair as run**, "recovered every difference" consumed both sketches, and the paper prices one. The
diagnosis is new since round 1. Round 1 ruled "the byte count does not depend on who decodes" when only a one-process
decoder existed. `kp10c` now makes the two-party structure explicit. The attacker disclosed the prior ruling correctly.
This is not a re-filing.

**Where the scope fails.**

1. **The draft row (C6) is wrong.** All 7,931 differences are on B's side (`run.log:20`: `A_minus_B 0`). In `decode()`
   with eA = 0, P_rev = 1:
   - `flag_roots` flags nothing on A, and nfa = 0 = eA holds whatever A holds.
   - `pP = 0`, so A's check (`-pQ_j = S_j`) and B's check (`-psum(fb) = S_j`) coincide once nfb = eB.

   B therefore reaches exactly kp10c's decision from A's sketch alone, and receives every difference as k-mers. The
   draft comparison is one message against one message: 17.3 against 20.5 is like for like. The prescribed
   "34.6 vs 20.5" would install a stronger false claim, a reversal that does not exist. I verified this by reading the
   code. It needs no run.
2. **"Reverses on every row" applies the protocol as run, not the sketch.** On GigaDB, a single B→A message gives A
   3,237,326 of the 3,238,101 differences (99.98%) as k-mers. Only B's 775 remain as roots of an unfactored Q, and a
   short reply of 775 k-mers would close them. That protocol was not built. So 22.5 is the cost of the symmetric
   protocol that ran, not a measured lower bound for the sketch.

   The attacker applies the opposite standard in F8: there the table's 38.5 is "an upper bound, not its cost" because
   a smarter protocol is conceivable. The same reasoning makes 22.5 and 41.5 upper bounds too. At D/8, 41.5 against
   38.5 compares two upper bounds.

**Rulings on the fixes.**

| fix | ruling | cost |
|---|---|---|
| Abstract and Results, GigaDB: "11.3 bytes per difference sent in each direction (22.5 in total)" | ENDORSED as a true description of the run | free |
| Add that the table's one message gives both sides to one party | ENDORSED | free |
| C2: label 20.8 as per direction (41.5 both ways) | ENDORSED; the paper must not say the reversal is intrinsic | free |
| "restrict the comparison explicitly to a one-way message" | OVERTURNED for GigaDB: one way, the sketch does not deliver every difference as k-mers, so the sentence would compare unlike outputs | none |
| Draft 34.6 vs 20.5 | OVERTURNED (point 1) | none |
| Optional for C6: "all 7,931 differences lie on the draft's side, so one sketch suffices" | ENDORSED (code reading, not a run) | free |
| Owed check, one-sided acceptance on a tiny pair | ENDORSED as a diagnostic. It needs a kp10c flag and a dry run, and it cannot by itself support a one-way "every difference" figure for GigaDB, because the 775 still need a reply or root-finding | small local run plus a code change |

**What survives.** For GigaDB as run, 11.3 is per direction and the two-way total is 22.5. The draft's 17.3 is a
genuine one-message figure. No reversal is established as a property of the sketch.

## F2: the "x and λ were separated" clause is stale. UPHELD

**Verified.**
- `8c34af9` (01:14) is an ancestor of `e76fad1` (13:47), the code the results file names (`:57`).
- `run.log` records the two seeds separately: `seed_x 0x5eed`, `seed_lam 0x1a3d` on all 14 fixed-seed records, and
  `os`/`os` on all 60 seed runs.
- kp10c takes them as separate flags.

**Whole set.** The same sentence also opens "Our experiments used fixed seeds". That is false for the 30-run results
behind "30 of 30" and "11.3", which used OS seeds. The attacker's replacement caveat fixes both halves.

**Fix.** ENDORSED: delete the clause and use the attacker's caveat. The fix is free.

## F3: the k = 63 provenance sentence. UPHELD

**Verified.**
- **The hashes come from C.** `aswritten.py:11-12` streams `hash/kp7e_stream.c`, and `s_of_strings` uses
  `kp7e_direct`.
- **The 20,000-k-mer check is out of scope.** It is the numpy ports of the k = 31 packed pipeline (`prior.md:31-32`).
- **The sizing rule differs.** The run JSON gives `bits 10` and `bits_for(exact value D, m)`, which uses m/2 = 32.
  Under the stated 0.8m rule, 21,382 / 512 = 41.8 ≤ 51.2, which would give 2^9.
- **The bytes check out.** 548,864 = 1,024 × 67 × 8.
- **The decoder.** `recon61.py:122-128` checks `psum(fa) - psum(fb)` in one process.

**Scope.** The m/2 point is a disclosure, not a refutation. The sentence says "sized from the true difference" and does
not invoke 0.8m. Whether 2^9 finishes in one round is unmeasured either way: a Poisson load of 41.8 against m = 64
predicts about 0.4 overfull buckets of 512, so one or two rounds.

**Fix.** ENDORSED with one amendment: in place of the dropped 20,000-k-mer clause, say the C roller is the one the tests
check bit for bit against the reference roller. The fix is free. The 2^9 rerun or a k = 63 kp10c is correctly flagged as
an owed laila run. It is not required for the corrected sentence.

## F4: "resent about 1,250" counts failures, not buckets sent. UPHELD

**Verified.**
- `run.log:26-27` records 1,263 failed buckets and 2,526 buckets sent.
- The 30-seed mean of splits is 1,253.97.
- 2,526 × 536 = 1,353,936 = 36,481,232 − 35,127,296.

**The power-of-two aside.** It is correct, but it does not refute "comes from the second round". Within this design
the only one-round alternative is 2^17, so the attribution holds as stated, and C4 needs no change.

**Fix.** ENDORSED for the caption and ledger row 178. The fix is free.

## F5: the extended-syncmer row comes from the superseded 5-trial run. NARROWED

**Upheld.**
- `prior.md:99-104` gives 5 trials, and the Methods sentence "each method ran 30 times with independent random seeds"
  reads as covering this row.
- `prior.md:104` does record "507,219 of 507,221 every trial", and the paper omits it.
- The 61% against 3% comparison holds as like for like, as the attacker found.

**Narrowed.** "Bears on 'Neither method returned a wrong difference'" is speculation. That sentence is about the sketch
and the k-mer table. The syncmer row already reports its own spurious and missing k-mers in the table. The exit status
is unrecorded, so a harness error is neither shown nor excluded, and round 1's gap 5 stays open.

**Fix.** The caption disclosure is ENDORSED and free. Use the source's literal wording: "listed 507,219 of 507,221
differing extended syncmers in every trial". The rerun in the b69c5c76 harness is owed and costs a laila run.

## F6: km-peeler's seeds were 1001-1030, not independent random draws. UPHELD

**Verified.** All 300 `iblt_trial` records used seeds 1001-1030, each seed 10 times, across all ratios and both pairs.
KP-10's 60 seed runs used `os`.

**Fix.** ENDORSED. The fix is free. The fixed seeds do not undermine the 9.5% bound in F9: 30 distinct seeds are 30
trials. Whether km-peeler's hash functions are independent across consecutive seeds is unknown, and the attacker
correctly left that open.

## F7: "about 66 multiplications per element" omits the decode. UPHELD

**Verified.**
- `flag_roots` (`kp10c.c:228-237`) runs Horner of degree eA on every own element of every bucket with a difference.
- All 65,536 buckets had a difference in round 1, and A's mean degree is 49.4.
- `timed.tsv`: decode took 110.3 CPU-s, against 49.5 s of power sums per side and 3.5 s for the IBLT's diff and list.
- The draft decode (67.2 s with B's degree about 31) fits the same model.

**Fix.** ENDORSED. The fix is free.

## F8: the table's failed listings were correct and were thrown away. UPHELD on the diagnosis; fix NARROWED

**Verified from `run.log`.**

| pair, ratio | failed trials | mean listed correctly | spurious |
|---|---|---|---|
| giga, 0.95 | 30 | 1,353,288 (41.8%) | 0 |
| giga, 1.0 | 12 | 2,021,206 (62.4%) | 0 |
| draft, 1.0 | 22 | 4,701 (59.3%) | 0 |

**Why "must" fails.** A second IBLT built from A's full set at a smaller size, with B cancelling its own set and the
already-listed items, is a valid continuation. "Must be resent whole" is therefore false as a property of the table,
and the 38.5 is an upper bound. The causal "because" in C2 does not separate the sketch's design from the harness's
protocol.

**Fix narrowed.**
- "is, in km-peeler, resent whole" misattributes. km-peeler implements no resend protocol; our harness does
  (`kp10_vs_iblt_job.sh:18`).
- Write instead: "in the protocol we ran, a table whose listing fails is resent whole at double size; its partial
  listings, which were correct, were not reused." Apply the same wording to the Introduction. The fix is free.
- Demote "Because" in the abstract to a description of the protocol. That is free.
- The subtract-and-resend IBLT run is owed, at the cost of a laila run, before any causal claim. Together with F1 point
  2, both D/8 figures are upper bounds under the protocols as run.

## F9: the 9.5% bound and "no wrong difference" reproduce. UPHELD

1 − 0.05^(1/30) = 0.0950. `missing_total` and `spurious_total` are 0, every failed IBLT listing had 0 spurious, and
`undetected_failures` is 0. This complies with `power.md`.

## F10: "almost all from unplaced scaffolds" is unmeasured. UPHELD

The only support is prose in the results file (`:40-42`, `:192-194`), and no artifact assigns k-mers to records. The
counts and the KMC independence are fine.

**Fix.** Either drop the clause or state what is measured: record and base counts from `sets.tsv`. That is free. The
per-record tally is owed: a rescan on laila, which is heavy locally.

## F11: "2.4 times" is a cross-program ratio. NARROWED

**Diagnosis upheld.** 24.67 (`bench_e`) / 10.24 (`bench_kept`) = 2.41, across two programs.

**Scope.** Both programs carry the same anchor, batched ntHash, and they agree on it within 0.5% (11.37 against 11.43).
Normalising each to its own ntHash gives 1.12 / 0.46 = 2.43. The 28% disagreement the attacker cites is a derived
per-kept-k-mer slope, which amplifies small differences. The raw sampler disagreement is at most 11% (17.12 against
15.43).

**Fix.** "About twice as fast" is OVERTURNED: it understates a number that survives normalisation. Endorse instead
"2.4 times (across our two benchmark programs, which agree on batched ntHash within 0.5 percent)". That is free. Timing
Eisenstein every k-mer inside `bench_kept` is the clean fix, and it costs a small benchmark run.

## F12: Figure 4B matches the tables. UPHELD

`iblt_protocol.tsv` gives mean_rounds = max_rounds for all 8 rows, so the plotted mean is every repetition's value. The
30-seed byte range is 36,415,840-36,547,696.

---

## Minority record (rule 8)

- **F1's draft-row error.** The attacker raised it and then priced it anyway ("a one-way A→B message might suffice
  there", yet 34.60 sits in the table). I record that the attacker's own hedge was the correct ruling.
- **The asymmetric standard between F1 and F8.** No attacker raised it. I record it here: the same "smarter protocol"
  argument cuts both ways.

## Verdict table

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | GigaDB as run used two sketches, so 11.3 is per direction (22.5 total; free fix endorsed); the draft 34.6 fix is OVERTURNED because one-way decode is identical there (all differences on B's side); "reverses on every row" is a property of the protocol as run, not the sketch |
| F2 | UPHELD | separation commit 8c34af9 is an ancestor of run code e76fad1, and run.log records separate seeds; "used fixed seeds" is also false for the 30 OS-seed runs; deletion plus the caveat is free |
| F3 | UPHELD | k = 63 hashes are C, the 20,000-k-mer check is out of scope, the run used m/2 (2^10) not 0.8m, and the decode used one process; corrected sentence free; the 2^9 or kp10c run is owed |
| F4 | UPHELD | 1,263 failures, 2,526 buckets sent; bytes right; caption and ledger row 178 fix free; C4's attribution stands |
| F5 | NARROWED | the syncmer row is the 5-trial prior run and "507,219 of 507,221" is omitted (caption fix free; rerun owed); the link to "neither method returned a wrong difference" is speculation |
| F6 | UPHELD | IBLT seeds 1001-1030, each reused 10 times; wording fix free |
| F7 | UPHELD | decode adds about 50 Horner multiplications per element of A on GigaDB, and the CPU times agree; fix free |
| F8 | UPHELD | failed listings were large and error-free, so "must be resent whole" and the causal "because" overreach; the fix's "in km-peeler" is NARROWED to "in the protocol we ran" (free); IBLT subtract-and-resend run owed |
| F9 | UPHELD | bound 0.0950 and zero wrong differences reproduce |
| F10 | UPHELD | "almost all from unplaced scaffolds" has no measurement; drop it or state record counts (free); rescan owed |
| F11 | NARROWED | 2.4 is cross-program but survives ntHash normalisation (2.43); "about twice" OVERTURNED; disclosure wording free |
| F12 | UPHELD | Figure 4B values and caption match the raw tables |
