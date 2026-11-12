# Round 2 adjudication of L1 (correctness)

Adjudicator: skeptic of the attack. Attack file `reports/rounds/round-2/attack-L1-correctness.md`, sha256
`a1b5b20e5614...` (matches the sealed digest). Raw outputs `results/kp10_vs_shibuya_2026-11-11/laila_b69c5c76/` (`L/`).
All local checks were tiny and ran under `nice -n 19` (table arithmetic, a Poisson tail, a 200k-shuffle permutation
test on 30 booleans, reading JSON). Nothing was hashed and nothing ran on laila.

**General defect in the attack: its AP line numbers are off by 16.** Every `AP:n` above line 16 points 16 lines early
in the sealed `ATTACKED-paper.md`. AP:119 is line 135, AP:51 is 67, AP:47 is 63, AP:130 is 146, AP:134 is 150 and
AP:136 is 152. The quoted text is right in every case, so no finding fails on this. The runner should cite the sealed
file's own line numbers when applying any repair.

**Method cards.** No `refs/methods/` directory exists in this repo or in `../library`, so no card governs F4's
independence question. Rule 6 therefore has nothing to apply. F4 is decided on the code.

**History (rule 5).** `scripts/round history` finds no earlier ruling on "every difference", "independent random seeds",
D/8, km-peeler seed overlap, the x/λ sentence or the k = 63 provenance. One earlier ruling bears on F1: round 1's
L5 F5 (adjudication-L5.md:95-112). It ruled that the byte count does not depend on who decodes, and that the two
methods return different things. F1 does not re-file that point. It is a new consequence that follows once the
decoder is two-party.

## F1: "every difference" needed both sketches. NARROWED, as filed. Fix endorsed (free)

**Premise checked.** `recon/kp10c.c` `decode` (l.257-307) builds `S = SA - SB`. It flags A's own roots of P and B's
own roots of Q. A bucket is accepted only when both flag counts and both c-sum checks pass. `reconcile` then pools
`oa` and `ob`. The byte count is `R.entries * 8`, which is one side's entries per round.

On giga, A∖B = 3,237,326 and B∖A = 775 (`truth.tsv:2`). Party A can turn into k-mers only the values it holds, so the
775 k-mers that A *lacks* never reach A in any form other than the roots of an unfactored Q. The scored 30/30 "missing
0" is therefore the union of two recoveries, each made from the other party's sketch. By contrast, a single table
listing gives one party both sides' k-mers (`iblt_trial` lines, `listed_A_minus_B` and `listed_B_minus_A`). I
recomputed 2 × 36,481,232 / 3,238,101 = 22.53.

**Scope.** The like-for-like test the claim sets itself ("one side's message") holds for the bytes, but not for the
output that "every difference" names. The paper body already concedes this (line 146, "Each party of the sketch
recovers only its own differing k-mers"). The abstract does not. On the draft pair, A∖B = 0, so one message A→B
delivers everything, and C6 is unaffected (I agree with the attacker there).

Two small items go uncounted in the measured exchange. Acceptance uses both parties' checks, and the round-2 split
set is decided jointly. In a real exchange that costs one accept bit per bucket, 8 KB on giga (0.0025 B/diff), which
is negligible. The paper should not claim that "each party checks with its own data only" yields a decision without
any return message. It needs one bit per bucket.

**Fix.** Endorsed: the attacker's first wording, "each party recovered its own differing k-mers in 30 of 30 runs from
the other's sketch of 11.3 bytes per difference over two rounds, against 20.5 for one table from which one party lists
both sides'". The alternative, "22.5 bytes per difference for the two sketches exchanged", is also correct. Either fix
is free. Not endorsed: the unbuilt "optimised return message" as anything more than a remark. It was never run.

## F2: the D/8 underestimate taxed both methods alike. UPHELD on the diagnosis. The fix is NARROWED (free)

**Premise checked against `kp10_dest.tsv` and `iblt_protocol.tsv`.**
- The splits add up: 8,192 + 16,384 + 34,026 = 58,602.
- So do the bytes: 536 × 8,192 + 536 × 16,384 + 54,039,520 = 67,212,256, which is exactly the D/8 row.
- The fixed seeds (0x5eed/0x1a3d) are shared by all sweep rows, so once rounds 1-2 overflow everywhere (loads of about
  395 and 198 against m = 64), the D/8 run replays the D/2 run.
- The degradation factors reproduce: sketch 1.48 / 1.72 / 1.84, table 1.50 / 1.75 / 1.875.
- At D/8 on giga, "resent alone" was vacuous for the first two rounds. The 20.8 vs 38.5 gap is the true-D gap carried
  forward.
- The abstract's "Because an overflowing bucket is resent alone, ..." asserts a cause that the GigaDB sweep does not
  isolate. The same goes for the Results sentence (line 146) and Figure 4's title, "Resending only overfull buckets
  keeps the sketch small when the difference is underestimated". That title captions a GigaDB-only panel 4B, so the
  attacker's C9 "SURVIVED" should carry this repair too.

**Where the attacker overreaches.**
1. The proposed text "so each round resent the whole sketch at double resolution" is false for the last round. In the
   D/8 run the final round sent about 2,516 halves, not 131,072 buckets: its D/2 replay ends exactly as the true-D run
   does, with about 1,260 failed buckets. That would install a new error. Write instead: "An eightfold underestimate
   overfilled every bucket in the first two rounds, so each was resent whole at double resolution; the cost rose 1.84-fold
   to 20.8, against 1.87-fold to 38.5 for the table resent at double size." (free)
2. "Rides on power-of-two granularity" understates the draft pair. The selective resend did operate there: draft D/2
   (128 buckets at a load of 62, 46 split) cost 14.87 B/diff, *below* the one-round true-D 17.30. Draft D/8 rose only
   1.23-fold against the table's 1.88-fold.
   - What survives is a causal statement on the **draft** pair only: "on the draft pair, resending only the failed
     buckets held the cost to 1.23 times its true-D value, against 1.88 for the table".
   - It needs the caveat that the true-D sizing rule left that pair at a load of 31. (free)

## F3: "needed 21.7" and "resent about 1,250". NARROWED. Fixes endorsed (free); the m = 88 run is optional and costs a run

**(a)** 70,254,592 B is exactly 2 × `bytes_round1` 35,127,296 (`kp10_dest.tsv:7`, `kp10_seeds.tsv:2`). I reproduced the
attacker's Poisson calibration. With 2^16 buckets, m = 64 and D = 3,238,101, the model gives 1,253.4 expected overfull
buckets against the observed mean of 1,254. The m = 88 projection gives 0.017 overfull at 14.73 B/diff.

What the attack shows is that "needed" overstates a measurement made at one power-of-two step. What it does not show
is that the attribution is false *within the paper's stated design* (m = 64 fixed, power-of-two ranges). There, 2^17 is
the only one-round option, so the sentence is a design-bound statement that the paper does not label as one.

The attacker's split of the gap (4-6 B/diff of one-round efficiency, 3-5 of second round) is a projection. It must not
enter the paper.

Fix endorsed: "With twice the buckets, enough for one round, it used 21.7; finer sizing was not tried." (free; it drops
the causal "so ... comes from the second round"). The `kp10c --m 88 --dest D` giga run is **owed only if** the paper
wants any attribution sentence back. It costs one laila run (hashing both genomes, about 6 minutes of CPU).

**(b)** `run.log:26-27` shows 1,263 failed and 2,526 sent in the fixed-seed run. "Resent about 1,250 of them" (of the
2^16 buckets) is ambiguous rather than false, because each of those buckets was resent as two halves. The fix "split
about 1,250 overfull buckets and sent their halves (about 2,500)" is endorsed and free.

## F4: km-peeler trials share hash functions. UPHELD. Wording fixes endorsed (free); the spaced-seed rerun costs a run

**Premise checked** at km-peeler b8ec39b:
- `ibflib.c:28` has `j = root_seed + i`, and `ibflib.c:331` has the per-block position.
- The harness gives trial t the root seed 1000 + t (`kp10_vs_iblt_job.sh`, `iblt_sweep`).
- Consecutive trials therefore share two of three key-to-cell maps, up to block offset, and the 30 trials span the 32
  Murmur seeds 1001-1032.
- "Independent random seeds" (line 135) is false for the table on both words: the seeds are fixed, and they overlap.

**Data check.**
- I reproduced the patterns: draft n = D `100000000000111111000001000000`, giga n = D 18/30.
- On draft, the longest run of successes (6) has p = 0.0010 by permutation, which matches the attacker.
- The attacker's "only three runs" is a miscount. There are six runs, three of them successes, and runs ≤ 6 has
  p = 0.0019.
- Giga shows nothing (p ≈ 0.66).
- With two strata tested, the code fact carries the finding, not the p-value. No p-value should be published.

**Scope limits the attacker did not state.**
1. The 10-rep protocol runs behind Figure 4B and C2 use seeds 100000k + 1000r + 1000f (`iblt_proto`), at least 125
   apart. They share no hash function, so **F4 does not reach C2 or C9**.
2. Choosing 1.05 as the smallest ratio that went 30/30 helps the table's byte count. It is not unfair *to* the table,
   but it does make the table's 9.5% a post-selected bound.

The sketch's 30 runs are OS-random (`run.log`: 60 `"seed_x": "os"` lines, 30 per pair), so the 9.5% stands for the
sketch.

**Fixes.** Endorsed and free:
- The attacker's Methods wording: "the table 30 times with root seeds 1001 to 1030, which km-peeler expands to
  overlapping hash functions".
- Restricting the 9.5% bound to the sketch.

The spaced-seed rerun (root seed 1000·t) is optional. It costs a laila run (30 trials × 5 ratios × 2 pairs of whole-set
builds) and could move the 30/30 ratio. It is not needed once the wording is fixed.

## F5: the x/λ sentence describes the superseded run. UPHELD. The fix is NARROWED (free)

**Premise checked.**
- The fixed-seed GigaDB C runs log x = 492247986318507 and λ = 834835136232447260. Those are the same values as the
  post-separation k = 63 run (`run_..._k63_sym_bauto_r12.json`), which was drawn from the domain-tagged streams of
  `results/round1_owed_2026-11-11.md` §1.
- The 30/30 runs used OS randomness for both.
- "The GigaDB sketch ran before x and λ were separated" is false for every GigaDB number in the paper.
- "Our experiments used fixed seeds ... so they test the construction rather than the bound" is false for the
  reconciliation's 30-run results. Those runs are in F_q, where OS mode meets the bound's premise.

**Repair narrowed.** The sentence at line 67 also covers the paper's other experiments (speed, sampling, spaced seeds).
I did not check their seeds, and the attacker did not either. Replacing the whole sentence could make the
fixed-seeds statement disappear where it is true. Write instead:

> "The 30 reconciliation runs per pair drew x and λ independently from the operating system; the size sweep and timed
> runs used one fixed pair of seeds."

Then keep "our other experiments used fixed seeds ..." only for the experiments it is true of. Delete the "before x and
λ were separated" clause. (free)

## F6: the k = 63 provenance sentence is wrong. UPHELD. Fix endorsed (free)

**Premise checked.**
- `recon/aswritten.py` hashes with `hash/kp7e_stream.c` and computes the truth with C `kp7e_direct`. numpy only
  reconciles.
- `kp7e_hash.py` `self_test` (l.85-98) runs at K = 31 on random 2-bit codes, which is not possible at k = 63. The
  ledger source is `prior.md:33`.

**Fix.** The attacker's wording "hashed with the same C code" is accurate. Both `kp10c` and `kp7e_stream` call
`kp7e_sym_seq` (kp7e_stream.c:6, 53), so I endorse the attacker's sentence as written (free). The hash-space note on
KMC agreement (`verify_k63_kmc.json`: KMC 63-mers hashed by C, compared with the recovered values) is accurate and
optional. Linearity and "needs no canonicalization" survive, as the attacker says.

## F7: the extended-syncmer listing was incomplete in 5/5 trials. NARROWED from REFUTED. Repair NARROWED (free)

**Facts confirmed** from `full_kmp_ext_giga.json`:
- All five trials have listed 507,219 of D_ext 507,221, `complete: false`, missing 7, spurious 1,986,807 and
  13,667,246 B.
- `complete` is `listed == D_ext` (`kmpeeler_run.py:179`).
- The row has 5 trials from the 2026-11-10 harness. The paper does not say so.

**Where the attacker overreaches.**
1. "'Which we did not trace' is contradicted by its own JSON" does not hold. The JSON has counts only. It does not
   identify the 7 missed k-mers, does not tie them to the 2 unlisted syncmers, and gives no cause. "The 7 missed k-mers
   are traced to the 2 unlisted syncmers by the run's own record" is an inference, not a trace.
2. "A table that failed to list" is not established either. With identical outcomes across five seeds, the 2 syncmers
   were probably never in the table, rather than caught in a peeling stopping set. One candidate is km-peeler's
   `pack2bit` skipping any fragment with a non-ACGT base (`ibflib.c`, comment above l.331). That is a hypothesis I did
   not test, and it should not be printed. The seeds 2000-2004 also overlap per F4, which weakens "five seeds" as five
   independent draws but does not change the conclusion that the cause is deterministic.

**Fix.**
- Endorsed: Table 4 reads "4.2, listing incomplete in 5 of 5 trials (2 of 507,221 syncmers unlisted, identically for
  every seed)". (free)
- Endorsed: line 150 states the incomplete listing and that the row comes from 5 trials of the earlier harness. (free)
- Overturned: "a probable harness defect". That is an unmeasured cause. Write "cause not investigated". Keeping
  "which we did not trace" is then accurate.

## F8: "66 multiplications per element" is the build. NARROWED. The fix is NARROWED (free)

**Diagnosis upheld.**
- `psums` does M = 66 multiplications per element, which is the build cost.
- `flag_roots` evaluates a degree-eA polynomial at every own element of every bucket with a difference. On giga that is
  all 65,536 buckets (`run.log:26`), so side A pays about 49 more per element. Side B, with eB ≈ 0, pays almost
  nothing.
- "3 cell updates" leaves out the three MurmurHash3 calls (`ibflib.c:67`), so the two counts are in unlike units.

**Repair narrowed (rule 3, whole set).** The attacker's replacement quotes decode CPU (110.3 against 3.5 s) and omits
the build. In the same `timed.tsv` row, KP-10's whole one-side build (hash, dedup and power sums) is 89.4 CPU-s against
**204.8** for km-peeler's text and build. Printing only the recovery gap would install a one-sided cost picture.

Endorsed (free): "Building the sketch takes m + c = 66 field multiplications per element; recovery evaluates a polynomial
of degree up to m at each element of every bucket that differs, which on the GigaDB pair was every bucket. On one
thread, one side's build took 89 CPU-seconds against 205 for the table, and decoding took 110 against 3.5."

## F9: ledger rows. NARROWED

Under C15's test ("cites a source line that does not say what the row says"), **no row miscites its source**:
- AL:115 matches `prior.md:33`.
- AL:109 matches `prior.md:102`.
- AL:183 is true, as the attacker concedes.
- AL:159 "1 of 37 runs" is literal. The collision was at OS seed 25 (results file l.128), so it is also 1 of 31
  distinct draws.

Three of the attacker's rows (AL:115, 109, 183) are the paper defects F6, F7 and F4 restated. Their rows should take
the same caveats once the paper changes, which follows from those findings and is not separate. What survives
independently is AL:150 and AL:178, where "buckets resent" names the failed buckets while 2,526 halves were sent.
Relabel them "buckets split (halves resent)". (free) Adding "31 distinct (x, λ) draws" to AL:159 is optional and
free.

## Minority / corroboration (rule 8)

I read no other round-2 attack, as instructed, so I cannot say whether F2 or F4 was raised elsewhere. Both rest on a
reproduced premise: the D/8 byte identity, and `root_seed + i` in the pinned source. If they turn out to be single-
attacker findings, record them as a correct minority.

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | numbers reproduce; "every difference" pooled two recoveries each needing the other's sketch, the abstract must say each party recovers its own (free) |
| F2 | UPHELD | the giga D/8 sweep degrades both methods ~1.85x, so the causal "because ... resent alone" is unsupported; attacker's "each round resent the whole sketch" is false for the last round, use "every bucket overflowed in the first two rounds" (free); a draft-pair causal claim survives |
| F3 | NARROWED | "needed 21.7" is one power-of-two step and is true only within m = 64 power-of-two design; drop the attribution (free); "resent 1,250" ambiguous, relabel as split with ~2,500 halves (free); m = 88 run optional, costs a laila run |
| F4 | UPHELD | km-peeler root_seed+i with seeds 1001-1030 makes the table's trials overlap; "independent random" false for the table, 9.5% valid for the sketch only (free); does not reach the 10-rep protocol (Fig 4B); spaced-seed rerun optional, costs a run |
| F5 | UPHELD | the GigaDB C runs post-date the x/λ separation and the 30/30 used OS seeds; repair narrowed to the reconciliation runs so the sentence is not dropped where it may be true (free) |
| F6 | UPHELD | k = 63 hashed by C kp7e_stream (kp7e_sym_seq), 20,000 check is K = 31 packed codes; attacker's sentence endorsed (free) |
| F7 | NARROWED | 5/5 identical incomplete listings confirmed, disclose them and the 5-trial earlier harness (free); "did not trace" is not contradicted and "probable harness defect" is unmeasured, write "cause not investigated" |
| F8 | NARROWED | 66 is the build cost and the units differ; the fix must add build CPU (89 vs 205 s) beside decode (110 vs 3.5 s) or it installs a one-sided picture (free) |
| F9 | NARROWED | no row miscites its source; only AL:150/178 "buckets resent" need relabelling (free), the rest follow F4/F6/F7 |
