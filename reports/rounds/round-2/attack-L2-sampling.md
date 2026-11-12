# Round 2 attack, lens L2: sampling

Attacker: sampling lens. Frozen text: `reports/rounds/round-2/ATTACKED-paper.md` and `ATTACKED-ledger.md`.
Raw outputs: `results/kp10_vs_shibuya_2026-11-11/laila_b69c5c76/` (abbreviated `L/`).

**What I checked first.**
- **Prior rounds.** I ran `scripts/round history` on `power-of-two`, `bucket count`, `sawtooth`, `0.8 m`, `IBLT seed`,
  `Figure 4` (no prior mention of any), and on `partition 0`, `calibration`, `independent`, `1001` and `20,000`.
  - Round 1 already upheld the point that KP-10's (2^16, m = 64) configuration was calibrated on partition 0 of this
    same GigaDB pair: round 1's attack-L2 F7, upheld by adjudication-L2 F7, lines 147-148. F1 below builds on that
    point and says what is new.
  - No round has examined the IBLT seeds, the position of D on the bucket grid, or the sampling grid of Figure 4B.
- **Method cards.** No appendix was declared, so I looked for them. They are in `../library/methods/`, at library
  commit 2601816. Two apply here:
  - `pseudo_replication.md`: "your effective sample size is the number of independent decisions, not the number of
    rows". It governs F3.
  - `heaping.md:63-65`: "a single seed's draw is not a result". It governs the single-run rows in F2 and F4.
- **Literature.** The km-peeler paper is in `../library/papers_md/shibuya_2022_genomic-reconciliation.md`, which I read
  in full text for F6.
- **Local reruns.** All were tiny, run under `nice -n 19 ionice -c3`, with scripts in my scratchpad:
  - `kp10_bytes_sim.py`: KP-10's byte count simulated from D and the bucket count alone.
  - `one_round_m.py`: Poisson bucket loads.
  - `runs_test.py`: a Wald-Wolfowitz permutation test.

## F1 -- 11.3 vs 20.5 is one point on a sawtooth, and the GigaDB D sits 3.5% below its cliff

**Claims:** C1, C6, and the C13 inference.

**Defect.** The bytes KP-10 spends at the true D are a deterministic function of where D falls between two powers of
two. Whole-genome data, the 30 seeds and the pair add nothing to that number beyond D itself. The GigaDB D sits at
the most favourable point of that function. The abstract reports this one point as the method's size against the
table.

**Evidence.**
- **The sizing rule.** `recon/kp10c.c:541` sets `while (dest / 2^bits > 0.8*m) bits++`, so every bucket starts
  with a load between 25.6 and 51.2. A bucket costs (m + c + 1) × 8 = 536 B (`kp10c.c:574`). One round therefore
  costs 536 / load bytes per difference: 20.9 at the bottom of an octave and 10.5 at the top.
- **Where GigaDB sits.** D/2^16 = 3,238,101 / 65,536 = **49.41**, or 96.5% of the 51.2 cap.
- **D and the bucket count reproduce every measured KP-10 row.** `kp10_bytes_sim.py` draws D uniform values, buckets
  them, and splits while a bucket's count exceeds m. Its results against the measured rows:

| run | simulated B/diff | measured B/diff | simulated splits | measured splits | measured source |
|---|---|---|---|---|---|
| giga, D_est = D | 11.25 | 11.27 | 1,226 | 1,263 (mean 1,254) | `kp10_dest.tsv:6`, `kp10_seeds.tsv:2` |
| giga, D_est = D/8 | 20.74 | 20.76 | 58,565 | 58,602 | `kp10_dest.tsv:2` |
| giga, D_est = D/4 | 19.39 | 19.40 | – | – | `kp10_dest.tsv:3` |
| giga, D_est = D/2 | 16.69 | 16.69 | – | – | `kp10_dest.tsv:4` |
| giga, D_est = 2D | 21.70 | 21.70 | – | – | `kp10_dest.tsv:7` |
| draft, D_est = D | 17.30 | 17.30 | – | – | `kp10_dest.tsv:11` |

  The 30 OS seeds span only 36,415,840 to 36,547,696 bytes, ±0.18% (`kp10_seeds.tsv:2`). The seeds sample hash
  noise, not the method's size.
- **The same rule at other D.** With the same rule and the same m, the simulated cost at other first-round loads is:

| first-round load | 26 | 30 | 35 | 40 | 45 | 49.4 (GigaDB) | 51 |
|---|---|---|---|---|---|---|---|
| B/diff | 20.62 | 17.87 | 15.31 | 13.40 | 11.97 | 11.27 | 11.20 |

  Over a log-uniform D within one octave, the median is **14.8** B/diff and the maximum **20.9**. The maximum exceeds
  the table's 20.53.
- **The draft pair shows the same effect.** Its load is 7,931 / 256 = 31.0, so 17.3 B/diff (C6). Its 1.19× ratio to
  the table, against GigaDB's 1.82×, is set by D's position in the octave and not by anything about the pair.
- **The cliff.** Any estimate above 1.036 D (51.2 / 49.41) moves GigaDB to 2^17 buckets, at **21.70** B/diff.
  `kp10_dest.tsv:7` shows those bytes. The table at the same 4% overestimate has n = 1.092 D, about 21.35 B/diff
  (19.55 B per unit of n, `iblt_trials.tsv:7-11`), and is cheaper.
- **What is new since round 1.** Round 1 established that (2^16, m = 64) was chosen as the cheapest of four
  calibrations on partition 0 of this pair (adjudication-L2 F7). The "0.8 m" rule is described as reproducing "the
  2^16 and 2^8 buckets of the 2026-11-10 runs at the true D; it is fixed before any run here" (`kp10c.c:30-31`). So
  the rule was fitted afterwards to a configuration already tuned on this pair. A 0.75 m rule (48 per bucket) would
  have given GigaDB 2^17 buckets and 21.7 B/diff. The paper's Methods present the rule without this history.
- **Under the sampling card.** By `pseudo_replication.md`, the independent unit for "how big is the sketch against
  the table" is the pair, or rather its D. There are two (n = 2), and they give 1.82× and 1.19×. The "30 of 30"
  replicates the seed, which barely moves the bytes.

**Correction.**
- Abstract: "...it recovered every difference in 30 of 30 runs with 11.3 bytes per difference over two rounds, against
  20.5 for a one-round invertible Bloom lookup table; because buckets come in powers of two, the sketch's cost at the
  true difference ranges from about 11 to 21 bytes depending on the difference's size, and this pair sits near the
  low end."
- Results: add that the sizing rule was set to reproduce a bucket count calibrated on this pair, and that a 4%
  overestimate of D doubles the sketch to 21.7.

## F2 -- "its smaller size comes from the second round" confounds the round with the power-of-two grid

**Claims:** C4, and the inference in C13.

**Defect.** The one-round comparison point is 2^17 buckets at m = 64. That doubles the bucket count, and the
attribution changes two things at once. Holding the bucket count at 2^16 and raising m sizes the sketch for one round
at far less than 21.7.

**Evidence.**
- `one_round_m.py` takes Poisson bucket loads with λ = 49.41. It is calibrated on this pair: it predicts 1,253
  overflowing buckets at m = 64, and the measured mean is 1,254. At 2^16 buckets with c = 2:

| m | expected overflowing buckets | P(any overflow) | B/diff, one round |
|---|---|---|---|
| 84 | 0.18 | 0.16 | 14.09 |
| 88 | 0.017 | 0.017 | **14.73** |
| 90 | 0.005 | 0.005 | 15.06 |

- `kp10c` accepts `--m` up to MAXM = 256 (`kp10c.c:51,515`).
- **The split.** A one-round sketch costs about 14.7 B/diff, against the table's 20.5. Of the 9.2 bytes the sketch
  saves against the table, about 5.8 (63%) are saved in one round, and the second round adds 3.4 (37%). "Comes from
  the second round" is therefore false as an attribution.
- **What the 21.7 is.** One fixed-seed run (`kp10_dest.tsv:7`). Its bytes are deterministic, because at a load of 24.7
  no bucket overflows. Single-run variance is not the problem here. Choosing a doubled bucket count is.

**Owed.** One `kp10c --m 88 --dest 3238101` run on giga, 30 OS seeds, on laila. I did not run it locally.

**Correction.** "Sized for a single round by doubling the buckets, the sketch needed 21.7. With the bucket count held
and m raised to 88, a one-round sketch is expected to need about 14.7 (modelled, not run), so the second round
accounts for roughly 3.4 of the 9.2 bytes saved."

## F3 -- km-peeler's 30 "independent random seeds" share two of three hash functions between neighbours

**Claims:** C5, C7, and C15 (ledger rows 152 and 179).

**Defect.** km-peeler derives its r hash functions from a root seed as root + i. The harness used consecutive root
seeds 1001 to 1030. Trial t therefore uses MurmurHash3 seeds {1000+t, 1001+t, 1002+t}, and trial t+1 uses
{1001+t, 1002+t, 1003+t}. Adjacent trials share two of their three hash functions. All 30 trials together draw on only
32 distinct functions. The trials are not independent, so the Clopper–Pearson 9.5% bound does not hold for the
table.

**Evidence.**
- **How km-peeler seeds its hashes.** `bench/data/kmer_poly_work/km-peeler`, at HEAD b8ec39b (the pinned commit;
  `KMP_REPO` in `recon/kp10_vs_iblt_launch.sh:16`):
  - `ibflib.c:27-29` sets `j = root_seed + i; /*FIXME use a better method?*/`.
  - `ibflib.c:67` computes `MurmurHash3_x64_128(..., res[i].seed, ...)`.
  - `ibflib.c:331` places a key at `hash.ls64b % chunk_size + j * chunk_size`.
- **How the harness seeds trials.** `recon/kp10_vs_iblt_job.sh:176` passes seed `$((1000 + t))` for t = 1..30, with
  the same 30 seeds at every ratio and on both pairs.
- **The outcomes cluster.** Success by seed 1001 to 1030, from `run.log` `@@ iblt_trial`:
  - draft at n = D: `100000000000111111000001000000`, 8 of 30. Six consecutive successes at seeds 1013-1018 give 6
    runs where independence predicts 12.7. The permutation P(runs ≤ 6) is **0.002** (`runs_test.py`, 200,000
    shuffles).
  - giga at n = D: 16 runs, P = 0.66, so no clustering is visible at that scale.
- **The KP-10 side is clean.** Its 30 seeds are OS-random, and all 60 (x, λ) pairs in `run.log` are distinct. The
  9.5% bound holds for KP-10's 30/30.
- **What the 30/30 rows actually are.** The IBLT's "30 of 30" at 1.05 D comes from 32 functions arranged in
  overlapping triples. By `pseudo_replication.md`, its effective n is unknown and below 30. The selection compounds
  this: 1.05 was chosen as the smallest grid ratio passing 30/30 on these same seeds, at 0.95 / 1.0 / 1.05 / 1.1 /
  1.25.

**Correction.**
- C7: "...each method ran 30 times; the sketch with independent OS-random seeds, the table with consecutive root seeds,
  which km-peeler turns into overlapping hash triples."
- C5 and ledger rows 152 and 179: restrict the 9.5% bound to the sketch.
- Fix to discharge it: rerun `iblt_sweep` with root seeds spaced at least r = 3 apart, or drawn at random. This is
  cheap on laila.

## F4 -- Figure 4B's estimate grid lands on the table's best case, and its lines interpolate step functions

**Claims:** C9, and C2 (which survives).

**Defect.** Both methods' costs are step functions of D_est/D. The plotted grid {1/8, 1/4, 1/2, 1} hits the same phase
of both:
- **The table.** Doubling from 1.05 · D/2^j ends at n = 1.05 D exactly: the cheapest final round the protocol can
  reach.
- **The sketch.** Each grid point sits at 96.5% of its plateau's right edge, where its cost is the plateau's constant
  value.

Straight lines between the points imply intermediate costs that exist for neither method. The measured sketch points
at 2D and 4D (21.70 and 43.39, `kp10_dest.tsv:7-8`) are dropped, although the x-axis runs to 1.5
(`make_figures.py:250, 304`).

**Evidence.** Off-grid costs, computed from the measured per-n cost of 19.55 B and the measured success
thresholds (0/30 at 0.95 D, 30/30 at 1.05 D):

| D_est / D | table protocol | table B/diff | sketch | sketch B/diff |
|---|---|---|---|---|
| 0.9 | n = 0.945 D fails, then 1.89 D | about **55.4** | still 2^16 | 11.27 |
| 1/6 | 0.175, 0.35, 0.7, 1.4 D | about **51.3** | 2^14, same as D/4 | 19.40 |

The sketch's advantage under underestimation is therefore larger off the grid than on it.

**C2 survives this lens.** The D/8 point favours the table, and the sketch's single fixed-seed run is not a weakness
here: the simulation reproduces 20.74 against 20.76, and seed spread at D is ±0.18%.

**Correction.**
- Figure 4B caption: "Estimates are powers of two, where the table's final resend is exactly 1.05 D, its cheapest
  case. Both costs are step functions, and the lines only join the measured points."
- Draw the curves as steps, or as markers only.
- Show the sketch's 2D point, or state that a 4% overestimate doubles the sketch.

## F5 -- the 20,000-k-mer check that vouches for the k = 63 run sampled neither k = 63 nor its code path

**Claim:** C3.

**Defect.** The paper says "The k = 63 run used a numpy prototype, whose hashes match the C code on 20,000 sampled
k-mers." That sample cannot vouch for the k = 63 run, for three reasons:
- **It is not k = 63.** The sample comes from `recon/kp7e_hash.py:85-97` `self_test`. It draws uniform random 2-bit
  codes at `K` = 31 (`prototype/kmers.py:9`), not genome k-mers, and checks the numpy port against the C
  `kp7e_direct`.
- **It is not the k = 63 hash path.** The k = 63 run hashed with C, `hash/kp7e_stream.c`
  (`results/round1_owed_2026-11-11.md:81-83`). Only the decode is numpy, in `recon61`.
- **The real check is elsewhere.** The truth is hashed by the C reference (`recon/aswritten.py:225-226`). The actual
  guarantee is that check, element for element over all 21,382 values, not the 20,000 sample.

The result survives. The sentence names the wrong sample.

**Correction.** "The k = 63 run hashed with the C stream hasher and decoded in a numpy prototype, with buckets sized
from the true difference. Its recovered values equal the hashes of KMC's canonical 63-mer difference, element by
element."

## F6 -- "61 percent against about 3 percent" compares two different populations

**Claim:** C8, and round 1's carried item on harness error.

**Defect.**
- **Where the authors' 3% comes from.** Their Table 1 has three datasets: 50 SARS-CoV-2 genomes (mean diff 325.29,
  err 11.03), and random sequences mutated at pm = 0.001 (1,708.78 / 60.03) and 0.01 (15,342.81 / 357.57)
  (`shibuya_2022_genomic-reconciliation.md:747-782`). These are small or repeat-free genomes.
- **What drives spurious k-mers.** Their Section 3.5 says spurious k-mers arise when "an extended syncmer can contain
  a k-mer which belongs to another extended syncmer common to both datasets". In "the extreme case where each k-mer
  occurs once, our method computes exact k-mer set differences" (lines 378-383).
- **This pair is the opposite case.** The GigaDB difference is the unplaced scaffolds, which "share much sequence with
  the chromosomes" (`results/kp10_vs_shibuya_2026-11-10_prior.md:106-107`). The authors' own mechanism predicts a
  high spurious rate here.
- **What the contrast does.** It reads as a sign of misuse or a harness error. The population, not the harness,
  explains it. It is also 5 trials, sized from the true extended-syncmer difference (`prior.md:100`).

**Correction.** "...returned spurious k-mers numbering 61 percent of the true difference, against about 3 percent in
the authors' tests on viral and random sequences; their method is exact only when k-mers are unique, and this
difference consists of scaffolds that share sequence with the chromosomes."

## F7 -- "almost all from unplaced scaffolds" has no measurement behind it

**Claim:** C10.

**Defect.**
- **No count exists.** No file attributes the 3,238,101 differing k-mers to scaffolds. I grepped
  `results/kp10_vs_shibuya_2026-11-11.md:42,193`, `prior.md:27,106`, `round1_owed_2026-11-11.md` and the ledger. All
  of them assert it, and none counts it.
- **What the counts do show.** `L/tables/sets.tsv:2-3` give 147 records and 914,397,426 bases against 10 records and
  904,080,375 bases. That fits dropped scaffolds (about 10.3 Mb). It does not show how many of the 3.24 M differing
  k-mers lie on them rather than on edited chromosomes.

**Owed.** Count A∖B k-mers by source record: rescan Cannbio-2's 137 non-chromosome records against the GigaDB set.
This is a moderate job on laila, and I did not run it locally.

**Correction.** Either report that count, or write "consistent with the unplaced scaffolds (about 10 Mb) that the
release drops".

## Claim-list steering note

C1's refutation conditions list the numbers, the byte accounting and fairness to the table. They do not name the one
variable that sets the sketch's bytes: where D falls on the power-of-two bucket grid (F1).

C7 asserts "independent random seeds" for both methods as settled fact. The IBLT side is where it fails (F3). Neither
omission looks deliberate, but both point an attacker away from the weakest point.

## Verdicts

| claim | verdict | finding |
|---|---|---|
| C1 | NARROWED | F1, F3: the 11.3 is near the best point of a sawtooth (median 14.8, max 20.9 over an octave), and the IBLT's 30/30 rests on dependent seeds |
| C2 | SURVIVED | F4: the grid favours the table; single-seed bytes are reproduced by simulation to ±0.1 |
| C3 | NARROWED | F5: the 20,000-k-mer check is at k = 31 on another code path |
| C4 | REFUTED | F2: one round at 2^16 with m = 88 is about 14.7 B/diff; the 21.7 is the doubled grid, not "one round" |
| C5 | REFUTED | F3: the IBLT's 30 trials share hash functions; runs test P = 0.002 on draft at n = D; the bound holds for KP-10 only |
| C6 | NARROWED | F1: 17.3 is the draft D's octave position (load 31) |
| C7 | REFUTED | F3: km-peeler seeds are 1001-1030, hashes root + i, not independent |
| C8 | NARROWED | F6: 61% vs 3% compares repeat-sharing plant scaffolds with viral and random sequences |
| C9 | NARROWED | F4: the grid hits the table's minimum, lines interpolate step functions, and the 2D and 4D sketch points are omitted |
| C10 | NARROWED | F7: "almost all from unplaced scaffolds" is unmeasured |
| C11 | UNCHALLENGED | not reached under this lens |
| C12 | UNCHALLENGED | not reached under this lens |
| C13 | NARROWED | F2: about 63% of the size advantage (5.8 of 9.2 B/diff) is available in one round, not from resending |
| C14 | UNCHALLENGED | one genome, medians of 5; nothing concrete found under this lens |
| C15 | NARROWED | F3: ledger rows 152 and 179 extend the 9.5% bound to the IBLT; rows 153 and 177 call the doubled grid "sized for one round" |
