# KP-10 vs Shibuya et al.'s IBLT (km-peeler), whole genomes, one host (laila job b69c5c76857f93ac, 2026-11-11)

**Verdict.**
- **Exact recovery.** On two real near-identical cannabis assembly pairs, both methods recovered the exact canonical
  31-mer difference whenever their protocol finished. KMC 3.2.4 confirmed the truth element for element.
- **Size at the true D.**
  - KP-10 used **11.27 bytes per difference** on the GigaDB pair (D = 3,238,101; D/|A| = 0.73%). It always needed **2
    rounds**: about 1,254 of 65,536 buckets were resent.
  - km-peeler's IBLT needed **20.53 B/diff** to succeed in one round in 30 of 30 trials (n = 1.05 D).
  - The ratio is **1.82×**, the same as the earlier partition-subset run. It now holds on the whole genome, over 30
    random seeds per method, on one host.
- **The size advantage comes from the second round.**
  - Sized to finish in one round (2^17 buckets for the same D), KP-10 took 21.70 B/diff, about the IBLT's 20.53.
  - On the small-difference draft pair (D = 7,931), KP-10 took 17.30 B/diff in one round, against 20.54 for the IBLT.
- **When D is underestimated, KP-10 degrades more gently.**
  - KP-10 resends only the buckets that fail, so at D/8 it still finished exactly at 20.76 B/diff (5 rounds).
  - The IBLT fails whole: at n = 0.95 D it succeeded 0 times in 30, and at n = D 18 times in 30. Resending at double
    size from D/8 cost 38.49 B/diff over 4 rounds.
- **Compute is comparable end to end, single thread**, though the work splits differently.
  - KP-10's decode is about 30× slower than the IBLT's diff and list.
  - km-peeler's build is slower than KP-10's whole sketch, because it parses k-mer text one character at a time.

Raw outputs: `results/kp10_vs_shibuya_2026-11-11/laila_b69c5c76/` (abbreviated `L/` below).
- `run.log` holds every result as an `@@ <kind> <json>` line.
- `counters.tsv` holds pcount's wall time and task-clock per phase.
- `meta.tsv` records the host; `job.json` the Galaxy job.

Tables: `L/tables/*.tsv`, written by `recon/kp10_vs_iblt_summarize.py`. Each row's `src` column gives its run.log
line(s). The previous run (2026-11-10, numpy, IBLT on 4 of 64 partitions) is preserved verbatim in
`results/kp10_vs_shibuya_2026-11-10_prior.md`, the version of this file at 69825c3.

## 1. The pair, and why

| pair | \|A\| | \|B\| | \|A∖B\| | \|B∖A\| | D | D/\|A\| | source |
|---|---|---|---|---|---|---|---|
| **giga**: Cannbio-2 vs Cannbio2_GigaDB | 442,151,894 | 438,915,343 | 3,237,326 | 775 | **3,238,101** | 0.73% | `L/tables/truth.tsv:2` (run.log:17,19) |
| **draft**: Cannbio-2 vs Cannbio-2_draft | 442,151,894 | 442,159,825 | 0 | 7,931 | **7,931** | 0.0018% | `L/tables/truth.tsv:3` (run.log:18,20) |

- **Why these pairs.** Both are real assemblies of one cannabis genome on bench (`hub-build/2bit`), and both meet
  D/|A| ≲ 1%.
- **giga is the primary pair.** It is the largest real difference on bench under 1%, so it tests size. The difference
  is almost entirely the unplaced scaffolds that the GigaDB release drops.
- **draft is the secondary pair.** It tests the small-D regime, and its contigs are written in other orientations
  (45.36% of bases on flipped contigs, `results/round1_owed_2026-11-11.md`, check (c)).
- **Not used.** ANC-1 vs Salk_ANCa has D = 0 (prior run).
- **Truth.** Distinct canonical 31-mers from `recon/kmtools.c kset` agree with KMC 3.2.4 (`recon/kmc_diff.sh`, `-ci1`)
  in every set count and element for element in both differences (`L/tables/truth.tsv:2-3`, last columns;
  `L/tables/sets.tsv:2-4`).

## 2. Methods

- **Host.** laila: 2 × Xeon Gold 6330, 112 threads, kernel 6.1.0-52, gcc 14.2.0 in `docker.io/library/gcc:14.2.0`
  (`L/meta.tsv`). It ran as a UDT job through cannabis-genome `scripts/laila_perf.py`.
  - The nine phases ran one after another, each timed by pcount.
  - **Timed runs are single-threaded and ran alone within the job.** laila is shared, so other Galaxy jobs may have run
    at the same time.
- **Code.** Repo commit e76fad1 (`recon/kp10c.c`, `recon/kmtools.c`, `recon/kp10_vs_iblt_job.sh`,
  `recon/kp10_vs_iblt_launch.sh`). The summarizer and the fetcher were added in b051309 and fixed in this commit.
- **KP-10 (`recon/kp10c.c`), on k-mers as written.**
  - **Hash.** Every 31-mer is hashed by the Eisenstein KP-7 hash S = (F + λ)(R + λ) over F_q, q = 2^61 − 1
    (`hash/kp7e.c`, `kp7e_sym_seq`). There is no canonicalization and no 2-bit packing. The distinct values are each
    side's set.
  - **Sketch.** Buckets are ranges of the value (its top bits). Each bucket holds its size plus power sums S_1..S_66
    (m = 64, c = 2), at 8 bytes per entry. The byte count is one side's message, summed over rounds.
  - **Decode.** Newton's identities, then a Padé step, as in `recon/recon61.py`. **Each side checks itself with its own
    data only**: its own flagged roots, plus power sums of the other polynomial's roots computed from that
    polynomial's coefficients. This answers round 1's "oracle decoder" finding (adjudication L5 F5). A failed bucket
    is split and resent.
  - **Back to k-mers.** Each side rescans its genome for the k-mers whose S it flagged. Those k-mers are canonicalized
    only to score them against the truth.
  - **Sizing rule.** The initial bucket count is the smallest 2^b with D_est / 2^b ≤ 0.8 m, fixed before the run. At
    the true D it gives the 2^16 and 2^8 buckets of the earlier runs.
- **Shibuya et al. (km-peeler).** https://github.com/yhhshb/km-peeler at
  **b8ec39b789366f5f5f6d22121dc409c65900a408**, source unmodified.
  - **Build.** `compile_options.h` is written byte for byte as their `configure.py sequences -k 31` writes it. Their
    `make ibltseq` is run with `-fcommon` added, as gcc ≥ 10 requires.
  - **Parameters.** Their defaults: r = 3 hash functions, ε = 0. n = ratio × the true D, one IBLT file (one side's
    message) per party.
  - **Input.** The SAME canonical k-mer sets, streamed as text into `ibltseq build`. This is the generous choice: their
    own `ibltseq kmers` is strand-specific, and their `2set.py` dedup is Python.
  - **Scoring.** A trial succeeds when the listing, after cancelling k-mers listed on both sides, equals the truth.
    Such net-zero pairs appeared in local smoke tests, because km-peeler cells (`int64 counter` + `uint8 keysum[8]`,
    `ibflib.h`) carry no check hash. They appeared in **no** whole-genome trial
    (`L/tables/iblt_trials.tsv`, column `trials_with_both_side_listings`), so the raw and cancelled scores agree
    throughout.
- **Protocol decisions never read the truth.** KP-10 splits on its own failed checks. The IBLT protocol resends at
  double size when `ibltseq list` exits non-zero. The truth only scores the outcome.

**Commands** (from the repo root):

```
bash recon/kp10_vs_iblt_launch.sh laila            # inputs tarball + laila_perf run of the 9 phases below
#   build; setup; kp10_timed giga; iblt_timed giga; kp10_timed draft; iblt_timed draft;
#   kp10_sweep; iblt_sweep; iblt_protocol; cleanup                (recon/kp10_vs_iblt_job.sh <phase>)
python3 recon/kp10_vs_iblt_fetch.py b69c5c76857f93ac <dir> --cleanup   # the driver lost the tunnel mid-poll
python3 recon/kp10_vs_iblt_summarize.py L L/tables
bash recon/kp10_vs_iblt_launch.sh dry A.fa giga.fa draft.fa <dir>     # the same phases locally, on tiny inputs
```

## 3. Size, success and rounds at the true D

| pair | method | bytes (one side) | B/diff | rounds | exact | source |
|---|---|---|---|---|---|---|
| giga | KP-10, fixed seeds | 36,481,232 | **11.27** | 2 (1,263 buckets split) | yes | `L/tables/kp10_dest.tsv:6` (run.log:44); run.log:28 |
| giga | KP-10, 30 OS-random seeds | 36,415,840–36,547,696 | **11.26** mean | 2 in 30/30 | **30/30** | `L/tables/kp10_seeds.tsv:2` |
| giga | KP-10, round 1 alone | 35,127,296 | 10.85 | – | – | `L/tables/kp10_seeds.tsv:2` |
| giga | KP-10 sized for one round (2^17 buckets) | 70,254,592 | 21.70 | 1 | yes | `L/tables/kp10_dest.tsv:7` (run.log:45) |
| giga | IBLT n = 1.05 D | 66,477,002 | **20.53** | 1 | **30/30** | `L/tables/iblt_trials.tsv:9` |
| giga | IBLT n = 1.0 D | 63,311,450 | 19.55 | 1 | 18/30 | `L/tables/iblt_trials.tsv:8` |
| giga | IBLT n = 0.95 D | 60,145,898 | 18.57 | 1 | 0/30 | `L/tables/iblt_trials.tsv:7` |
| draft | KP-10, 30 OS-random seeds | 137,216 | **17.30** | 1 in 30/30 | **30/30** | `L/tables/kp10_seeds.tsv:3` |
| draft | IBLT n = 1.05 D | 162,938 | **20.54** | 1 | **30/30** | `L/tables/iblt_trials.tsv:3` |
| draft | IBLT n = 1.0 D | 155,162 | 19.56 | 1 | 8/30 | `L/tables/iblt_trials.tsv:2` |

- **Equal footing.** Each method has 30 independent trials on the whole genome, sized from the true D. 30/30 bounds the
  failure rate below **9.5%** (one-sided 95% Clopper–Pearson; `failure_upper95` columns) for both. That is not a
  measured equality of reliability, and the paper should not call it one.
- **What KP-10's 2 rounds mean.** On giga a mean of 1,254 buckets (1.9% of 65,536) failed their checks in round 1. They
  were split in two and resent (`kp10_seeds.tsv:2`, `mean_splits`), so the bytes above include round 2. The IBLT
  rows are one round.
- **The 1.82× gap is a two-round sketch against a one-round table.** With one round enforced by sizing, KP-10 took
  21.70 B/diff and the IBLT 20.53. The 21.70 is one fixed-seed run, not 30.
- **Undetected errors: none, on either side.**
  - KP-10: 0 missing and 0 spurious k-mers over 30 + 30 seeded runs (`kp10_seeds.tsv:2-3`).
  - IBLT: 0 trials exited 0 with a wrong listing (`iblt_trials.tsv`, `undetected_failures`).
  - KP-10's detection is still probabilistic: round 1's 140-difference counterexample to C16 stands. It was simply
    never hit here.
- **Hash collisions.** In 1 of 37 giga KP-10 runs (OS seed 25), two classes shared by both genomes collided. A and B
  each lost one distinct value, 442,151,893 and 438,915,342, and the difference was still recovered exactly
  (run.log:96; `L/tables/collisions.tsv:2`). The draft runs saw none (`collisions.tsv:3`). This fits the bound's
  allowance of a few colliding pairs per assembly (round-1 adjudication L3): a random function would give about
  0.04 per 4.4 × 10^8-element set per draw.

## 4. When D is underestimated (and overestimated)

KP-10 used fixed seeds and one run per D_est (`L/tables/kp10_dest.tsv`). The IBLT used the resend protocol: start at
n = 1.05 D_est, double n after each failed listing, 10 repetitions each (`L/tables/iblt_protocol.tsv`).

| pair | D_est / D | KP-10 rounds | KP-10 B/diff | IBLT rounds (all 10 reps) | IBLT B/diff (mean) | sources |
|---|---|---|---|---|---|---|
| giga | 1/8 | 5 | 20.76 | 4 | 38.49 | kp10_dest.tsv:2; iblt_protocol.tsv:6 |
| giga | 1/4 | 4 | 19.40 | 3 | 35.93 | kp10_dest.tsv:3; iblt_protocol.tsv:7 |
| giga | 1/2 | 3 | 16.69 | 2 | 30.79 | kp10_dest.tsv:4; iblt_protocol.tsv:8 |
| giga | 1 | 2 | 11.27 | 1 | 20.53 | kp10_dest.tsv:6; iblt_protocol.tsv:9 |
| giga | 2 | 1 | 21.70 | – | – | kp10_dest.tsv:7 |
| giga | 4 | 1 | 43.39 | – | – | kp10_dest.tsv:8 |
| draft | 1/8 | 4 | 21.36 | 4 | 38.55 | kp10_dest.tsv:5; iblt_protocol.tsv:2 |
| draft | 1/4 | 3 | 19.19 | 3 | 35.97 | kp10_dest.tsv:9; iblt_protocol.tsv:3 |
| draft | 1/2 | 2 | 14.87 | 2 | 30.82 | kp10_dest.tsv:10; iblt_protocol.tsv:4 |
| draft | 1 | 1 | 17.30 | 1 | 20.54 | kp10_dest.tsv:11; iblt_protocol.tsv:5 |
| draft | 2 | 1 | 34.60 | – | – | kp10_dest.tsv:12 |
| draft | 4 | 1 | 69.20 | – | – | kp10_dest.tsv:13 |

- **Both methods ended exact in every row.** KP-10: `exact` column, 0 unresolved buckets. IBLT: 10/10 cancelled-exact
  in every row.
- **KP-10.**
  - **Behaviour.** An underestimate overfills the buckets. Most fail their checks and are split, and only failed
    buckets are resent. Bytes therefore grow slowly: at D/8, 1.84× its own true-D cost on giga.
  - **More rounds.** It costs about one round per halving of the estimate.
  - **Overestimates.** These cost bytes in proportion: one round, but 2× or 4× the buckets.
- **IBLT.**
  - **Behaviour.** An underfilled table fails as a whole: 0/30 at 0.95 D, 18/30 at D on giga
    (`iblt_trials.tsv:7-8`). km-peeler's partial listing cannot be reused, so each retry resends a full, larger table.
  - **Cost.** 1.50× its own true-D cost at D/2 and 1.87× at D/8 on giga, and 1.8–2.1× KP-10's bytes at every
    underestimate here.
  - **Not tested.** A smarter IBLT protocol (subtract what was listed, then resend a small table) is not what km-peeler
    implements.

## 5. Wall and CPU time, same host, one thread each

| pair | step | KP-10 (kp10c) CPU s | km-peeler CPU s | sources |
|---|---|---|---|---|
| giga | per side, FASTA → message | 6.9 hash + 33.0 dedup + 49.5 power sums = **89.4** | 35.5 canonical set (kmtools kset, user) + **204.8** text→build = 240.3 | run.log:28; `L/tables/timed.tsv:2,4`; `L/tables/sets.tsv:2` |
| giga | decode | **110.3** (both rounds, both sides' root flagging) | **3.5** (diff + list) | `timed.tsv:2,4` |
| giga | back to k-mers | 36.7 (rescan both genomes) | – (the IBLT lists k-mers) | `timed.tsv:2` |
| giga | whole phase (pcount wall / task-clock) | 326.1 / 326.0 | 287.4 / 389.9 (excludes kset: 41.5 + 38.9 s) | `L/counters.tsv` run_2, run_3; `sets.tsv:2-3` |
| draft | whole phase (pcount wall / task-clock) | 271.5 / 271.4 | 197.7 / 277.0 (excludes kset: 41.5 + 41.6 s) | `L/counters.tsv` run_4, run_5; `sets.tsv:2,4` |

- **Totals.** From FASTA to recovered k-mers, KP-10 took about 326 CPU-s on giga. km-peeler's pipeline took about 460
  CPU-s: 390 for the phase, less a 10-second text-generation measurement that is not part of the pipeline, plus
  about 80 for the two canonical sets.
- **This is not a speed claim for either algorithm.**
  - km-peeler's build reads text through `fgetc`, with 41 s of system time per side.
  - KP-10's decode evaluates degree-~50 polynomials at every element of every differing bucket, single-threaded. With
    8 threads its whole reconcile took 27.5 s on giga (`kp10_dest.tsv:6`).
  - Per element, KP-10 still costs m + c = 66 field multiplications against r = 3 cell updates.
- **The C port changed the cost picture, not the result.** The numpy prototype took about an hour per pair. The C
  kp10c reproduces its sizes to within 0.04% (36,481,232 vs 36,494,096 B there: 1,263 vs 1,275 splits under different
  seeds, canonical-code buckets there and value-range buckets here), in minutes.

## 6. What this run does not show (limitations)

- **Two pairs, one genome.** Both differences are dominated by one structural event: GigaDB drops the unplaced
  scaffolds, and the draft differs by 7,931 k-mers, all on its side. Neither pair is a population-scale variant
  set.
- **Who recovers what.**
  - Each KP-10 party recovers its OWN differing k-mers. The other side's differences appear only as roots of an
    unfactored polynomial, which is not implemented.
  - One km-peeler listing gives both sides' k-mers. The comparison counts one side's message for each method, and both
    parties' sketches are the same size.
  - The decode here ran both parties in one process, but each acceptance test used only one side's data (section 2).
- **Encodings and parameters are as shipped.**
  - Bytes are as the implementations store them: 8-byte entries for KP-10, 16-byte cells for km-peeler. Each could be
    packed tighter: 61-bit entries, a 62-bit key plus a few counter bits.
  - KP-10's m = 64 and the 0.8 m load rule were not tuned here.
  - km-peeler's r = 3 is its default and its smallest-threshold choice (c_3 ≈ 1.222).
- **Statistical reach.** 30/30 bounds failure below 9.5%. The 2^17 one-round KP-10 row and the D_est sweep are single
  fixed-seed runs. The IBLT protocol rows are 10 repetitions each.
- **The extended-syncmer method (their Section 3.5) was not rerun.** Its numbers (4.2 B per k-mer difference,
  1,986,807 spurious, 7 missed, 48,293,748 on draft vs final) and the round-1 critic's gap 5 (61% against the
  authors' 3%) stand as in `results/kp10_vs_shibuya_2026-11-10_prior.md:99-116`.
- **Shared host.** Timings are from a shared machine. The sweep phases ran up to 48 trials at once (task-clock 202,534
  CPU-s over 3,957 s, run_7), above the job's 16-slot grant on an uncapped host.
- **Driver outage.** The laila_perf driver lost its tunnel to laila's Galaxy before the job ended. Galaxy records the
  job as ok, last updated 22:49 UTC. Its outputs were fetched afterwards by `recon/kp10_vs_iblt_fetch.py` and are byte
  for byte the job's datasets.

## 7. Changes from the 2026-11-10 run (`results/kp10_vs_shibuya_2026-11-10_prior.md`)

| | 2026-11-10 | this run |
|---|---|---|
| IBLT scope | partitions 0–3 (6.3% of D), 30 trials | whole genome, 30 trials per ratio, 5 ratios, both pairs |
| IBLT at n = D, giga | 15/30 | 18/30 |
| IBLT at 1.05 D, draft | 26/30 (partitions 0–3, D = 488) | 30/30 (whole genome, D = 7,931) |
| KP-10 runs | one per field (numpy) | 30 OS-random seeds per pair, plus a D_est sweep (C) |
| KP-10 decode | both sets in one process | each side checks with its own data only |
| KP-10 input | canonical 2-bit codes | k-mers as written (S is the canonicalizer) |
| timing | numpy vs C, different hosts' loads | C vs C, one host, one thread |
| D underestimated | not tested | KP-10 and an IBLT resend protocol, D/8 to D |
