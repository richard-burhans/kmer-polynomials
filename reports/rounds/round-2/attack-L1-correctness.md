# Round 2 attack, lens L1: correctness

Attacker: correctness lens. Frozen text: `reports/rounds/round-2/ATTACKED-paper.md` (AP) and `ATTACKED-ledger.md` (AL).
Raw outputs: `results/kp10_vs_shibuya_2026-11-11/laila_b69c5c76/` (`L/`). km-peeler source read at the pinned commit
b8ec39b in `bench/data/kmer_poly_work/km-peeler` (the clone `recon/kp10_vs_iblt_launch.sh:16` bundles).

Searched before filing: `scripts/round history` for "one side's message", "own differing", "both sides", "own data",
"oracle", "unplaced", "km-peeler", "syncmer", "D/8", "canonical". Only one prior ruling touches these findings: round 1's
L5 F5 (adjudication-L5.md:95-112). It is disclosed under F1, the only finding it bears on. Nothing about D/8, seed
overlap in km-peeler, the x/λ sentence, or the k = 63 provenance was ever filed.

Method cards: the brief's appendix declares none. Rule 3 governs C5 only, and F4 uses no card. F4's premise is the
independence that a one-sided Clopper–Pearson bound assumes.

Local checks were all tiny and run under `nice -n 19` (twoBitInfo on two 2bit headers; Poisson tails; a
permutation test on 60 booleans). No genome was hashed.

## F1 -- "recovered every difference ... with 11.3 bytes" prices one message, but the scored recovery used both parties' sketches (C1)

**Claim.** AP:14 and AP:130: the sketch "recovered every difference in 30 of 30 runs with 11.3 bytes per difference
over two rounds, against 20.5 for a one-round" table.

**What the code does.** `recon/kp10c.c:257-307` (`decode`) forms `S = SA - SB`, which needs BOTH sketches. It flags A's
own roots of P (`flag_roots(pr, ..., ownA, ...)`, l.294) and B's own roots of Q (l.295). It then accepts a bucket only
if **both** parties' flag counts and **both** c-sum checks pass (l.296, l.301-304). `reconcile` (l.319-399) scores
`oa ∪ ob` against the truth. The byte count is `R.entries * 8` (l.555), which is one side's entries
(l.330: `R.entries += np * (M + 1)`, once per round).

On the GigaDB pair, A∖B = 3,237,326 and B∖A = 775 (`L/tables/truth.tsv:2`). A recovers its 3,237,326 only from B's
sketch, and B recovers its 775 only from A's. The "missing 0, spurious 0" of all 30 runs (`L/tables/kp10_seeds.tsv:2`;
e.g. `run.log:28`, `kmers_recovered_A_minus_B 3237326`, `kmers_recovered_B_minus_A 775`) is therefore the outcome of an
exchange in which **both** 36.5 MB sketches crossed. That is 2 × 36,481,232 / 3,238,101 = **22.5 B/diff**, not 11.3.
One message buys less. B→A gives A its own 3,237,326 k-mers. The 775 reach A only as roots of an unfactored Q, and only
as S values, never as k-mers, because S is not invertible (kp10c.c:24-25 says so). The table's 20.5 is one message,
from which one party lists **both** sides' k-mers (`L/tables/iblt_trials.tsv:9`; `kmtools listscore`,
`listed_A_minus_B` and `listed_B_minus_A`).

**What changed since round 1.** Round 1's L5 F5 adjudication ruled "the sketch's byte count does not depend on who
decodes" (adjudication-L5.md:103). That was correct for the one-process numpy decoder it judged. The C decoder now makes
each party recover its own k-mers, and "every difference" sums the two parties' recoveries. "Who decodes" now
determines which message the scored result needed. The paper body concedes the different outputs (AP:130, "Each party
of the sketch recovers only its own differing k-mers"). It never says that the 30/30 "every difference" needed both
messages, and the abstract carries neither caveat.

**Correction.** An optimised return message would close most of this, but it was not built or run. For example, A
learns deg Q per bucket and could send back only the ≤ 775 buckets with eB > 0. So the fix is wording, not a new
number. Abstract: "...each party recovered its own differing k-mers in 30 of 30 runs from the other's sketch of 11.3
bytes per difference over two rounds, against 20.5 for one table from which one party lists both sides'." Or price the
measured exchange: "22.5 bytes per difference for the two sketches exchanged". **NARROWED**: the numbers reproduce, and
the like-for-like test the claim sets itself ("one side's message") fails for "every difference".

## F2 -- at D/8 every bucket overflowed: the sketch resent everything at double resolution, just like the table, and both degraded by the same factor (C2, C9, C13)

**Claim.** AP:14: "Because an overflowing bucket is resent alone, an eightfold underestimate of the difference raised
this to 20.8, against 38.5 for the table." AP:130 repeats the causal "because".

**Arithmetic from `L/tables/kp10_dest.tsv`.** The giga rows are 2-4 and 6. At D/8, bits0 = 13 (8,192 buckets) for
D = 3,238,101, so about 395 differences per bucket against m = 64. At D/4 the load is about 198. Every bucket overflows
in rounds 1 and 2. The split counts show the protocol replaying the D/2 run exactly:

| D_est | splits | = |
|---|---|---|
| D/2 | 34,026 | |
| D/4 | 50,410 | 16,384 + 34,026 |
| D/8 | 58,602 | 8,192 + 16,384 + 34,026 |

The bytes agree to the byte: 536 × 8,192 + 536 × 16,384 + 54,039,520 = **67,212,256**, which is the D/8 row. So the
D/8 sketch is two whole-sketch transmissions, each at double the last resolution, followed by the D/2 run. D/2 itself
is one near-whole resend (32,768 buckets at a load of about 99) followed by the true-D run: 536 × 32,768 + 36,481,232 =
54,044,880, within 10 buckets of 54,039,520. That is the same geometric doubling as the table's resend-at-double
protocol (`recon/kp10_vs_iblt_job.sh:187-200`).

**The factors are the same.** The cost relative to each method's own true-D run is:

| D_est | sketch | table |
|---|---|---|
| D/2 | 16.69 / 11.27 = 1.48 | 30.79 / 20.53 = 1.50 |
| D/4 | 1.72 | 1.75 |
| D/8 | 1.84 | 1.87 |

The sketch-to-table ratio at D/8 (38.49 / 20.76 = 1.85) is the ratio at the true D (20.53 / 11.27 = 1.82). On the
abstract's pair, the underestimate taxed both methods alike. The 20.8 vs 38.5 gap is inherited whole from the true-D
gap, which F3 shows is partly one-round efficiency, not resending alone. Selective resend did act on the draft pair
(21.36 / 17.30 = 1.23 against 1.88), but there it rides on power-of-two granularity: the true-D draft sketch sat at a
load of 31.

**Correction.** "An eightfold underestimate made every bucket overflow, so each round resent the whole sketch at double
resolution; the cost rose 1.84-fold to 20.8, against 1.87-fold to 38.5 for the table resent at double size." Drop
"Because an overflowing bucket is resent alone" from the abstract and from AP:130. Drop "degrades more gently" for
giga. Figure 4A's "split and resend this bucket only" may stay as the mechanism, but the 4B caption should not imply
that it produced the gap. **REFUTED** (the causal clause). The numbers themselves reproduce.

## F3 -- "sized for a single round, the sketch needed 21.7" is the next power of two, not a need; and "resent about 1,250" is half the buckets sent (C4, C13)

**(a) 21.7 is a granularity artifact.** `kp10c.c:30-31` permits only 2^b buckets. The only one-round configuration run
was 2^17 (`kp10_dest.tsv:7`), which is exactly twice round 1's 10.85 (`kp10_seeds.tsv:2`, `bytes_round1`). The IBLT
was given 5% steps (0.95, 1.0, 1.05, 1.1, 1.25 D; `kp10_vs_iblt_job.sh:31`). A Poisson-load model reproduces the
observation:

- **Calibration.** It predicts **1,253** overfull buckets at 2^16 and m = 64. The observed mean is **1,254** over 30 OS
  seeds (`kp10_seeds.tsv:2`, `mean_splits`).
- **With m raised and 2^16 buckets kept** (m is a CLI parameter; `MAXM` = 256), the expected overfull buckets in one
  round are:

  | m | expected overfull | B/diff |
  |---|---|---|
  | 80 | 1.5 | 13.4 |
  | 88 | 0.017 | 14.7 |
  | 96 | 1e-4 | 16.0 |

- **With m = 64 and about 100,000 value-range buckets** (if ranges need not be powers of two): 0.03 expected overfull,
  16.6 B/diff.

So a one-round sketch is projected at about 15-16.5 B/diff, below the table's 20.5. Of the 9.2 B/diff gap, roughly 4-6
is one-round efficiency and 3-5 is the second round. "Its smaller size comes from the second round" overstates that.
This is a calibrated projection, not a measurement. **Owed**, cheap, laila only: one `kp10c --m 88 --dest D` run on
giga.

**(b) Buckets resent.** AP:128 says "resent about 1,250 of them in its second round". AL:178 says "about 1,250 |
buckets resent in round 2". But 1,254 is the mean number of failed buckets. Each is split and **both halves** are sent
(`kp10c.c:388-395`). In the fixed-seed run, round 2 is 36,481,232 − 35,127,296 = 1,353,936 B = **2,526 buckets** of
536 B (`run.log:28`). Write "split about 1,250 overfull buckets and sent their halves (about 2,500)".

**Correction.** Change "Sized for a single round, the sketch needed 21.7, so its smaller size comes from the second
round" to "With twice the buckets, enough for one round, it used 21.7; finer sizing was not tried." **NARROWED.**

## F4 -- km-peeler's 30 trials share hash functions: seeds are consecutive and each trial's three hashes are seed, seed+1, seed+2 (C5, C7)

**Code.** `km-peeler/ibflib.c:23-31`, `ibf_hash_init`, has `j = root_seed + i; /*FIXME use a better method?*/`. Block i
places a key at `hash(seed_i) mod chunk_size + i*chunk_size` (ibflib.c:331). The harness gives trial t the root seed
1000 + t (`recon/kp10_vs_iblt_job.sh:176`). So trial 1001 hashes with Murmur seeds {1001, 1002, 1003}, and trial 1002
with {1002, 1003, 1004}. At a fixed n, consecutive trials share **two of three** key-to-cell maps (up to block order,
which peeling ignores), and trials two apart share one. The 30 trials at one ratio draw on 32 hash functions, not 90.
The seeds are fixed integers, not random.

**The dependence shows in the data.** Success by seed order from `run.log` `iblt_trial` lines:

- draft at n = D: `100000000000111111000001000000`. A permutation test on that 8-of-30 pattern gives p = 0.0011 for a
  run of six successes and p = 0.0065 for only three runs.
- giga at n = D: no clustering detected (p ≈ 0.67).

**Consequences.**
- C7's "each method ran 30 times with independent random seeds" (AP:119) is false for the table. KP-10's 30 runs did
  use OS-random, independent x and λ (`run.log:48`, `"seed_x": "os"`).
- C5's 9.5% (AP:130) is computed correctly (1 − 0.05^(1/30) = 0.0950). But it assumes 30 independent trials, which
  the table's trials are not, so for the table it is not a valid 95% bound.
- 1.05 was chosen as the smallest grid ratio that went 30/30 on these same 30 seeds (Table 4 caption, "sized to
  succeed in 30 of 30"). The bound is then also post-selection.
- "Neither method returned a wrong difference" stands (`iblt_trials.tsv`, `undetected_failures` 0;
  `kp10_seeds.tsv`, `spurious_total` 0).

**Correction.** "...the sketch 30 times with independent random seeds, and the table 30 times with root seeds 1001 to
1030, which km-peeler expands to overlapping hash functions." Restrict the 9.5% to the sketch, or rerun the table with
root seeds spaced at least r apart (e.g. 1000·t), which is cheap. **NARROWED** for C5 and C7.

## F5 -- the Implementation's "the GigaDB sketch ran before x and λ were separated" and "our experiments used fixed seeds" describe the superseded numpy run, not job b69c5c76 (C11)

AP:51 says: "Our experiments used fixed seeds for reproducibility, so they test the construction rather than the bound,
and the GigaDB sketch ran before $x$ and $\lambda$ were separated." The paper's GigaDB numbers come from the C run:

- The fixed-seed runs log `"seed_x": "0x5eed", "seed_lam": "0x1a3d"`, giving x = 492247986318507 and
  λ = 834835136232447260 (`run.log:28, 44`). These are drawn from the separate, domain-tagged streams introduced in
  round 1 (`results/round1_owed_2026-11-11.md:15-26`). The same x and λ appear in the post-separation k = 63 run
  (`results/round1_owed_2026-11-11/run_..._k63_sym_bauto_r12.json`).
- The 30/30 headline runs used OS randomness for both x and λ (`run.log:48`, `"seed_x": "os", "seed_lam": "os"`). So
  "our experiments used fixed seeds" is false for the result the abstract quotes. Those runs do meet the bound's
  premise, contrary to "test the construction rather than the bound".

The first sentence of C11 holds for the C run. It hashes as written (`kp10c.c:149-160`, `kp7e_sym_seq`) with buckets
as top-bit ranges of S (`slice`, l.171-176). **REFUTED** (the second sentence). Replace it with: "The 30 runs per pair
drew x and λ independently from the operating system. The size sweep used one fixed pair of seeds."

## F6 -- the k = 63 run did not hash in numpy, and the 20,000-k-mer check cannot cover it (C3)

AP:47 says: "The $k = 63$ run used a numpy prototype, whose hashes match the C code on 20,000 sampled k-mers."

- **Who hashed.** The k = 63 values were hashed by the C streaming hasher `hash/kp7e_stream.c`
  (`recon/aswritten.py:12`, `run_aswritten_pair.sh:28`), and the truth was hashed by C `kp7e_direct`
  (`aswritten.py:226-233`). numpy only reconciled.
- **What the 20,000 check covers.** It is `recon/kp7e_hash.py:85-98`, `self_test`. It runs at the module's K (31) on
  random **2-bit packed codes**, which cannot exist at k = 63 (packing is "exact for k ≤ 32", AP:55). The ledger's
  source for it is the superseded run (AL:115, `kp10_vs_shibuya_2026-11-10_prior.md:33`), under the old seed
  derivation.

So the sentence cites a check of a different program, at a different k, under different seeds, for a step numpy did not
perform. Also, "matching KMC element by element" was checked in hash space: S of KMC's 63-mers against recovered S
values (`verify_k63_kmc.json`), not k-mers. That is adequate under the bound, but should be said.

**Correction.** "The $k = 63$ run hashed with the same C code and reconciled in a numpy prototype, with buckets sized
from the true difference. Its recovered values equal the hashes of KMC's 63-mer difference." The core of C3 ("linear",
"needs no canonicalization": the C run hashes as written, and canonical codes appear only in scoring, `kp10c.c:27-28,
440`) **survives** my check. **NARROWED** on the k = 63 provenance.

## F7 -- the extended-syncmer table never listed completely, identically across five seeds; "which we did not trace" is contradicted by its own JSON (C8)

`results/kp10_vs_shibuya_2026-11-11/full_kmp_ext_giga.json` shows:

- **Trials.** All five trials (seeds 2000-2004, `recon/kmpeeler_run.py:167`) report `"listed": 507219`,
  `"complete": false`, `"missing": 7`, `"spurious": 1986807` against `"D_ext": 507221`.
- **Sizing.** n = 1.05 × D_ext.
- **What `complete` means.** It is `listed == D_ext` (kmpeeler_run.py:179). Every trial's listing was **incomplete**: 2
  differing extended syncmers were never listed.

A peeling failure depends on the seed. Five different seeds giving byte-identical outcomes points to a deterministic
cause in the input or harness, not a sizing property. Candidates include strings km-peeler cannot represent at its
fixed 51-character payload, and lines `comm -3` counts but `ibltseq` drops. That is exactly round 1's open "harness
error not ruled out" (critic gap 5). The 7 missed k-mers are traced to the 2 unlisted syncmers by the run's own record.

**Consequences.**
- Table 4's "4.2 bytes per difference" and AP:134's "2.7 times smaller" are the size of a table that failed to list in
  every trial, not of a working one.
- This row is 5 trials from the 2026-11-10 numpy-era harness, set beside 30-trial rows from job b69c5c76. The paper does
  not say so (the results file does, §6).

**Correction.** Table 4: "4.2, listing incomplete in 5 of 5 trials (2 of 507,221 syncmers unlisted, identically for
every seed)". AP:134: drop "which we did not trace" and state the incomplete listing. Keep the row only with "a
probable harness defect, not investigated". **REFUTED** (as a measured size of a working method).

## F8 -- "m + c = 66 multiplications per element" counts the sketch build, not recovery, and against the wrong unit (C12)

AP:136 says: "Recovery costs more computation than a table listing. The sketch needs about $m + c = 66$ field
multiplications per element, against 3 cell updates."

- **The build.** `psums` (kp10c.c:178-198) does 66 multiplications per element. That is building the sketch.
- **Recovery.** `flag_roots` (l.228-237) evaluates P, of degree eA ≤ m, at **every own element of every bucket with a
  difference**. On giga that is all 65,536 buckets (`run.log:28`, `buckets_with_difference_round1: 65536`) at about 49
  differences per bucket. Side A therefore pays about 66 + 49 ≈ 115 multiplications per element, plus round 2's
  repeats.
- **Measured.** Decode took 110.3 CPU-s against 49.5 for A's sketch (`L/tables/timed.tsv:2`).
- **The table's unit.** Each of the table's "3 cell updates" first costs a MurmurHash3_x64_128 of the 31-byte key
  (`ibflib.c:67`), which is several 64-bit multiplications. "66 field multiplications against 3 cell updates" compares
  unlike units.

**Correction.** "Building the sketch takes $m + c = 66$ field multiplications per element; recovery evaluates a
polynomial of degree up to $m$ at every element of each bucket that differs, which on the GigaDB pair was every bucket
(110 against 3.5 CPU-seconds for listing the table)." **NARROWED.**

## F9 -- ledger rows that do not say what their source says (C15)

- **AL:115** "20,000 | sampled k-mers on which the numpy hashes equal the C". The source line is right, but the paper
  attaches it to the k = 63 run, which it cannot cover (F6).
- **AL:150** "1,275 | GigaDB buckets resent in round 2" and **AL:178** "about 1,250 | buckets resent". These are split
  (failed) buckets. Buckets resent are twice that: 2,526 in the fixed-seed C run (F3b). The cited line
  (`kp10_vs_shibuya_2026-11-11.md:118`) itself says "split in two and resent".
- **AL:109** "4.2 | extended-syncmer IBLT bytes per true k-mer difference". It omits that the listing was incomplete in
  every trial (F7). `prior.md:102-105` gives "507,219 of 507,221", and the raw JSON gives `complete: false`.
- **AL:183** "3 hash functions | km-peeler default r = 3". This is true, but it is the row a reader would check for
  seed independence, and it omits root+i (F4).
- **AL:159** "1 of 37 | KP-10 runs with an S collision". The 7 fixed-seed runs share one (x, λ) (`run.log:28, 44-46`),
  so the distinct draws are 31, not 37. This is ledger-only and not in the paper.

Every other row checked against its cited line matched: AL:100-105, 152-156, 174-184 and the Figure 4B row 214. Every
number in AP:117-136 has a row. **NARROWED.**

## Note on the claim list (rule "who is running this round")

The list frames C2 and C13 around "resent alone" and tests whether the table's protocol is a straw man. It does not ask
the question F2 answers: whether the sketch at D/8 was doing anything other than resending everything.

## Verdicts

| claim | verdict | finding |
|---|---|---|
| C1 | NARROWED | F1 (two sketches exchanged for "every difference"; 22.5 B/diff for the scored exchange) |
| C2 | REFUTED | F2 (causal "because" refuted; both degrade 1.84× vs 1.87×; the D/8 run is whole resends then the D/2 run) |
| C3 | NARROWED | F6 (k = 63 provenance sentence wrong); linearity and no-canonicalization survive my check |
| C4 | NARROWED | F3 (21.7 is a power-of-two artifact, projected one-round ~15-16.5; "resent 1,250" is 2,526 halves) |
| C5 | NARROWED | F4 (bound arithmetic right; IBLT trials not independent; post-selected ratio) |
| C6 | SURVIVED | 137,216 / 7,931 = 17.30 and 162,938 / 7,931 = 20.54 (`kp10_dest.tsv:11`, `iblt_trials.tsv:3`); A∖B = 0 on draft, so one message A→B delivers all 7,931 and the one-message count is like for like here |
| C7 | NARROWED | F4 (table seeds consecutive, hash functions shared; "independent random" false for the table) |
| C8 | REFUTED | F7 (incomplete listing in 5/5 trials, identical across seeds; "did not trace" contradicted) |
| C9 | SURVIVED | Figure 4B values match `kp10_dest.tsv:2-4,6` and `iblt_protocol.tsv:6-9` (`make_figures.py:249-258` asserts exact and constant rounds); the caption is accurate, but the causal title inherits F2 |
| C10 | SURVIVED | counts reproduce (`truth.tsv:2-3`, KMC identical); twoBitInfo: the 10 chromosomes have identical lengths in both releases, the 137 dropped records total 10,317,051 bp, and B∖A = 775 bounds chromosome edits |
| C11 | REFUTED | F5 (the x/λ sentence and "fixed seeds" describe the superseded run) |
| C12 | NARROWED | F8 (66 is the build; recovery is about 115 per element on giga side A; unit mismatch) |
| C13 | NARROWED | F2, F3 (the advantage is not "due to resending only overfull buckets" at D/8, and partly one-round efficiency at D); "must be resent whole" is a property of km-peeler's CLI, not of IBLTs, whose cells are also linear |
| C14 | SURVIVED | "with a proof" checked: `hash/kp7kept.h:31-33`, multiply-shift with a random odd 128-bit a, collision ≤ 2^-63 (Dietzfelbinger 1997); timings not re-examined (another lens) |
| C15 | NARROWED | F9 |
