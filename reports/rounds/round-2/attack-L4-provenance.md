# Round 2, attack L4: provenance

Lens: does each published number and sentence come from the run, code and source it is attributed to?

Prior rounds searched with `scripts/round history` ("km-peeler seed", "separated", "independent", "syncmer", "20,000",
"numpy", "unplaced", "element by element", "x and lam"). Where a round-1 ruling bears on a finding it is named in that
finding. Method cards: none were declared in the brief's appendix. I found them in `../library/methods/`, and F3
cites `pseudo_replication.md`.

Run provenance established first. `run.log:3-7` of job b69c5c76 records sha256 of `recon/kp10c.c`, `recon/kmtools.c`,
`hash/kp7e.c`, `hash/kp7e.h` and `recon/kp10_vs_iblt_job.sh`. All five equal `git show e76fad1:<f> | sha256sum` and
HEAD. The 2bit md5s in `run.log:10-12` equal the files in `bench/data/hub-build/2bit`. (`hash/kp7seed.h`, which sets
how x and λ are drawn, is not in the hashed list. It is a small gap and does not change any number.)

**A steering note, as rule "Who is running this round" asks.** The claim list's C5 and C7 ask whether the 30 runs are
independent "(seeds, host, shared inputs)". They do not mention that km-peeler derives its hash functions from the
seed. That derivation is where the independence fails (F3).

## F1 -- "the GigaDB sketch ran before x and λ were separated" describes a superseded run

**Claim:** C11 (second sentence). **Verdict: REFUTED.**

- Paper, Implementation: "the GigaDB sketch ran before $x$ and $\lambda$ were separated."
- The x/λ separation is commit 8c34af9 (2026-11-11 01:14). The paper's GigaDB numbers come from kp10c at e76fad1
  (2026-11-11 13:47; source hashes above). `git merge-base --is-ancestor 8c34af9 e76fad1` succeeds.
- The run itself uses two seeds. `run.log:28` (timed_giga) and `:40-46` (sweep) record `"seed_x": "0x5eed",
  "seed_lam": "0x1a3d"`. The 30 seeded runs used `--seed-x os --seed-lam os` (`kp10_vs_iblt_job.sh:163`), and the
  x values in the 60 `kp10_seed` lines are all distinct.
- The sentence was true of the 2026-11-10 numpy run (`results/kp10_vs_shibuya_2026-11-10_prior.md`). That run predates
  8c34af9. It is now cited only through the ledger's "both fields in the numpy run" (ledger row for 11.3).
- It is probably a leftover of round 1's L1 F1 repair (NARROWED: "state that ... the experiments use fixed seeds").
  When the numbers moved to the C run, the sentence was not updated.
- The first sentence of C11 holds. kp10c hashes as written and keys buckets on `v >> (61 - bits)` (`kp10c.c:1-13`).
  Canonicalization happens only to score against the truth (`kp10c.c:27-28`).
- **Correction:** delete the clause. Or write: "The GigaDB and draft runs drew x and λ from separate seeds (fixed
  seeds 0x5eed and 0x1a3d for the sweep, OS randomness for the 30 runs per pair)."

## F2 -- the k = 63 run was not hashed by the numpy prototype, and the 20,000-k-mer check is about other code

**Claim:** C3. **Verdict: NARROWED.**

- Paper, Methods: "The $k = 63$ run used a numpy prototype, whose hashes match the C code on 20,000 sampled k-mers,
  with buckets sized from the true difference."
- **Who hashed.** Every element value at k = 63 came from the C streaming roller `hash/kp7e_stream.c`
  (`recon/aswritten.py:12,47,53`; `recon/run_aswritten_pair.sh:28`). The KMC truth was hashed by the C reference
  `kp7e_direct` + `kp7e_sym` (`aswritten.py`, `s_of_strings`). Numpy is used only for the sketch and the decode
  (`recon61`).
- **What the 20,000 check covers.** The check (ledger row "20,000" → `results/kp10_vs_shibuya_2026-11-10_prior.md:33`)
  is `recon/kp7e_hash.py:self_test`. It hashes random **2-bit packed uint64 codes at K = 31** in numpy and compares
  them with C. A 63-mer does not fit that code. So the check tests a numpy path the k = 63 run never used for its
  elements, and nothing checks the numpy code the run did use (the decoder) except the final comparison with KMC.
- **The bucket rule differs from the stated one.** The k = 63 run used `bits_for`, which requires at most **m/2**
  differences per bucket (`aswritten.py:136-142`; run JSON `"bits_rule": "auto: bits_for(exact value D, m)"`). The
  Methods state 0.8m. At D = 21,382, 0.8m gives 2^9 buckets (41.8 per bucket ≤ 51.2) and 274 KB. The run used 2^10
  and 548,864 B (`D/run_..._k63_sym_bauto_r12.json`, `bits: 10`). So "549 KB" and the derived 25.7 B/diff come from
  a sizing rule that is twice as conservative as the one in the paper. Whether 2^9 would finish in one round is
  untested.
- **"Matching KMC element by element"** is a match in S-value space. `verify_k63_kmc.json` compares the recovered
  values with S of KMC's 63-mers (`true_distinct_hash_values: 21382`). No recovered value was mapped back to a 63-mer.
- The numbers themselves reproduce: 21,382 exact, 548,864 B, 1 round, 0 missed and 0 spurious
  (`round1_owed_2026-11-11/verify_k63_kmc.json`, `kmc_k63_summary.json`).
- **Correction:** "The $k = 63$ run hashed every k-mer with the C roller and decoded with the numpy prototype, with
  $2^{10}$ buckets from the true difference (at most $m/2$ per bucket). Its recovered values equal $S$ of KMC's
  63-mer difference." Drop the 20,000-k-mer clause, and its ledger row, from this sentence.

## F3 -- km-peeler's 30 "independent random seeds" share hash functions, so the IBLT trials are not independent

**Claims:** C7 and C5 (also C6, C1 for the IBLT's 30 of 30). **Verdict: C7 REFUTED; C5 NARROWED.**

- Paper: "each method ran 30 times with independent random seeds on one host"; and "Each 30 of 30 bounds the failure
  rate only below 9.5 percent".
- **The seeds are consecutive integers, not random.** `iblt_sweep` passes `-s $((1000 + t))`, t = 1..30
  (`recon/kp10_vs_iblt_job.sh:176`). `run.log` shows `"seed": 1001` ... `1030`.
- **km-peeler derives its r = 3 hash seeds as root_seed + i.** See
  `bench/data/kmer_poly_work/km-peeler/ibflib.c:27-29` at the pinned b8ec39b: `j = root_seed + i; /*FIXME use a better
  method?*/`. The table is partitioned into r chunks of equal size (`ibflib.c:121-122`).
  - Trial t uses Murmur seeds {t, t+1, t+2}, and trial t+1 uses {t+1, t+2, t+3}.
  - So adjacent trials share two of their three hash functions, and trials two apart share one.
  - The 30 trials draw on only 32 distinct functions.
- **The dependence shows in the data.** These are the success strings at n = D, in seed order (`run.log`
  `iblt_trial` lines; my tally):
  - draft: `100000000000111111000001000000` (8 of 30). Six consecutive seeds, 1013 to 1018, all succeed. A runs test
    gives 6 runs against 12.7 expected (z = −3.2). A permutation test (200,000 shuffles) gives p = 0.0019.
  - giga: `011101011101011111000100011011` (18 of 30). No clustering (p = 0.66).
  - The rows the paper quotes (1.05·D, 30 of 30) cannot show clustering. But they are built from the same overlapping
    function windows.
- **The bound.** Clopper–Pearson's 9.5% assumes 30 independent Bernoulli trials. The project's own card says "Your
  effective sample size is the number of independent decisions, not the number of rows"
  (`../library/methods/pseudo_replication.md:8`). About 10 of the 30 trials have disjoint function triples (seeds
  1001, 1004, ...). With 10 trials, 10/10 bounds failure only below 25.9%.
- **KP-10's side is clean.** It used 30 OS-random (x, λ) draws (distinct x in all 60 `kp10_seed` lines), so 9.5% is
  valid for the sketch.
- **The same pattern elsewhere.** The extended-syncmer run used seeds 2000 to 2004 (`recon/kmpeeler_run.py:167`). The
  protocol sweep's seeds (100000k + 1000r + 1000f, `job.sh:191`) are spaced by at least 125 and do not overlap.
- No prior round raised this (`scripts/round history "km-peeler seed"`: searched absence).
- **Correction:** "km-peeler ran 30 trials with consecutive root seeds, and because it derives its three hash seeds by
  adding 0, 1 and 2, adjacent trials share two hash functions. The 9.5 percent bound applies to the sketch, whose 30
  runs drew independent OS-random parameters." Or rerun `iblt_sweep` with seeds spaced at least 3 apart, or random.
  That is a laila job, owed rather than run here.

## F4 -- the extended-syncmer row comes from a different run, never listed completely, and its 7 misses are deterministic

**Claim:** C8. **Verdict: NARROWED** (on "2.7 times smaller" and on round 1's peeling explanation).

- `results/kp10_vs_shibuya_2026-11-11/full_kmp_ext_giga.json`, the 2026-11-10 run by `recon/kmpeeler_run.py`, has 5
  trials. All 5 give **identical** outcomes: `listed 507219` of `D_ext 507221`, `complete: False`, `missing 7`,
  `spurious 1986807`, `bytes 13667246`.
- The trials used root seeds 2000 to 2004, and by F3 seeds 2000 and 2004 share no hash function. An IBLT peeling
  failure depends on the hash functions, so five identical shortfalls of exactly 2 syncmers point to something
  deterministic (input, encoding or harness), not to random peeling.
- **Round 1's ruling.** Round 1's completeness adjudication (Gap 5, NARROWED) wrote "The misses are explained by
  peeling, not by a harness defect". Its evidence was "the listing returned 507,219 of 507,221 every trial". The
  seed derivation was not known then. With it, that same evidence argues against random peeling. The ledger's
  carried-items row already says "harness error not ruled out" (OPEN). This finding gives a concrete reason to
  suspect one.
- **Not like for like with the other rows.**
  - The listing never completed in any trial. Under the resend rule applied to the k-mer table (resend when the
    listing fails), this configuration would need another round. Yet Table 4 shows it as "1" round, and "2.7 times
    smaller" (11.27 / 4.2) compares the sketch's complete recovery with a listing that was incomplete.
  - n was set from the true extended-syncmer difference (`D_ext`, `kmpeeler_run.py:156-164`).
  - It ran 5 trials, in the superseded harness, not in job b69c5c76. Yet the paragraph above Table 4 says "each method
    ran 30 times ... on one host", and the results file says this method "was not rerun" (§6, :209-211).
- **Correction.**
  - Table 4's row should read "no: listing incomplete in 5 of 5 trials, 1,986,807 spurious and 7 missed", with a
    footnote that it comes from the earlier 5-trial run.
  - Drop "2.7 times smaller", or qualify it as "at a size at which its listing did not complete".
  - Before the 61%-against-3% comparison is used against the competitor, identify the 2 unlisted extended syncmers
    (a cheap local check on the saved `kmp_ext_*` files).

## F5 -- the sketch's size parameters were selected on the test pair, like the table's

**Claims:** C7 and C1. **Verdict: C7 NARROWED; C1 NARROWED.**

- Paper: "The k-mer table was sized to succeed in 30 of 30 trials". The sketch, by contrast, is described by a rule:
  "the smallest power of two of buckets that holds at most 0.8m estimated differences".
- **Where m and the bucket count came from.** `results/kp10_vs_shibuya_2026-11-10_prior.md:47-54` is a "Partition 0
  calibration, primary pair". It tried four (buckets, m) settings on a slice of the same GigaDB pair: 2^19/16,
  2^18/24, 2^17/40 and 2^16/64. It kept the one with the lowest B/diff, 2^16 with m = 64, at **11.3**.
- **How the 0.8m rule was set.** `recon/kp10c.c:31-32` says: "this rule reproduces the 2^16 and 2^8 buckets of the
  2026-11-10 runs at the true D; it is fixed before any run here".
  - Any threshold t with 49.4 ≤ t < 62.0 differences per bucket gives both 2^16 (3,238,101 / 2^16 = 49.4) and 2^8
    (7,931 / 2^8 = 31.0, and 2^7 would give 62.0).
  - 0.8m = 51.2 sits in that window. The rule was fitted to the two pairs it is then tested on.
- **"Not tuned here".** The results file says "m = 64 and the 0.8 m load rule were not tuned here" (:205). That is
  true only of job b69c5c76.
- **Not symmetric.** Both methods' sizes were chosen with outcomes from the test pairs. The table was chosen for
  success, and the sketch for minimum bytes on a 1.6% slice of the GigaDB difference. The paper discloses only the
  table's selection.
- **Both grids favour exact multiples.** The D/8 to D estimates are exact powers of two of D, so the table's resend
  doubling lands on exactly the 1.05·D it was tuned to succeed at (`iblt_protocol.tsv`: 4 rounds in all 10 reps, with
  identical bytes). The sketch's power-of-two buckets land exactly as well. This favours the table, so it does not
  reverse the 20.8-against-38.5 comparison. It does mean that neither curve is a draw over arbitrary estimates.
- **Correction:** add to Methods: "We chose m = 64 and $2^{16}$ buckets in a calibration on one 64th of the GigaDB
  pair, and set the 0.8m load so that it reproduces those bucket counts. Both methods' sizes were thus chosen on these
  pairs." A held-out pair would be needed to call either rule a priori.

## F6 -- "66 field multiplications per element" is cited from the superseded prototype and omits the decode

**Claim:** C12. **Verdict: NARROWED.**

- Ledger row "m + c = 66 | field multiplications per element" cites `results/kp10_vs_shibuya_2026-11-10_prior.md:8`,
  the numpy prototype's verdict. It does not cite the C code behind the paper's numbers.
- In kp10c each side also **evaluates its degree-d polynomial at every element it holds in every differing bucket**
  (`kp10c.c:226-228` `flag_roots`, called at `:294-295`).
  - On giga, 3.2 million differences fall in 65,536 buckets (about 49 per bucket), so nearly every bucket differs and
    nearly every element pays about d ≈ 49 more multiplications.
  - Elements in the 1,263 failed buckets pay again in round 2.
- The measured split agrees: 49.5 CPU-s of power sums against **110.3** CPU-s of decode, versus 3.5 s for the IBLT's
  diff and list (`results/kp10_vs_shibuya_2026-11-11.md:173-174`).
- So "about m + c = 66" understates the per-element work at this difference by about 2×. The paper omits the times
  that would show it.
- **Correction:** "about m + c = 66 multiplications per element to build the sketch, and about as many again to find
  the roots when most buckets differ. On the GigaDB pair the decode took 110 CPU-seconds against 3.5 for the table's
  listing." Cite `kp10_vs_shibuya_2026-11-11.md:173-174`.

## F7 -- ledger rows that do not say what the paper says

**Claim:** C15. **Verdict: REFUTED.**

1. Row "20,000 | sampled k-mers on which the numpy hashes equal the C | prior.md:33". The source is a k = 31,
   packed-code check from the 2026-11-10 run. The paper attaches it to the k = 63 run (F2).
2. Row "m + c = 66 | prior.md:8". This is the superseded source, and it omits the decode (F6).
3. The Implementation sentence "the GigaDB sketch ran before x and λ were separated" has no row, and is false for the
   run the rows cite (F1).
4. Row "0.8m, D/8 to 4D | sizing rule | section 2". Section 2's own text (and `kp10c.c:31-32`) says the rule was set
   to reproduce earlier bucket counts. The row presents it as a design rule (F5).
5. "about 3 percent in the authors' tests" has no number row pointing to Shibuya et al. It appears only as quoted text
   in the carried-items table. The value does exist in the held summary
   (`../library/papers_md/summaries/shibuya_2022_genomic-reconciliation.md:27-29`, the paper's lines 782-783: "within
   about 3% of the true set difference size"). The ledger should cite it.
6. Row "2.4× | 24.67 / 10.24 ... across the two programs" (`paper_bench_2026-11-12.md:33, 51`). The row discloses
   that the ratio crosses programs. Its own source forbids that: "compares against the baselines only within
   bench_kept" (`paper_bench_2026-11-12.md:77`), and so does the Figure 3 caption ("The panels are not compared with
   each other"). See F8.
7. Row "9.5% ... both KP-10 (30 OS seeds) and the IBLT (30 seeds)". The IBLT half rests on non-independent trials
   (F3).

**Correction:** fix rows 1, 2, 4 and 7, and add rows for 3 (or delete the sentence) and 5.

## F8 -- the Discussion's 2.4× compares across two harnesses, which the paper's own caption rules out

**Claim:** C14. **Verdict: NARROWED.**

- "2.4 times as fast as our hash" is 24.67 (bench_e, Eisenstein every k-mer, 2 lanes) / 10.24 (bench_kept, universal
  hash every k-mer). The ledger row says "across the two programs".
- Figure 3's caption refuses that comparison because the two programs time the same sampler code at 28.7 and 22.5
  cycles per kept k-mer. `paper_bench_2026-11-12.md:70-77` says baselines are compared "only within bench_kept".
- bench_kept has no Eisenstein every-k-mer row.
- The size of the error is probably small. The two harnesses' ntHash rows agree to 0.5% (11.37 and 11.43), and the
  ratio via ntHash is 2.17 × 1.12 ≈ 2.4. But the number's provenance breaks the paper's stated rule.
- The 1.12× (11.43 / 10.24, both in bench_kept) is clean.
- **Correction:** either time Eisenstein every-k-mer inside bench_kept, or write "about 2.2 times (through each
  program's own ntHash)" and say so.

## F9 -- "almost all from unplaced scaffolds" is asserted, not counted

**Claim:** C10. **Verdict: NARROWED** (the counts and KMC agreement survive).

- **The counts.** 3,238,101 = 3,237,326 + 775, and KMC agrees element for element (`tables/truth.tsv:2-3`). KMC ran
  as an independent tool in the same job (`kp10_vs_iblt_job.sh:90-95`).
- **The scaffold attribution was never counted.** `results/kp10_vs_shibuya_2026-11-11.md:42` and `prior.md:27`
  assert it from record counts (147 against 10, `tables/sets.tsv`). No file assigns the A∖B k-mers to records.
- **My spot check is consistent with the claim but does not prove it.**
  - twoBitInfo gives the 10 GigaDB chromosomes exactly the lengths of CM028014.1 to CM028023.1, and Cannbio-2's
    other 137 records total 10,317,051 bases.
  - One pair, CM028022.1 against Cs_Cb2_09 (71,141,854 bp), is byte-identical after uppercasing (md5 e2c1697c... for
    both).
  - B∖A = 775 shows that not every chromosome is identical, so the release does more than drop scaffolds.
- **Correction:** count A∖B k-mers by Cannbio-2 record (one pass of the rescan kp10c already does). Or write "most of
  it on the 137 unplaced scaffolds (10.3 Mb) that the release drops" with that count.

## Other claims

- **C2 -- SURVIVED.**
  - 20.76 and 38.49 are at `tables/kp10_dest.tsv:2` (run.log:40) and `iblt_protocol.tsv:6`.
  - The protocol resends only on `list_exit` (`job.sh:186-199`) and never reads the truth.
  - The IBLT's protocol bytes are deterministic across its 10 reps, so its mean is not a smoothed value.
  - F5's grid caveat applies, and it favours the table.
- **C4 -- SURVIVED.**
  - 21.70 is at `kp10_dest.tsv:7` (one fixed-seed run, as stated). "About 1,250" is the 1,254 mean of
    `kp10_seeds.tsv:2`.
  - The round-2 bytes check out. 36,481,232 − 35,127,296 = 1,353,936 = 2 × 1,263 × 536, so each failed bucket was
    sent again as two halves, consistent with "split and resent".
- **C6 -- SURVIVED** on provenance (`kp10_seeds.tsv:3` 17.30 at 30/30; `iblt_trials.tsv` draft 1.05 20.54). F3 and F5
  qualify the IBLT's 30/30 and the 2^8 choice.
- **C9 -- SURVIVED.** `make_figures.py:249-258` reads `kp10_dest.tsv` and `iblt_protocol.tsv` directly. It plots only
  exact rows and asserts constant IBLT rounds. The caption matches the tables (10/10 `exact_cancelled` in every
  protocol row).
- **C13 -- UNCHALLENGED** under this lens. It concerns the literature and inference.

| claim | verdict |
|---|---|
| C1 | NARROWED |
| C2 | SURVIVED |
| C3 | NARROWED |
| C4 | SURVIVED |
| C5 | NARROWED |
| C6 | SURVIVED |
| C7 | REFUTED |
| C8 | NARROWED |
| C9 | SURVIVED |
| C10 | NARROWED |
| C11 | REFUTED |
| C12 | NARROWED |
| C13 | UNCHALLENGED |
| C14 | NARROWED |
| C15 | REFUTED |
