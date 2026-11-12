# Round 2 adjudication of L4 (provenance)

Adjudicator: an independent skeptic. The target is the attack. I read the attack first and the claims and their evidence
second, as the brief orders.

**Seal.** `sha256sum attack-L4-provenance.md` starts `4b1d1ee5555b`, and `round.json` records the same full digest. The
file I ruled on is the sealed one.

**What I re-ran.** Everything here is small. It reads committed raw outputs, the pinned km-peeler source, or one saved
0.3-second listing, run under `nice -n 19 ionice -c3`. Nothing was hashed or counted over a whole genome.

- I re-hashed the five sources in `run.log:3-7` at e76fad1 and at HEAD. All five match.
- I re-derived the IBLT success strings and their runs and permutation tests from the `iblt_trial` lines in `run.log`.
- I read the km-peeler source at b8ec39b (`ibflib.c`, `build_main.c`).
- I ran `ibltseq list` on the one saved extended-syncmer difference table (`kmer_poly_work/kmp_ext_Cannbio-2__Cannbio2_GigaDB/d.iblt`).
- I checked the aswritten run JSON, the verify JSON and the code paths.

**Prior rounds.** I ran `scripts/round history` for "syncmer", "km-peeler seed", "0.8m", "66", "unplaced" and "2.4"
(adjudications). Only "syncmer" (round 1 completeness, Gap 5) bears on a finding here, and it is addressed under F4. The
other searches returned no prior ruling on these points.

---

## F1: "the GigaDB sketch ran before x and λ were separated" (C11, second sentence)

**Premise tested, and it holds.**

- Commit 8c34af9 ("Draw x and lam independently, from two seeds") is an ancestor of e76fad1.
- All five hashed sources in `run.log` equal e76fad1.
- `run.log:28` (timed_giga) and `:40` (sweep) record `"seed_x": "0x5eed", "seed_lam": "0x1a3d"`.
- The 60 `kp10_seed` lines record `"seed_x": "os", "seed_lam": "os"` with distinct x.
- The paper's 11.3 is `kp10_dest.tsv:6` / `kp10_seeds.tsv:2` from this job. Neither comes from a pre-separation run.
- The as-written k = 31 and k = 63 runs also used two separate seeds (`run_aswritten_pair.sh:15`, and run JSON
  `seed_x 24301 = 0x5eed`, `seed_lam 6717 = 0x1a3d`).

The clause describes no run behind any number in the paper.

**Whole-set check (rule 3).** The sentence before it, "Our experiments used fixed seeds for reproducibility", is also
false for the 30 runs per pair, which drew both parameters from the OS. The attacker's alternative wording fixes the
clause but leaves that sentence contradicting it.

**Repair, endorsed, free.** Delete the clause. Then write: "Our experiments used fixed seeds (0x5eed for x, 0x1a3d for
λ) for reproducibility, except the 30 runs per pair of the set-sketch comparison, which drew both from the operating
system." This changes text only and needs no re-measurement.

**UPHELD.** What survives: the first sentence of C11 holds (kp10c hashes as written, and buckets are ranges of S). The
second sentence is false and should be deleted.

## F2: what the k = 63 sentence attributes to the numpy prototype (C3)

**Premise tested, and it holds.**

- **Who computed the hashes.** `aswritten.build` streams values from the C roller `hash/kp7e_stream`. The KMC truth
  is hashed by `s_of_strings`, which uses the C `kp7e_direct` and `kp7e_sym`. Numpy (`recon61`) only builds and
  decodes the sketch.
- **What the 20,000 check covers.** `kp7e_hash.self_test` hashes random 2-bit packed codes at `K` (31) in numpy and
  compares them with C. No 63-mer can pass through that path.
- **Number of buckets.** The run JSON says `bits: 10, bits_rule: "auto: bits_for(exact value D, m)"`. `bits_for`
  keeps at most m/2 differences per bucket. At D = 21,382, a 0.8m load gives 2^9 (41.8 per bucket). The m/2 rule
  gives 2^10 (20.9), which is the 548,864 bytes reported.

**Two corrections to the attacker.**

1. **The paper never says 0.8m for k = 63.** The 0.8m sentence is about the C whole-genome comparison. The k = 63
   sentence says only "buckets sized from the true difference". The defect is an ambiguity that a reader will resolve
   wrongly, not a false statement.
2. **The hashes the run used are checked.** "Nothing checks the numpy code the run did use except the final
   comparison with KMC" undersells that comparison: it checks all 21,382 values, the strongest check available. The
   k = 63 hash path is also checked by the aswritten self-test (`c5`: stream values equal the C reference on strings
   and are equal across orientation, at k = 63).

The attacker is right that "matching element by element" is a match in S-value space.

**Repair, endorsed, free.** Adopt the attacker's wording: "The k = 63 run hashed every k-mer with the C roller and
decoded with the numpy prototype, with 2^10 buckets from the true difference (at most m/2 per bucket). Its recovered
values equal S of KMC's 63-mer difference." Drop the 20,000 clause and its ledger row (ledger:115) from this sentence.
Optionally add "the roller matches the C reference at k = 63 in its self-test". Whether 2^9 buckets finish in one round
is untested, and the paper need not claim it.

**UPHELD** on the provenance error, narrowed on the bucket-rule leg as above. What survives: the k = 63 result
reproduces (21,382, 548,864 B, 1 round, 0 missed, 0 spurious), but the 20,000 check is about other code and the bucket
rule is m/2.

## F3: km-peeler's 30 "independent random seeds" share hash functions (C7, C5)

This is the attack's strongest finding, and I re-derived every link.

- **The seeds.** `kp10_vs_iblt_job.sh:176` passes `$((1000 + t))`. `build_main.c:85-91` passes `-s` unchanged to
  `ibf_sketch_init`.
- **The hash functions.** `ibf_hash_init` sets seed_i = root + i (`ibflib.c:27`, "FIXME use a better method?").
- **The table layout.** Cell = `hash_i.ls64b % chunk_size + i·chunk_size` (`ibflib.c:331, 421`), and every chunk
  has the same size. So trial t's hypergraph has the coordinates (h_t, h_{t+1}, h_{t+2}), and trial t+1 reuses two
  of those three columns exactly, only in different chunks. Adjacent trials are not independent. The 30 trials draw
  on 32 functions.
- **The data.** I recomputed the success strings at n = D:
  - draft: `100000000000111111000001000000` (8/30), runs 6 against 12.7 expected, z = −3.23, permutation
    p = 0.0017 (200,000 shuffles; the attacker's 0.0019 is within Monte Carlo error).
  - giga: no clustering (p = 0.66).
  - This is one post-hoc test on one of two strings. The mechanism, not the p-value, carries the finding.
- **The method card.** `../library/methods/pseudo_replication.md` governs: "Your effective sample size is the number
  of independent decisions, not the number of rows." The 10 trials with disjoint triples (seeds 1001, 1004, …, 1028)
  are a valid, conservative independent subset. 1 − 0.05^(1/10) = 25.9%, which matches the attacker.
- **KP-10's side is clean.** It used 30 OS draws of (x, λ).

**Scope.** C7 is false in one clause only ("independent random seeds", for km-peeler). Its other clauses hold:

- km-peeler is unmodified (`git status` shows only untracked logs).
- r = 3 is its default (`build_main.c:38`).
- It ran on canonical sets on one host.

The 30/30 counts in C1 and C6 remain true as counts. Only the 9.5% bound and the word "independent" fail for the
table. The same dependence weakens "sized to succeed in 30 of 30": the 1.05·D choice was read off correlated trials.

**Repairs.**

- **The attacker's text fix: endorsed, free.** Add "a valid bound for the table, from its 10 trials with disjoint
  hash triples, is 25.9 percent", or say that no bound is given for it. The paper's conclusion ("the runs do not
  compare reliability") only gets stronger.
- **Rerunning `iblt_sweep` with spaced or random seeds: endorsed as optional.** It costs a laila run and is owed,
  not required, because the text fix is honest without it.

**Correct minority.** As far as this adjudication can see, only this attacker found the seed derivation. The brief's
own C5/C7 wording ("seeds, host, shared inputs") steered away from it. No earlier round raised it ("km-peeler seed": a
searched absence).

**UPHELD.** What survives: KP-10's 9.5% stands, and the table's 9.5% and "independent random seeds" do not.

## F4: the extended-syncmer row (C8)

**Provenance legs: upheld.**

- The row comes from `full_kmp_ext_giga.json`, a 2026-11-10 run of `recon/kmpeeler_run.py`, not job b69c5c76.
- That run had 5 trials, not 30, and the paragraph above Table 4 says "each method ran 30 times".
- n was set from the true `D_ext` (`kmpeeler_run.py:164`).
- All 5 trials (seeds 2000–2004; 2000 and 2003 or 2004 share no hash function) give identical outcomes.

**The "listing incomplete" leg: OVERTURNED by a premise the attacker did not test.**

- In this harness `"complete"` is `listed == D_ext`, a line count compared with `comm -3`. It is not the IBLT's
  peeling status. `build_diff_list` discards the exit code.
- I listed the saved difference table from the last trial (seed 2004): `ibltseq list` **exits 0**, writes nothing to
  stderr, and emits 507,219 lines (507,099 `i,` + 120 `j,`), all of payload length 51.
- So the table peeled completely. Two entries of the `comm -3` difference never reached the IBLT as differences:
  they cancelled inside the table or were not distinct as 51-character payloads.

**Consequences for the attacker's argument.**

- Under the resend rule the k-mer table used (resend when `list` exits nonzero), this configuration needs no second
  round. The paper's "1 round" is right.
- "2.7 times smaller" compares against a listing that completed but is not exact. The paper's own sentence already
  says this ("but returned spurious syncmers numbering 61 percent… and missed 7").
- The attacker's proposed row, "no: listing incomplete in 5 of 5 trials", and the qualifier "at a size at which its
  listing did not complete" would both publish a false statement. I reject both.

**What the finding gets right, more strongly than it argued.** Round 1's completeness ruling (Gap 5, NARROWED) said
"the misses are explained by peeling". Peeling succeeded, so that explanation is now refuted. The 2-element shortfall
and the 7 missed k-mers are deterministic, from the input or the encoding. Whether that is a harness defect stays open,
as the ledger's carried-items row already records.

**Repairs.**

- **Free, endorsed.** Add a footnote to the Table 4 row: "from an earlier 5-trial run, with the table sized from the
  true extended-syncmer difference. Its listing completed in every trial but returned 507,219 of the 507,221
  differing extended syncmers." Do not change "1 round" or "no, 1,986,807 spurious and 7 missed".
- **Cheap but not tiny, so owed.** Diff the listing against `comm -3` of the two saved `.ext.txt` files (5.7 GB
  sequential) to name the 2 cancelled entries. Run it on bench or laila, not on the laptop. This has to be done
  before the 61%-against-3% comparison is used against the competitor.

**NARROWED.** What survives: the row's provenance (an earlier run, 5 trials, sized from the true D_ext) must be
disclosed, and round 1's peeling explanation falls. The listing completed, so "incomplete" and "needs another round"
are wrong.

## F5: the sketch's size parameters were chosen on the test pair (C7, C1)

**Premise tested, and it holds.**

- `prior.md:47-54`: four (buckets, m) settings on partition 0 of the GigaDB pair (D = 50,644, 1/64 of the pair). The
  run kept the one with the fewest bytes, 2^16 with m = 64.
- `kp10c.c:31-32` and the results file at `:71-72`: the 0.8m rule "reproduces the 2^16 and 2^8 buckets of the
  2026-11-10 runs".
- The window arithmetic checks out. 2^16 at giga needs a threshold t ≥ 49.41. 2^8 rather than 2^7 at draft needs
  t < 61.96. 0.8m = 51.2 sits inside that window.
- The older m/2 rule would have given giga 2^17. So 0.8m was chosen after the fact, and the paper discloses only the
  table's selection.

The attacker rightly notes that the power-of-two grid favours the table and does not reverse 20.8 against 38.5.

**Repair, endorsed, free.** "We chose m = 64 and 2^16 buckets in a calibration on one 64th of the GigaDB pair, and
set the 0.8m load so that it reproduces those bucket counts. Both methods' sizes were thus chosen on these pairs."

**UPHELD.** What survives: neither sizing rule is a priori. The sketch's choice was tuned for bytes on a slice of the
test pair.

## F6: "about m + c = 66 multiplications per element" (C12)

**The diagnosis holds.**

- The ledger (row 113) cites the numpy-era `prior.md:8`.
- In kp10c, `flag_roots` evaluates each side's degree-e polynomial at every element it holds in every differing
  bucket (`kp10c.c:226`, called at `:294-295`). Sizing keeps about 0.8m differences per bucket, so almost every
  bucket differs on any pair.

**The attacker's scope needs one correction.** The degree is that side's own root count. On giga, A holds 3,237,326 of
the differences and B holds 775. So A's elements pay about 49 extra multiplications (about 115 in all, 1.75 times 66),
and B's pay almost none. "Nearly every element pays about 49 more" is true of A only. The CPU split (49.5 s per side for
power sums, 110.3 s for decode summed over both sides, almost all of it A's) is consistent with roughly 2 times on A's
side.

**Repairs.**

- **The multiplication sentence: endorsed, free.** "About m + c = 66 multiplications per element to build the
  sketch, and up to about m more on the side holding the differences, to find the roots." Re-cite the ledger row to
  `kp10c.c` and `kp10_vs_shibuya_2026-11-11.md` §5.
- **The sentence "the decode took 110 CPU-seconds against 3.5 for the table's listing": OVERTURNED.** It picks the
  one phase that favours the table. In the same table, km-peeler's build costs 204.8 CPU-s per side against
  KP-10's 89.4, and the whole pipeline costs about 460 against 326. The results file itself says "This is not a speed
  claim for either algorithm". If the paper quotes times, it should quote the totals with that caveat.

**NARROWED.** What survives: 66 understates the per-element work by up to about m on the side holding the differences,
and the ledger cites a superseded source. The timing contrast the attacker proposes would mislead.

## F7: ledger rows (C15)

| item | ruling |
|---|---|
| 1 (20,000) | Upheld, per F2. Free fix. |
| 2 (66) | Upheld, per F6. Free fix. |
| 3 ("separated" has no row) | Upheld, per F1. Delete the sentence. Free. |
| 4 (0.8m) | **Overturned as a ledger defect.** The cited source, section 2 (`:71-72`), does say the rule "gives the 2^16 and 2^8 buckets of the earlier runs". The row's label "sizing rule" is accurate. The defect lies in the paper's Methods, already covered by F5. |
| 5 ("about 3 percent") | Upheld. There is no number row. The source exists (`shibuya_2022` summary, lines 27-29 → paper lines 782-783). Adding the row is free. |
| 6 (2.4×) | Upheld that the row crosses programs, against `paper_bench_2026-11-12.md:77` and the Figure 3 caption. The repair is under F8. |
| 7 (9.5% for both) | Upheld, per F3. Free. |

**NARROWED.** Six of the seven items stand, and item 4 is misfiled.

## F8: 2.4× crosses two harnesses (C14)

**The diagnosis holds.** 24.67 (bench_e) / 10.24 (bench_kept) = 2.41. The paper's own rule (Figure 3 caption;
`paper_bench:77`) forbids comparing across the two programs.

The size of the error is nil:

- 24.67/11.37 = 2.170 (bench_e against its own ntHash), and 11.43/10.24 = 1.116 (bench_kept).
- Chained, these give 2.42.
- The two harnesses agree where both measure the same thing (ntHash 11.37 against 11.43; sampler at 1/1000 9.43
  against 9.39). The 28.7 against 22.5 discrepancy is a per-kept slope at 1/4, not an every-k-mer cost.

**The attacker's second repair, "about 2.2 times (through each program's own ntHash)": OVERTURNED.** 2.2 is
Eisenstein against ntHash, not Eisenstein against the universal hash. The attacker's own line computes the chained
ratio as 2.17 × 1.12 ≈ 2.4 and then prescribes 2.2. Applying it would publish a wrong number.

**Repairs.**

- **Free, endorsed.** "2.4 times (each program's ratio to its own batched ntHash, 2.17 × 1.12)". Update the ledger
  row to this chain.
- **Costs one benchmark job.** Timing Eisenstein every-k-mer inside bench_kept is also acceptable, but it is not
  needed.

**NARROWED.** What survives: the provenance breaks the paper's stated rule. The number is right, and the attacker's
replacement number is wrong.

## F9: "almost all from unplaced scaffolds" is not counted (C10)

**Premise tested, and it holds.**

- No committed file assigns A∖B k-mers to records.
- The basis is record counts (147 against 10, `sets.tsv`) and bases. Cannbio-2's 914,397,426 bases minus the
  GigaDB release's 904,080,375 is exactly 10,317,051.
- B∖A = 775 shows that the chromosomes are not all identical.

The counts and KMC's element-by-element agreement survive (`truth.tsv`). KMC is an independent tool run in the same job.

**Repairs.**

- **Costs one pass over the genome, owed on bench or laila, not the laptop.** Count A∖B by Cannbio-2 record, as the
  attacker proposes.
- **Free, endorsed in the meantime.** State only what is measured: "differ in 3,238,101 canonical 31-mers, 0.73
  percent of the set. The release keeps the 10 chromosomes and drops 137 unplaced scaffolds (10.3 Mb)." Drop "almost
  all from" until the count exists. The attacker's alternative, "most of it on the 137 unplaced scaffolds", still
  asserts the attribution and needs the same count.

**UPHELD.** What survives: the counts and KMC's agreement stand, and the attribution is unmeasured.

## The attacker's "Other claims"

I spot-checked C4's arithmetic and it holds: 36,481,232 − 35,127,296 = 1,353,936 = 2 × 1,263 × 536. These rulings
on C2, C4, C6 and C9 raise no finding, so they get no row.

The attacker noted that `hash/kp7seed.h` is missing from the hashed list. That is a minor provenance gap: adding it to
`run.log`'s hash list in future jobs is free, and no number changes.

## Verdict table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | x and λ were separate in every run behind the paper's numbers. Delete the clause and also amend "our experiments used fixed seeds" (the 30 runs per pair were OS-seeded); free |
| F2 | UPHELD | k = 63 was hashed by the C roller, the 20,000 check is k = 31 packed codes, and the buckets used m/2 (2^10). The paper never stated 0.8m for k = 63. Attacker's rewording endorsed; free |
| F3 | UPHELD | km-peeler uses seeds root+0,1,2 with equal chunks, so adjacent trials share two hash columns (draft runs p = 0.0017). The table's 9.5% and "independent" fail and KP-10's stand. Text fix free; reseeded rerun owed (laila) |
| F4 | NARROWED | provenance holds (earlier 5-trial run, n from the true D_ext, deterministic outcome, round 1's peeling explanation falls), but `ibltseq list` exits 0 on the saved table, so "listing incomplete / needs another round" is false. Footnote free; tracing the 2 entries owed |
| F5 | UPHELD | 2^16 and m = 64 were calibrated for bytes on 1/64 of the GigaDB pair, and 0.8m was set afterwards to reproduce them. Disclosure free |
| F6 | NARROWED | root flagging adds up to about m multiplications per element on the side holding the differences, and the ledger cites a numpy-era source. Wording fix free; the proposed 110-vs-3.5 s contrast is cherry-picked (totals 326 vs about 460 CPU-s) and rejected |
| F7 | NARROWED | items 1, 2, 3, 5, 6 and 7 stand. Item 4's source does disclose that the rule reproduces the earlier counts, so it belongs to F5. Fixes free |
| F8 | NARROWED | the cross-program provenance breaks the paper's rule, but the chained ratio is 2.42, so 2.4 is right. The attacker's "about 2.2" is wrong; write "2.4 (2.17 × 1.12)"; free |
| F9 | UPHELD | the scaffold attribution was never counted. State the measured records and bases (free), or count A∖B per record (one pass, owed) |
