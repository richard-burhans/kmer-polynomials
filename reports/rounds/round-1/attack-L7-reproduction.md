# Round 1 attack, slot L7, lens: reproduction

Attacker: L7 (reproduction). Target: `reports/rounds/round-1/ATTACKED-paper.md` (= manuscript/paper.md at 2dea163).

What I re-ran, not just re-read (all local; laila timings were re-derived from the committed raw counters, not re-timed):

- `prototype/kp7_symmetric_check.py` -> byte-identical to `results/kp7_symmetric_check_2026-11-10.txt`.
- A fresh copy of `hash/` (in `bench/data/kmer_poly_work/round1_L7/hash`): `test_kp7`, `test_kp7m31`, `test_kp7pack`,
  `test_kp7e` all pass; `make mutation`, `mutation_pack`, `mutation_e` kill all 14 mutants. This machine has no
  AVX-512, so `test_kp7m31` skips impl 2 locally (it was run inside an earlier laila job, `results/kp7_m31_simd_2026-11-10.md:32`).
- `ideas/eisenstein` spaced-seed tests and mutant: pass / killed.
- `tests/` (26 pytest tests) pass; reverting the Padé negative-index fix in `recon/recon61.py` or `prototype/recon.py`
  makes exactly the matching one-sided test fail.
- My own exhaustive Eisenstein enumeration over F_q, q = 2^61 - 1, two seeds (`tmp/L7/eis.py`).
- `manuscript/paper_bench_table.py results/paper_bench_2026-11-11` -> identical to the committed table.
- A full-genome element-by-element ntbatch-vs-library check (`round1_L7/nteq.cpp`).
- `ideas/eisenstein/kp1_tstv.py` on the committed pair and on four more pairs of the same WF-I MAF directory.
- The repository's own `recon61.decode` on toy fields (`tmp/L7/overflow.py`).
- `ideas/kp5_kp11/ir` (palindromes) and `mismatch` (KP-11), `ideas/eisenstein/dbg_e31` (KP-3), `hash/bench --nthash-classes`.
- A direct SeedNtHash reverse-complement test and a GF(2)-linearity test of the ntHash library (`round1_L7/seedrc.cpp`, `gf2.cpp`).

`scripts/round history` reports no closed round (this is round 1), so no finding below re-files an adjudicated point.
The repository holds no method cards (searched `--include=*.md` for "method card"; none outside the round briefs).

## F1 -- the abstract's "1/20 and below" fails at 1/20 by the paper's own fit, and the break-even is not resolved by two repeats

Claim C11. The break-even comes from a straight line through two batched-sampler rows (1/16, 1/4) of job
08653d9e6a7453a2; no rate between 1/16 and 1/100 was measured. I recomputed it from `results/paper_bench_2026-11-11/counters.tsv`
and `run.log` (`tmp/L7/be.py`), using the measured kept fractions (45,892,012 and 183,856,486 of 736,547,799 input k-mers):

- Mean fit: 28.53 cycles per kept k-mer, intercept 9.96. At exactly r = 1/20 the model gives **11.388 cycles against
  batched ntHash's 11.365**, so at the rate the abstract names, the sampler is predicted to be *slower*. The break-even is
  1/20.3 (`results/paper_bench_2026-11-11.md:47` says so itself).
- Pairing the individual repeats (ntHash 11.327 / 11.403; 1/16 11.691 / 11.787; 1/4 17.059 / 17.107) moves the break-even
  over **1/19.1 to 1/21.7**. Two repeats cannot place it at 1/20 versus 1/22.
- The measured rows are: 1/16 slower (0.97x), 1/100 faster (1.14x branchy, 1.11x batched).

Correction: abstract "faster than batched ntHash at sampling rates of 1/20 and below" -> "faster than batched ntHash at
1/100 and sparser (measured), with an interpolated break-even near 1/20". The Results' "about 1/20" is defensible only
with the interval stated. Verdict **NARROWED**.

## F2 -- "Every timing in the paper comes from one job" and "from hardware performance counters" are false for the paper as written

Claim C10. The 18 table numbers reproduce exactly (`paper_bench_table.py` output is identical to
`results/paper_bench_2026-11-11/table.md`; spreads <= 1.3% for every row in the paper). But Methods says "Every timing
in the paper comes from one job, in which each method ran in its own counted process", and the paper quotes timings that
did not:

- Spaced seeds, 0.76x: `bench_spaced_e`'s own wall-clock rates, three rounds inside ONE process, not counters and not
  one process per method (`results/paper_bench_2026-11-11.md`, "Spaced seeds (run_51; bench_spaced_e's own wall-clock
  rates"). It is also a time ratio over different hash counts (Eisenstein 42,088,079 hashes in 1.953 s, SeedNtHash
  42,667,269 in 1.482 s): per hash it is 21.55/28.79 = 0.75x, not 0.76x.
- Supplement S1: "four to six times slower" (KP-3) is laila job af1eb52f65c6d471 (`results/ideas_eisenstein_2026-11-10.md:7, 75-83`);
  "1,370 times faster" (KP-11) and the palindrome speed statement (KP-5) are job 9f793f18b3112141
  (`results/ideas_kp5_kp11_2026-11-11.md:22, 74`). The ledger's own rule (manuscript/ledger.md:5-6) is that numbers from
  different jobs must not be compared because between-job variance was 3-6%.
- KP-10 "about an hour per pair" is a numpy prototype, not laila.

Correction: "Every timing in the main text's speed table comes from one job..." and label the spaced-seed and
Supplement timings with their own jobs and method (wall clock). Verdict **NARROWED**.

## F3 -- the bit-for-bit claim for the batched ntHash has no committed evidence, though it does hold

Claim C12. The only check in the code is `hash/bench_pack.cpp:274-296`: an XOR-of-all-hashes plus count over the genome,
and element-by-element (pos, hash) on the first 2,000,000 k-mers of the first record only. No result file records its
output (grep for "positions and hashes identical" / "identical to ntHash" over `results/` finds nothing), and the paper
benchmark job (`manuscript/run_paper_bench.sh`) runs `bench_e`, which does not compare ntbatch with the library at all.

I closed the gap: `round1_L7/nteq.cpp` compares `ntbatch::hash_seq` against `nthash::NtHash::roll()` / `hashes()[0]`
(ntHash c26bd45, the pinned library), position and value, at every k-mer of cs10: **221 records, 736,547,799 k-mers on
both sides, 0 mismatches**. 23.53/11.37 = 2.07, so "about twice" holds. Verdict **SURVIVED** (by my check, not by the
repository's; the paper should cite a committed full-genome check).

## F4 -- the Ts/Tv range is one plant's two haplotypes, and mixes an all-columns row with a repeat-masked row

Claim C9. `kp1_tstv.py` reproduces `results/ideas_eisenstein/kp1_tstv_NC_044378.1.json` exactly (JSON-equal). But:

- The "pairwise alignments of cannabis assemblies" are two MAFs, KOMPa and KOMPb, the two haplotypes of ONE plant, out of
  216 in the same directory (`bench/data/wfi_laila/NC_044378.1/bare/`); the third row is not an alignment but a
  projection of the two through cs10 coordinates.
- 1.58 is an all-columns row (cs10 vs KOMPa/KOMPb); 1.85 is the uppercase-only (outside soft-masked repeats) row of the
  projected KOMPa-vs-KOMPb comparison. The all-columns range is 1.58-1.67 (61.3-62.5%).
- The range does not survive other assemblies. Running the same script on four more pairs from the same directory:

| comparison | all columns Ts/Tv (transitions) | uppercase only |
|---|---|---|
| cs10 vs ACDC | 1.511 (60.2%) | 1.635 (62.1%) |
| cs10 vs ANC-1 | 1.546 (60.7%) | 1.695 |
| cs10 vs ASM151000v1 | 1.745 (63.6%) | 1.871 (65.2%) |
| ASM1303036v1 vs ASM151000v1 | 1.810 (64.4%) | **1.936 (65.9%)** |

Across these the range is Ts/Tv 1.51-1.94, 60-66% transitions. Correction: "In alignments of cannabis assemblies to
NC_044378.1, 60 to 66 percent of substitutions are transitions (Ts/Tv 1.5 to 1.9)", or keep 61-65% and say it is the two
haplotypes of one plant, with and without soft-masked repeats. The qualitative point (transitions are the commonest
class) stands. Verdict **NARROWED**.

## F5 -- overflow detection is NOT deterministic: the repository's decoder accepts wrong answers for overflowing buckets

Claims C16 (and the Discussion's "detected its own overflow deterministically"). `recon/recon61.py:100-136` decodes from
S_1..S_m, flags roots among each side's own elements, then checks ONLY S_{m+1}..S_{m+c} and the count against the
flagged elements. It never re-checks S_1..S_m. Whether a wrong answer passes depends on the hash values, so detection is
probabilistic over the hash.

Demonstration with the repository's own `decode` (I set `recon61.Q` to a small prime and replaced only the vectorized
`_flags` by the same Horner loop mod Q; `tmp/L7/overflow.py`), on random buckets that ALL hold more than m differences
plus 5 shared elements, m = 2:

| q | c | trials | detected | accepted with a WRONG answer |
|---|---|---|---|---|
| 31 | 1 | 50,000 | 49,833 | **167** |
| 31 | 2 | 50,000 | 49,994 | **6** |
| 61 | 1 | 50,000 | 49,921 | **79** |
| 101 | 1 | 50,000 | 49,958 | **42** |
| 101 | 2 | 50,000 | 50,000 | 0 |

A concrete one, q = 31, m = 2, c = 1: true A\B = {18, 29, 30}, B\A = {13, 22, 23}; the decoder returns an EMPTY
difference and accepts it, although S_1 = 19 and S_2 = 15 are nonzero. It passes because S_3 = 0 and the counts balance.
With c = 2: true {10, 22, 26} / {3, 9, 11, 19} accepted as {} / {15}, leaving S_2 = 14 unexplained.

The rate scales as about 1/q^c, so at q = 2^61 - 1 and c = 2 it is negligible (order 10^-37 per overflowing bucket), and
no real run is in doubt. But "deterministically" is false, and the decoder throws away a free check: comparing the flagged
elements against all m + c sums would cost nothing. Correction: "is detected except with probability about 1/q^c per
overflowing bucket" (Methods), and drop "deterministically" from the Discussion and from the table's "deterministic
overflow" row (results/kp10_vs_shibuya_2026-11-11.md:97, 124). Verdict **REFUTED** (the wording; the practical claim
survives).

## F6 -- "of equal reliability" compares a two-round sketch with a one-shot IBLT

Claim C15. Every number reproduces from `results/kp10_vs_shibuya_2026-11-11/*.json`: 36,494,096 bytes / 3,238,101 =
11.27; IBLT 4,158,506 / 202,557 = 20.53; extended syncmers 13,667,246 / 3,238,101 = 4.22, 2.67x smaller; 1,986,807
spurious, 7 missed. Bytes are accounted the same way for both (one side's sketch; KP-10 includes its round-2 resends:
4,561,762 elements = 65,536 x 67 + 1,275 x 134). But the regimes are not equal:

- KP-10's round 1 was NOT exact: 1,275 of 65,536 buckets overflowed and were split and resent (`rounds_max: 2`,
  `full_kp10_eis_giga.json`). The 11.27 is a two-round, interactive protocol. The IBLT's 30/30 at 1.05 x D is one round.
  An IBLT also knows when it fails to peel and could likewise be resent; the IBLT was given no second round.
- KP-10's "exact in every run" is one run per field at whole-genome scale (two runs); the IBLT's 30/30 is 30 trials on
  partitions 0-3 only (D = 202,557), then extrapolated x16 ("66 MB ... extrapolated", results md:95).
- At the same 1.05 x D the IBLT failed 4 of 30 on the draft pair, so 20.5 B is "equal reliability" only on those four
  GigaDB partitions.
- Both sizes were set from the TRUE difference (IBLT n = 1.05 x D; KP-10's 2^16 x m = 64 calibrated on partition 0 with
  D known, and 2^8 buckets for the draft pair), which the paper does not say.

Correction: "11.3 bytes per difference over two rounds, against 20.5 for a one-round invertible Bloom lookup table that
succeeded in 30 of 30 trials on 4 of 64 partitions, both sized from the true difference." Same for the Discussion's
"0.55 times the bytes of an invertible Bloom lookup table of equal reliability". Verdict **NARROWED**.

## F7 -- ntHash's strand-symmetric value is not GF(2)-linear; only its per-strand hashes are

Claim C18. The paper: "Constant-time rolling, constant-time single-base substitution and hashing of concatenations ...
all hold for ntHash too, which is linear over GF(2)." ntHash's canonical value is fwd + rev mod 2^64
(`hash/ntbatch.h:75`, `buf[n++] = f + r`, bit-identical to the library per F3), an integer sum with carries. Test
(`round1_L7/gf2.cpp`, ntHash c26bd45): apply the same substitution (A -> G at a random position) in two random
backgrounds and compare the XOR differences, 10,000 trials:

- forward hash: equal 10,000/10,000; reverse hash: 10,000/10,000 (GF(2)-affine, as claimed for the strand hashes);
- **canonical hash (what is compared with ours): equal 0/10,000.**

The operations do hold for ntHash via its two strand states, so the conclusion "not advantages" is right, but the stated
reason is wrong for the value being compared. Correction: "which hold for ntHash's two strand hashes, each linear over
GF(2)". The 8-mer test exists and reproduces: `hash/bench --nthash-classes` -> "32896 distinct values over 65536 k-mers
(32896 double-stranded classes)". Verdict **NARROWED**.

## F8 -- the CM31 bound 1.4e-17 is wrong, and the stated 2(k-1)/q is 2.6e-17, not 2.7e-17

Claim C1 (bound arithmetic). Methods states the bound as 2(k-1)/q and then "At k = 31 that is about 2.7e-17". 2(k-1)/q =
60/(2^61 - 1) = **2.60e-17**; 2.7e-17 is 2k/q = 2.69e-17 (the abstract's form). Pick one. Implementation: "its collision
bound is about 1.4 x 10^-17 per pair at k = 31" for F_{p^2}, p = 2^31 - 1: 2k/p^2 = 62/(2^31-1)^2 = **1.34e-17** and
2(k-1)/p^2 = 1.30e-17; neither rounds to 1.4 (1.39e-17 is 64/p^2). The error originates in
`results/kp7_m31_simd_2026-11-10.md:13` and propagates to the ledger (line 17) and `results/kp10_vs_shibuya_2026-11-11.md:43,129`.

The rest of C1 reproduces: the exhaustive k = 7/8 counts (8,192 / 32,896) reproduce byte-identically for the fourth-root
script and in `test_kp7e` for Eisenstein; my independent Eisenstein enumeration over F_q gives F injective (16,384 /
65,536), S = 8,192 / 32,896 at two different (x, lambda) seeds. Total degree 2(k-1) is right (lambda^2 cancels; lambda(F+R)
has degree <= k <= 2(k-1) for k >= 2). Verdict **NARROWED** (arithmetic only; the injectivity argument survives my checks).

## F9 -- testing claims hold on re-run, with two overstatements

Claim C20. All test binaries pass and every named mutant dies (14 C mutants: ROLL, LANE, TAB among them; spaced-seed
mutant; both Padé reversions each fail exactly the one-sided test). The pre-fix Gaussian draft run is 3 rounds, 257 splits
(`full_kp10_gau_draft.prefix.json`), the post-fix one 1 round. Overstatements:

- "block sizes from 1 to 4,096" reads as a range; the tests use discrete sets ({1, 3, 5, 64, 4096} in `test_kp7e.c:149`,
  {1, 3, 64, 4096} in `test_kp7pack.c`; {1, 3, 5, 7, 8, 9, 15, 16, 17, 64, 4096} only for CM31).
- "Every variant is checked": the paper-benchmark job built `bench`/`bench_e` and ran no test (`run_paper_bench.sh`
  `--build` line), so the AVX-512 binary that produced 25.86 was not tested in that job; its only AVX-512 test run is an
  earlier job. Locally (no AVX-512) the 16-lane path is skipped.

Correction: "block sizes 1, 3, 5, 64 and 4,096". Verdict **SURVIVED** (the substance reproduces; the wording should be fixed).

## F10 -- C17's orientation cause, checked against the assemblies

Claim C17. 137,216 bytes = 17,152 elements x 8 = 2^8 buckets x 67, one round (`full_kp10_eis_draft.json`): reproduces.
48,293,748 is in the results md (line 110); I did not rerun km-peeler (it died of memory). The abstract's "saw 48 million"
counts extended syncmers (51-mers) against 7,931 k-mers: different units, which the abstract does not say.

I tested the cause directly (`tmp/L7/orient.py`): 31-mers sampled by a strand-independent content filter (mix64 of the
canonical code, 1 in 1,024), kept only if single-copy in Cannbio-2 (323,937), each recording which strand it is written
on; each draft contig is then called same / flipped when at least 90% of its sampled k-mers agree / disagree with
Cannbio-2. Result over the draft's 8,477 contigs: **4,110 same (48.7% of bases), 3,971 flipped (45.4%), 168 mixed (5.5%),
228 with no sampled k-mer (0.4%)**. With ~55 M extended syncmers per assembly, flipping ~45-50% of the sequence predicts a
strand-specific symmetric difference of roughly 2 x 0.45-0.5 x 55 M = 50-55 M; 48.3 M is in line with that. So the
48 M IS an orientation effect. The abstract says "a strand-specific method", which is accurate, though it counts a
different object. Verdict **SURVIVED** (with the units caveat: "48 million differing extended syncmers").

## F11 -- smaller reproduction notes on the remaining claims

- C2: the table reproduces byte-identically from `prototype/kp7_symmetric_check.py`, which computes over Z[i]/(q) =
  F_{q^2}, q = 2^31 - 1. The caption's "q = 2^31 - 1" omits that the field is F_{q^2} and reuses the letter q, which the
  paper defines as 2^61 - 1. The F.R count is encoding-specific: under the Eisenstein encoding over F_q I get 4,160 at
  k = 7 but **16,500** at k = 8 (fourth-root: 16,362). Both are below the 16,512 reverse/complement orbits, so F.R merges
  more than plain reverses. SURVIVED; fix the caption to "F_{p^2}, p = 2^31 - 1".
- C3: both recurrences hold algebraically and `test_kp7e` checks rolling == direct at every window for k = 31, 8, 1.
  SURVIVED.
- C4: q mod 3 = 1, q mod 4 = 3 (computed); 47.37/24.69 = 1.92. SURVIVED.
- C8: `test_kp7e` prints "transitions 1, A-C/G-T 3, complements 4". SURVIVED.
- C13: recomputed from both TSVs: mean max bit deviation 0.00090 (Eisenstein) / 0.00088 (ntHash); max density deviation
  0.40% (ntHash, w = 50); long gaps 341 + 363 + 371 + 332 + 367 = 1,774, outside repeats 2; (19/21)^200 x ~70.1 M samples =
  0.14. One overstatement: "Every method left tens of ... gaps per seed" -- ntHash seed 0 left **8**
  (`kp7_stats_seeds_2026-11-10.eisenstein.tsv`, ntHash row seed 0). NARROWED (wording: "all but one method-seed left tens").
- C14: numbers reproduce from the job log (0.74722 = exact, SeedNtHash 0.81794). I added the missing direct test
  (`round1_L7/seedrc.cpp`): SeedNtHash(s) == SeedNtHash(rc s) at **0 of 199,970** windows for the non-palindromic mask
  1101100111010111011100101110101, and 199,970/199,970 for a palindromic mask and for all-ones. ntHash2 prints the
  warning in the log. SURVIVED; the paper should state the 0/199,970 result, since a reader sees 0.818 > 0.747 as better.
- C21: 13,608,457 edges reproduce (`dbg_e31`, 0 false, 0 missed), but on "the first 20 Mb of cs10", a prefix of
  NC_044371.1 (about 101 Mb), not "a cs10 chromosome's graph"; locally the slowdown is 6.8x (1,004 vs 148 ns/node), outside
  "four to six"; laila gave 5.9x. 7,760 palindromic 20-mers reproduce (`ir ... palin 20`, hash and pack, same checksum).
  KP-11 at k = 12, d = 3 reproduces exactly (`mismatch ... 12 3 both`: identical sums), but the speed-up is not portable:
  locally, two runs give direct 172.4 s / 156.7 s against NTT 0.361 s / 0.383 s, i.e. **478x and 409x**, against the
  paper's 1,370x from ONE laila run (direct 1,108 s, `results/ideas_kp5_kp11_2026-11-11.md:74`). The ratio depends on
  how slow the direct baseline runs on the host (6.6x slower on laila than here), so "1,370 times faster" should read
  "400 to 1,400 times faster, depending on the machine" or carry its single-run, single-host qualifier. NARROWED.
- C5, C6, C7, C19: priority and literature readings; nothing in them can be re-run. C5's identity and "no fourth-root map"
  are elementary (A perpendicular to C gives Re(A conj C) = 0). C5 SURVIVED on the algebra; C6, C7, C19 UNCHALLENGED by
  this lens.

## Verdicts

| C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 | C9 | C10 | C11 | C12 | C13 | C14 | C15 | C16 | C17 | C18 | C19 | C20 | C21 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| NARROWED | SURVIVED | SURVIVED | SURVIVED | SURVIVED | UNCHALLENGED | UNCHALLENGED | SURVIVED | NARROWED | NARROWED | NARROWED | SURVIVED | NARROWED | SURVIVED | NARROWED | REFUTED | SURVIVED | NARROWED | UNCHALLENGED | SURVIVED | NARROWED |
