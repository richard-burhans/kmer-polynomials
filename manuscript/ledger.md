# Numeric ledger (pass 1)

Every number the manuscript uses, with its committed source. Paths are relative to the repository root, and
`file:line` points at the line that states the number. "Job" is the laila job that produced a timing. Numbers from
different jobs must not be compared directly in the paper, because between-job variance was 3–6%
(results/kp7_kept_eval_2026-11-10.md:46–47). The consolidated-benchmark list is at the end.

## Algebra and exactness (deterministic, no job)

| number | meaning | source |
|---|---|---|
| 729 / 6,561 | distinct values of F + R on all 7-mers / 8-mers (Gaussian encoding) | results/kp7_symmetric_check_2026-11-10.txt:2, 7 |
| 4,160 / 16,362 | distinct values of F·R, k = 7 / 8 (Gaussian) | results/kp7_symmetric_check_2026-11-10.txt:3, 8 |
| 8,192 / 32,896 | distinct values of (F+λ)(R+λ), k = 7 / 8, equal to the number of double-stranded classes (Gaussian) | results/kp7_symmetric_check_2026-11-10.txt:1, 4, 6, 9 |
| 8,192 / 32,896 | the same for the Eisenstein hash, with F taking 4^k distinct values on the k-mers | results/kp7_eisenstein_2026-11-10.md:33–34 |
| 2(k−1)/q ≈ 2.6 × 10⁻¹⁷ | collision bound of the Eisenstein symmetric hash per pair, k = 31, q = 2⁶¹ − 1 (corrected in round 1; was 2.7, which is 2k/q) | results/kp7_eisenstein_2026-11-10.md:29–30 |
| 2(k−1)/p² ≈ 1.3 × 10⁻¹⁷ | collision bound of the CM31 (x complex) hash per pair, k = 31 (corrected in round 1; was 1.4) | results/kp7_m31_simd_2026-11-10.md (Method paragraph) |
| 1 / 3 / 4 | squared distances in ℂ: transitions, A–C/G–T transversions, complements (Eisenstein) | results/kp7_eisenstein_2026-11-10.md:35–36 |
| 3 multiplies per k-mer | Eisenstein rolling cost | results/kp7_eisenstein_2026-11-10.md:19–20 |
| 16 KB | Eisenstein v2 table size (8-byte entries) | results/kp7_eisenstein_2026-11-10.md:22 |
| θ = π/3 | Das, Nanda & Das 2017 fix θ explicitly | literature/eisenstein_wantlist.md (Resolution section) |
| eq. (5), eq. (4) | Cristea 2002: the fourth-root map (up to sign) and the θ = π/4 point | library papers_md/summaries/cristea_2002_genomic-signals.md:25–30; literature/eisenstein_wantlist.md:103–104 |
| rationale echoes Cristea's eq. (3) remark | Das 2017's printed rationale ("purines -G and pyrimidines-C have identical imaginary" parts, false at θ = π/3) + complements sum to zero, no distances, θ asserted; Cristea's remark sits under his conjugation map (round-1 adjudication L4 F13) | library papers_md/summaries/das_2017_trigonometric-mapping.md:23–27, 63–64 |
| sensitivity = specificity formula; 3 exon ranges | Das 2017 reliability caveats | library papers_md/summaries/das_2017_trigonometric-mapping.md:34–39, 48–51 |
| A and G swapped, complement is conjugation | Das, Nanda & Das 2019, the mapping's source cited by Gulocak 2022 | library papers_md/summaries/das_2019_adaptive-kaiser-window.md |
| binary codes with transition-aware distance | Yu, Li & Yu 2018 (narrows the novelty claim) | library papers_md/summaries/yu_2018_encoding-survey.md:24–31 |

## Speed (laila, Xeon Gold 6330, cycles per input k-mer, cs10, k = 31)

Every timing in the paper comes from ONE job, c84a395ba84b315e, at ONE commit (5bee54c), median of 5 interleaved repeats, run by
`manuscript/run_paper_bench.sh`; the table is `manuscript/paper_bench_table.py` on its raw counters
(results/paper_bench_2026-11-12/). Input: 736,547,799 k-mers (results/paper_bench_2026-11-12.md:14).

| number | meaning | job | source |
|---|---|---|---|
| 11.37 | batched ntHash (bench_e) | c84a395ba84b315e | results/paper_bench_2026-11-12.md:30 |
| 23.92 (0.48×) | ntHash library roll() | c84a395ba84b315e | results/paper_bench_2026-11-12.md:31 |
| 24.67 (0.46×) | Eisenstein, every k-mer, 2 lanes | c84a395ba84b315e | results/paper_bench_2026-11-12.md:33 |
| 27.85 (0.41×) | Eisenstein, every k-mer, 1 lane | c84a395ba84b315e | results/paper_bench_2026-11-12.md:32 |
| 47.19 (0.24×) | Gaussian 61-bit, every k-mer, 2 lanes | c84a395ba84b315e | results/paper_bench_2026-11-12.md:34 |
| 26.44 (0.43×) | CM31 AVX2, every k-mer, 8 lanes | c84a395ba84b315e | results/paper_bench_2026-11-12.md:35 |
| 25.94 (0.44×) | CM31 AVX-512, every k-mer, 16 lanes | c84a395ba84b315e | results/paper_bench_2026-11-12.md:36 |
| 9.43 (1.21×) | Eisenstein sampler, branchy, 1/1000 | c84a395ba84b315e | results/paper_bench_2026-11-12.md:37 |
| 9.94 (1.14×) | Eisenstein sampler, branchy, 1/100 | c84a395ba84b315e | results/paper_bench_2026-11-12.md:38 |
| 11.74 (0.97×) | Eisenstein sampler, batched, 1/16 | c84a395ba84b315e | results/paper_bench_2026-11-12.md:39 |
| 17.12 (0.66×) | Eisenstein sampler, batched, 1/4 | c84a395ba84b315e | results/paper_bench_2026-11-12.md:40 |
| 28.7 | Eisenstein cycles per kept k-mer, bench_e (slope 1/16 to 1/4) | c84a395ba84b315e | results/paper_bench_2026-11-12.md:42 |
| 1/20.2 | Eisenstein sampler break-even against batched ntHash, bench_e | c84a395ba84b315e | results/paper_bench_2026-11-12.md:42 |
| 2.1× ("about twice") | library roll() / batched ntHash, 23.92 / 11.37 | c84a395ba84b315e | results/paper_bench_2026-11-12.md:30–31 |
| 11.43 | batched ntHash in bench_kept (reference for the second table) | c84a395ba84b315e | results/paper_bench_2026-11-12.md:50 |
| 10.24 (1.12×) | canonical code + universal hash, every k-mer | c84a395ba84b315e | results/paper_bench_2026-11-12.md:51 |
| 1.22 / 1.15 | Eisenstein sampler, branchy, 1/1000 / 1/100, bench_kept | c84a395ba84b315e | results/paper_bench_2026-11-12.md:52 |
| 1.18 / 1.12 | gate + ntHash on kept, branchy, 1/1000 / 1/100 | c84a395ba84b315e | results/paper_bench_2026-11-12.md:55 |
| 1.25 / 1.20 | gate + universal hash on kept, branchy, 1/1000 / 1/100 | c84a395ba84b315e | results/paper_bench_2026-11-12.md:56 |
| 22.5 / 19.3 / 5.9 | cycles per kept k-mer, bench_kept: Eisenstein / gate + ntHash / gate + universal hash | c84a395ba84b315e | results/paper_bench_2026-11-12.md:60 |
| 22.5 vs 28.7 | per-kept cost of the same sampler code in the two programs, unexplained | c84a395ba84b315e | results/paper_bench_2026-11-12.md:70–77 |
| 2.4× | 24.67 / 10.24: universal hash every k-mer vs Eisenstein 2 lanes (Discussion; across the two programs) | c84a395ba84b315e | results/paper_bench_2026-11-12.md:33, 51 |
| 2.2× | 24.67 / 11.37: every-k-mer Eisenstein time over batched ntHash (abstract) | c84a395ba84b315e | results/paper_bench_2026-11-12.md:30, 33 |
| 2.1–3.1 GHz | measured core clock range over all 165 runs (cycles / ref-cycles × 2.0) | c84a395ba84b315e | results/paper_bench_2026-11-12.md:counters.tsv |
| median of 5; 6 load runs | statistic and load baselines | c84a395ba84b315e | results/paper_bench_2026-11-12.md:17–24 |
| 0.75× per hash (0.76× wall-clock) | Eisenstein spaced seeds vs SeedNtHash, one mask, median of 5 rounds | c84a395ba84b315e | results/paper_bench_2026-11-12.md:83 |
| 0.747 / 0.818 | windows unchanged at 1% substitutions, ours vs SeedNtHash | c84a395ba84b315e | results/paper_bench_2026-11-12.md:86 |

Superseded (job 08653d9e6a7453a2, mean of 2, results/paper_bench_2026-11-11.md): 23.53 / 24.69 / 27.82 / 47.37 / 26.43 /
25.86 / 9.40 / 10.01 / 17.08 (0.67×) / 28.5 / 1/20.3, and the CM31 per-kept 64.3 (no longer quoted). Earlier: 11.41 /
23.44 / 24.59 / 27.95 / 47.51 / 25.83 / 9.41 / 9.97 / 11.87 / 17.21 from 7fc2a700073c1f4f, and the break-even "≈ 1/19".

## Sampling quality (local, full cs10, 5 seeds per method, no job)

| number | meaning | source |
|---|---|---|
| 0.14 | geometric expectation of gaps > 200 k-mers per seed, rate 2/21 | results/kp7_next_steps_2026-11-10.md:50–51 |
| ~~1,774 long gaps, all but 2 in repeats~~ SUPERSEDED (round 1): 1,583 gaps at 151 loci, 147 in repeats, after reseeding | across 5 methods × 5 seeds | commit 0db2b6e; results/kp7_stats_seeds_2026-11-11.gaps.tsv |
| ~~367, 0 outside repeats~~ SUPERSEDED: dropped from the paper (a zero with no power; after reseeding 2 Eisenstein gaps lie outside repeats) | Eisenstein provable hash's long gaps | results/round1_owed_2026-11-11.md |
| 0.0009 to 0.0010 (0.00088–0.00103) | mean max per-bit deviation from ½, per method | results/kp7_stats_seeds_eisenstein_2026-11-11.md:19–23 |
| 0.00090 / 0.00088 | Eisenstein provable / ntHash mean max bit deviation | results/kp7_stats_seeds_eisenstein_2026-11-11.md:19, 23 |
| within 0.4% (max 0.40%, ntHash) | minimizer density deviation, w = 10, 20, 50, every method | results/kp7_stats_seeds_eisenstein_2026-11-11.md:19–23 |
| 47.7% | soft-masked fraction of cs10 | results/kp7_next_steps_2026-11-10.md:50 |

## Spaced seeds (frozen job c84a395ba84b315e; numbers unchanged from 08653d9e6a7453a2)

| number | meaning | source |
|---|---|---|
| 0.75× per hash (0.76× wall-clock; median of 5 rounds) | Eisenstein spaced-seed hash against ntHash2 SeedNtHash, one mask, first 60 Mb of cs10 (was 0.78× in job af1eb52f65c6d471) | results/paper_bench_2026-11-12.md:83 |
| 0.747 (0.74722) | fraction of windows unchanged at 1% SNPs, non-palindromic mask: exact double-stranded value, matched by ours | results/paper_bench_2026-11-12.md:86 |
| 0.818 (0.81794) | SeedNtHash's value on the same test | results/paper_bench_2026-11-12.md:86 |

## Transition frequency (KP-1)

| number | meaning | source |
|---|---|---|
| 1.58–1.67 (61–63% transitions) over all columns; 1.77–1.85 (64–65%) outside soft-masked repeats | Ts/Tv, cannabis NC_044378.1 lastz alignments, three pairs, the haplotype pair at shared cs10 positions (strata kept apart since round 1) | results/ideas_eisenstein_2026-11-10.md:23–26 |
| 61.3–62.5% / 63.9–64.9% | fraction of substitutions that are transitions, all columns / outside soft-masked repeats | results/ideas_eisenstein/kp1_tstv_NC_044378.1.json (all_columns, uppercase_only) |

## Set reconciliation (KP-9, KP-10)

| number | meaning | source |
|---|---|---|
| 0.42 \|A\| | KOMPa vs KOMPb chr7 difference (Supplement) | results/kp9_kompa_kompb_2026-11-10.md (section 1 table) |
| 3,238,101 | canonical 31-mer difference, Cannbio-2 vs Cannbio2_GigaDB (KMC 3.2.4 agrees element for element) | results/kp10_vs_shibuya_2026-11-11.md:36 |
| 7,931 | difference, Cannbio-2 vs its draft (KMC 3.2.4 agrees element for element) | results/kp10_vs_shibuya_2026-11-11.md:37 |
| 11.3 (11.27) | KP-10 bytes per difference, GigaDB pair, 2 rounds: Eisenstein in C, as written (mean 11.26 over 30 OS seeds); both fields in the numpy run | results/kp10_vs_shibuya_2026-11-11.md:104–105; results/kp10_vs_shibuya_2026-11-10_prior.md:60–61 |
| 36.5 MB | KP-10 sketch size, GigaDB pair (36,481,232 B in C; 36,494,096 B in the numpy run) | results/kp10_vs_shibuya_2026-11-11.md:104; results/kp10_vs_shibuya_2026-11-10_prior.md:60 |
| 137 KB, 1 round | KP-10 on the draft pair, exact (137,216 B, 30/30 OS seeds) | results/kp10_vs_shibuya_2026-11-11.md:111 |
| 20.5 | km-peeler k-mer IBLT bytes per difference at 1.05·D, 30/30, whole genome, one round (20.53) | results/kp10_vs_shibuya_2026-11-11.md:108 |
| ~~1.8×~~ RETIRED (round 2) | IBLT / KP-10 bytes (20.53 / 11.27): one of two crossing sketches against one table message stored with 8-byte counters; no longer published | results/kp10_vs_shibuya_2026-11-11.md:10, :121 |
| 15/30 | SUPERSEDED: IBLT successes at n = D on partition 0 (numpy-era run); whole genome now 18/30 | results/kp10_vs_shibuya_2026-11-10_prior.md:90; results/kp10_vs_shibuya_2026-11-11.md:109 |
| 26/30 (4 of 30 failed), D = 488 | SUPERSEDED: IBLT on draft partitions 0–3 at 1.05·D; whole genome (D = 7,931) now 30/30 | results/kp10_vs_shibuya_2026-11-10_prior.md:92; results/kp10_vs_shibuya_2026-11-11.md:112 |
| 4.2 | extended-syncmer IBLT bytes per true k-mer difference | results/kp10_vs_shibuya_2026-11-10_prior.md:102 |
| 2.7× | KP-10 / extended syncmers (11.27 / 4.2) | results/kp10_vs_shibuya_2026-11-10_prior.md:102 |
| 1,986,807 spurious (61% of the true difference D, 38% of the returned set), 7 missed | extended-syncmer recovery, GigaDB | results/kp10_vs_shibuya_2026-11-10_prior.md:105 |
| 48,293,748 | differing extended syncmers, draft vs final | results/kp10_vs_shibuya_2026-11-10_prior.md:110 |
| m + c = 66 | field multiplications per element to BUILD the sketch (recovery adds about 50 per element on GigaDB, round-2 adjudication L7 F7) | results/kp10_vs_shibuya_2026-11-11.md:63–64 |
| 6.9 million (6,913,034) | elements in the mergeability test | results/kp10_vs_shibuya_2026-11-10_prior.md:68 |
| 20,000 | sampled k-mers on which the numpy hashes equal the C, at k = 31 on packed codes only; NOT the k = 63 path, which hashed in C (no longer published, round 2) | results/kp10_vs_shibuya_2026-11-10_prior.md:33 |
| 3 rounds, 257 splits | draft Gaussian run before the Padé index fix (771092b) | results/kp10_vs_shibuya_2026-11-10_prior.md:73–74 |

## Supplement: negative and niche results

| number | meaning | source |
|---|---|---|
| 0.258 | energy in 8 of 31 DFT coefficients at chance; observed 0.257–0.296 | results/ideas_spectral_2026-11-10.md:22, results/ideas_eisenstein_2026-11-10.md (KP-4) |
| 0.210 / 0.766 | KP-4 pruned / prunable at d = 3, seed-sharing candidates | results/ideas_spectral_2026-11-10.md:29 |
| 0.971 / 0.635 → 0.991 / 0.576 | KP-13 precision / recall, quarter → Eisenstein | results/ideas_eisenstein_2026-11-10.md:64 |
| 13,608,457 edges, 0 false / 0 missed | KP-3 adjacency on cs10 31-mers | results/ideas_algebra_2026-11-10.md:19 |
| 1,039 vs 140 ns/node | KP-3 hash vs packed navigation (Gaussian) | results/ideas_algebra_2026-11-10.md:22 |
| 7,760 / 34,736 | KP-5 palindromic 20-mers / inverted repeats (arm 20, spacer ≤ 100) on NC_044378.1 | results/ideas_kp5_kp11_2026-11-11.md:27–28 |
| 1,370× | KP-11 mismatch-counting transform against direct, k = 12, d = 3 | results/ideas_kp5_kp11_2026-11-11.md:74 |

## Consolidated benchmark: SUPERSEDED by the frozen job c84a395ba84b315e (2026-11-12, see Speed); first done as job 08653d9e6a7453a2 (2026-11-11)

Items 1–6 of the earlier requirement (ntHash batched and library; Eisenstein scalar, 1 and 2 lanes; Gaussian 61-bit;
CM31 AVX2 and AVX-512; Eisenstein and CM31 samplers at 1/1000 to 1/4; spaced seeds with one and two masks) ran in
one job, two interleaved repeats, three load baselines per binary. Item 7, the KP-10 sketch build and recovery, is
a numpy prototype and is not timed on laila; the paper quotes its times only as prototype times.

## Round 1 additions (2026-11-11; results/round1_owed_2026-11-11.md)

| number | meaning | source |
|---|---|---|
| 7,931 / 137,216 B / 1 round | sketch of k-mers as written, k = 31, buckets from bits of S, exact | results/round1_owed_2026-11-11.md |
| 398,217,248 | differences faced by a forward-only hash in the numpy prototype of the same sketch (2^8 buckets, m = 64, c = 2), which recovered none and overflowed every bucket | results/round1_owed_2026-11-11.md |
| 21,382 / 548,864 B / 1 round | k = 63, S from the sequence, exact, KMC 3.2.4 agrees element by element | results/round1_owed_2026-11-11.md |
| 0 of 736,547,799 | batched ntHash vs library mismatches, full cs10 | results (nteq), commit 9b1347e |
| 0 of 199,970 | SeedNtHash windows agreeing with their reverse complement, non-symmetric masks | results (seedrc), commit ef2fd8c |
| 45.36% | draft bases on flipped contigs | results/round1_owed_2026-11-11/draft_orientation.json |
| 1,583 gaps / 151 loci / 147 in repeats / 8–96 per seed | long gaps after reseeding x and λ independently | commit 0db2b6e |
| 16,512 | reversal orbits of 8-mers (F·R merges further, 16,362) | commit ec9ec66 |
| 38% / 61% | spurious share of the returned set / of D | results/kp10_vs_shibuya_2026-11-10_prior.md:105 |
| 1,275 | GigaDB buckets SPLIT in round 1, their halves sent in round 2 (numpy run; 1,263 in the C run, so 2,526 halves; mean 1,254 over 30 OS seeds) | results/kp10_vs_shibuya_2026-11-10_prior.md:60; results/kp10_vs_shibuya_2026-11-11.md:104–105 |
| 11.6% | exact 95% upper bound on IBLT failure from 30/30 | adjudication-L2.md (round 1) |
| 9.5% | one-sided 95% Clopper–Pearson bound from 30/30, KP-10 only (30 OS seeds); the IBLT's 30 seeds 1001–1030 share hash functions (seed + 0, 1, 2), so the bound does not hold for it (round 2) | results/kp10_vs_shibuya_2026-11-11.md:116 |
| 21.70 | KP-10 B/diff sized for one round (2^17 buckets), GigaDB, one fixed-seed run | results/kp10_vs_shibuya_2026-11-11.md:107 |
| 0/30, 18/30 | IBLT successes at n = 0.95·D and n = D, GigaDB, whole genome | results/kp10_vs_shibuya_2026-11-11.md:109–110 |
| 17.30 vs 20.54 | KP-10 (1 round, 30/30) vs IBLT at 1.05·D (30/30), draft pair B/diff | results/kp10_vs_shibuya_2026-11-11.md:111–112 |
| 20.76 vs 38.49 | B/diff from D_est = D/8, GigaDB: KP-10 (5 rounds) vs IBLT resend-at-double (4 rounds, 10/10) | results/kp10_vs_shibuya_2026-11-11.md:141 |
| 326 vs about 460 CPU-s | FASTA to recovered k-mers, GigaDB, one thread, laila: KP-10 (kp10c) vs km-peeler (+ canonical sets) | results/kp10_vs_shibuya_2026-11-11.md (section 5) |
| 110.3 vs 3.5 CPU-s | decode: KP-10 vs IBLT diff + list, GigaDB | results/kp10_vs_shibuya_2026-11-11.md (section 5) |
| 1 of 37 | KP-10 runs with an S collision (shared classes, difference still exact), GigaDB | results/kp10_vs_shibuya_2026-11-11.md (section 3, hash collisions) |
| 2.5 pairs | expected colliding pairs allowed by the bound over 442M classes | adjudication-L3.md (round 1) |
| ~~2.2–3.0 GHz~~ SUPERSEDED: 2.1–3.1 GHz in job c84a395ba84b315e (see Speed) | measured core clock in job 08653d9e6a7453a2 (turbo on) | adjudication-L4.md (round 1) |
| 8 vs 3 | multiplications per k-mer, fourth-root vs Eisenstein | adjudication-L1.md F8 (round 1) |
| 35% | haplotype 31-mer difference over the union (42% over one haplotype) | adjudication-L4.md F12 (round 1) |
| 0.6% | max minimizer-density deviation after reseeding (v1gate seed 1, w = 50) | results/kp7_stats_seeds_2026-11-11.gaps.tsv |
| 62 | loci hit by all five methods after reseeding | results/round1_owed_2026-11-11.md |

## KP-10 vs km-peeler, whole genome (paper section "Set sketches", 2026-11-12; results/kp10_vs_shibuya_2026-11-11.md, b082197)

The paper's sketch numbers now come from the C run (laila job b69c5c76857f93ac), not the numpy-era run. Rows above that
cite results/kp10_vs_shibuya_2026-11-10_prior.md for 11.6%, 15/30, 26/30 and the numpy-only hour are superseded here.

| number | meaning | source |
|---|---|---|
| 0.73% (3.2 million), all but 775 in Cannbio-2 | GigaDB difference over \|A\|; \|B∖A\| = 775 | results/kp10_vs_shibuya_2026-11-11.md:36 |
| 11.3, 30/30, 2 rounds | KP-10 B/diff on GigaDB, 30 OS seeds (11.26 mean; 11.27 fixed seed) | results/kp10_vs_shibuya_2026-11-11.md:105 |
| 20.5, 30/30, 1 round | IBLT at 1.05·D, GigaDB | results/kp10_vs_shibuya_2026-11-11.md:108 |
| 21.7 | KP-10 sized for one round (2^17 buckets), one fixed-seed run | results/kp10_vs_shibuya_2026-11-11.md:107 |
| about 1,250 / about 2,500 | buckets split in round 1 (mean 1,254 over 30 seeds) / halves sent in round 2 | results/kp10_vs_shibuya_2026-11-11.md:118 |
| 9.5% | one-sided 95% bound from 30/30, the sketch only; the table's trials share hash functions (round 2) | results/kp10_vs_shibuya_2026-11-11.md:115–117 |
| 20.8 vs 38.5 | B/diff from D_est = D/8, GigaDB: KP-10 per party vs km-peeler's table resent at double size (Results only since round 2; no longer in the abstract) | results/kp10_vs_shibuya_2026-11-11.md:141 |
| 17.3 vs 20.5 | draft pair, one round each, 30/30 each | results/kp10_vs_shibuya_2026-11-11.md:111–112 |
| 0.8m, D/8 to 4D | sizing rule and the D_est sweep | results/kp10_vs_shibuya_2026-11-11.md (section 2, sizing rule; section 4) |
| 3 hash functions | km-peeler default r = 3, source unmodified | results/kp10_vs_shibuya_2026-11-11.md (section 2) |
| 45% | draft bases on flipped contigs (45.36%) | results/round1_owed_2026-11-11/draft_orientation.json |

Carried caveats, all in the paper: no "equal reliability" (30/30 bounds failure only below 9.5%); overflow detection is
probabilistic; each KP-10 party recovers only its own differing k-mers, one IBLT listing gives both sides'; one genome,
two pairs. Since round 2 the paper quotes build, recovery and end-to-end CPU for both methods together (89/110 vs 205/3.5; 326 vs about 460), never one phase alone.

## Round 1's carried open items: what the 2026-11-11/12 results address (reports/adversarial_log.md, "Still open")

| open item | status after this pass |
|---|---|
| as-written k = 31 sketch (7,931 exact, 137 KB) and the forward-only control (398,217,248) | ATTACKED in round 2 (claim C3, NARROWED: the control ran in the numpy prototype; wording "faced"). Before that, REPRODUCED, not attacked: the C run hashed k-mers as written on both pairs, 30/30 OS seeds each, with each side checking only its own data (results/kp10_vs_shibuya_2026-11-11.md). The forward-only control was not rerun. |
| k = 63 sketch (21,382 exact, 549 KB), buckets sized from the true difference | ATTACKED in round 2 (C3): hashed in C, 2^10 buckets at most m/2 each, so one round could not fail; a 0.8m rerun is owed. Before that, OPEN: no new k = 63 run. The D_est sweep (sizing from D/8 to 4D) was at k = 31 only. |
| the reframing (k ≤ 32: canonical code + universal hash is the cheaper provable baseline) | SUPPORTED by measurement: universal hash on every k-mer faster than batched ntHash (results/kp7_kept_eval_2026-11-11.md, re-timed in results/paper_bench_2026-11-12.md). Not attacked. |
| reseeded sampling statistics (1,583 gaps, 151 loci, 147 in repeats) | OPEN: unchanged, not re-examined. |
| unmeasured baselines: canonical code + universal hash row, and gate + ntHash on kept k-mers (critic gaps 2–3) | MEASURED: results/kp7_kept_eval_2026-11-11.md, and in the frozen job (results/paper_bench_2026-11-12.md). |
| km-peeler extended syncmers, 61% spurious vs about 3% in the authors' tests; harness error not ruled out (critic gap 5) | ATTACKED in round 2 (C8): an earlier 5-trial run, the same 7 misses in every trial, the authors' 3% from viral and simulated genomes; rerun owed. Before that, OPEN: the extended-syncmer method was not rerun (results/kp10_vs_shibuya_2026-11-11.md, section 6). The paper still quotes it as measured. |

## Figures (2026-11-12; `manuscript/figures/make_figures.py`, run by `make -C manuscript figures`)

Every plotted value is read by the script from the file below, never typed into it. Figure 3 recomputes its medians
from the raw counters through `paper_bench_table.compute`, the same code that wrote the table it is checked against.

| figure | values | source |
|---|---|---|
| 1A | A = 1, C = ω, G = −ω, T = −1; squared distances 1 / 3 / 4 | paper.md definitions (diagram); results/kp7_eisenstein_2026-11-10.md:35–36 |
| 1B | AACG and CGTT: F and R swap (the script asserts F(rc(s)) = R(s) before drawing) | paper.md definitions (diagram) |
| 2 | 16,384 7-mers, 8,192 classes; F + R 729, F·R 4,160, (F+λ)(R+λ) 8,192 | results/kp7_symmetric_check_2026-11-10.txt:1–4 |
| 3A | batched ntHash 11.37; Eisenstein every k-mer, 2 lanes 24.67; sampler 9.43 / 9.94 / 11.74 / 17.12 at 1/1000 / 1/100 / 1/16 / 1/4; break-even 1/20.2 (slope 28.7, intercept 9.95) | job c84a395ba84b315e, results/paper_bench_2026-11-12/counters.tsv + run.log; agrees with results/paper_bench_2026-11-12.md:30, 33, 37–40, 42 |
| 3B | batched ntHash 11.43; canonical code + universal hash every k-mer 10.24; sampler 9.39 / 9.90 / 11.21 / 15.43; same filter + universal hash 9.15 / 9.52 / 10.25 / 11.36 | job c84a395ba84b315e, same raw files; agrees with results/paper_bench_2026-11-12/table.md (bench_kept rows `lib-e 8 br`, `lib-e 8 bat`, `uh 0`) |
| 4A | protocol (diagram): buckets of S, count and power sums 1..m + c, overflowing bucket split and resent | paper.md Methods, "Set sketches from power sums" |
| 4B | GigaDB pair, D_est/D = 1/8, 1/4, 1/2, 1: sketch 20.76 / 19.40 / 16.69 / 11.27 B/diff in 5 / 4 / 3 / 2 rounds, table 38.49 / 35.93 / 30.79 / 20.53 in 4 / 3 / 2 / 1 rounds | results/kp10_vs_shibuya_2026-11-11/laila_b69c5c76/tables/kp10_dest.tsv:2–4, 6 and iblt_protocol.tsv:6–9; agrees with results/kp10_vs_shibuya_2026-11-11.md:141–144 |

Figure 3's panels are the two benchmark programs and are never compared with each other (per-kept cost 28.7 vs 22.5,
unexplained). No optional spaced-seed figure: the text already says what it would show.

## Round 2 corrections to the set-sketch section (reports/rounds/round-2, 2026-11-12)

| number | meaning | source |
|---|---|---|
| 11.3 per party, 22.5 in all | GigaDB as run: both parties send a sketch and acceptance is joint, so 11.3 is one of two crossing messages (2 × 36,481,232 / 3,238,101) | recon/kp10c.c decode (SA − SB, both checks); round-2 adjudication L7 F1 |
| 17.3, one sketch | draft pair: all 7,931 differences on the draft's side (A∖B = 0), so one sketch suffices | results/kp10_vs_shibuya_2026-11-11/laila_b69c5c76/run.log:20; adjudication L7 F1 |
| 16-byte cells, 8-byte counter; 1.28 cells per difference | km-peeler cell = int64 counter + 8-byte key sum; n = 1.05 D × c_3 = 1.222 → 66,477,002 / 16 / D = 1.283 | km-peeler ibflib.h, ibflib.c:121; adjudications L6 F3, completeness Gap 1 |
| (not published) | packed sketch vs packed table: about equal per message on GigaDB, table smaller on the draft pair; arithmetic only, neither packing built; the paper says only that half the table's bytes are 64-bit counters | round-2 adjudication-completeness.md Gap 1 |
| 3.5% | GigaDB D below the next power-of-two bucket boundary (D / (2^16 · 0.8m) = 0.965) | round-2 adjudications L2 F1, L3 F1 |
| seeds 1001–1030, seed + 0, 1, 2 | km-peeler trial seeds and its derived hash functions; neighbouring trials share two | km-peeler ibflib.c:27–29; recon/kp10_vs_iblt_job.sh; adjudications L1 F4, L4 F3 |
| 1/64 and 4/64 | slices of the GigaDB pair on which the sketch's m and 2^16 (then the 0.8m rule) and the table's 1.05 were chosen | results/kp10_vs_shibuya_2026-11-10_prior.md:47–54 (partition 0), :84–92 (partitions 0–3) |
| 1.84 and 1.87 | each method's cost at D/8 over its own cost at D (20.76/11.27, 38.49/20.53) | tables/kp10_dest.tsv, tables/iblt_protocol.tsv |
| about 20.5 | a doubling IBLT that sends only new cells, at D/8 (arithmetic, not run) | mitzenmacher_2017_multiparty-reconciliation; adjudication L6 F1 |
| 0.05% | set sizes alone estimate D: (\|A\| − \|B\|) / D = 3,236,551 / 3,238,101 on GigaDB; exact on the draft pair | results/kp10_vs_shibuya_2026-11-11.md:36–37 |
| below 0.04 per run | expected silent losses from a differing k-mer colliding with one the other party holds, GigaDB | round-2 adjudication-completeness.md Gap 2 |
| about 50 per element | multiplications recovery adds per element on GigaDB (root flagging) | round-2 adjudication L7 F7 |
| 89 / 110 vs 205 / 3.5 CPU-s; 326 vs about 460 | build per side / decode: KP-10 vs km-peeler text→build / diff + list; FASTA to recovered k-mers | results/kp10_vs_shibuya_2026-11-11.md section 5 |
| 2.4 = 2.17 × 1.12 | universal hash vs Eisenstein, chained through each program's own batched ntHash (24.67/11.37 × 11.43/10.24) | results/paper_bench_2026-11-12.md; adjudication L4 F8 |
| 2^10, mean 20.9 per bucket | k = 63 bucket count, from the true difference (mean load ≤ m/2) | results/round1_owed_2026-11-11.md section 3; adjudications L3 F5, L4 F2 |
| 2^16 | GigaDB first-round buckets (65,536) | results/kp10_vs_shibuya_2026-11-11.md:118 |
| about 3 percent | spurious k-mers relative to the true difference in Shibuya et al.'s tests (viral and simulated genomes) | library papers_md/summaries/shibuya_2022_genomic-reconciliation.md:29 |
| 2,526 | halves sent in round 2 on GigaDB (1,263 splits) | results/kp10_vs_shibuya_2026-11-11.md:104; adjudication L7 F4 |
