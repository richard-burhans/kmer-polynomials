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
| 8,192 / 32,896 | the same for the Eisenstein hash, with F injective on all 4^k k-mers | results/kp7_eisenstein_2026-11-10.md:33–34 |
| 2(k−1)/q ≈ 2.6 × 10⁻¹⁷ | collision bound of the Eisenstein symmetric hash per pair, k = 31, q = 2⁶¹ − 1 (corrected in round 1; was 2.7, which is 2k/q) | results/kp7_eisenstein_2026-11-10.md:29–30 |
| 2(k−1)/p² ≈ 1.3 × 10⁻¹⁷ | collision bound of the CM31 (x complex) hash per pair, k = 31 (corrected in round 1; was 1.4) | results/kp7_m31_simd_2026-11-10.md (Method paragraph) |
| 1 / 3 / 4 | squared distances in ℂ: transitions, A–C/G–T transversions, complements (Eisenstein) | results/kp7_eisenstein_2026-11-10.md:35–36 |
| 3 multiplies per k-mer | Eisenstein rolling cost | results/kp7_eisenstein_2026-11-10.md:19–20 |
| 16 KB | Eisenstein v2 table size (8-byte entries) | results/kp7_eisenstein_2026-11-10.md:22 |
| θ = π/3 | Das, Nanda & Das 2017 fix θ explicitly | literature/eisenstein_wantlist.md (Resolution section) |
| eq. (5), eq. (4) | Cristea 2002: the fourth-root map (up to sign) and the θ = π/4 point | library papers_md/summaries/cristea_2002_genomic-signals.md:25–30; literature/eisenstein_wantlist.md:103–104 |
| rationale restates Cristea's | Das 2017's purine/pyrimidine + complements-sum-to-zero rationale, no distances, θ asserted | library papers_md/summaries/das_2017_trigonometric-mapping.md:23–27, 63–64 |
| sensitivity = specificity formula; 3 exon ranges | Das 2017 reliability caveats | library papers_md/summaries/das_2017_trigonometric-mapping.md:34–39, 48–51 |
| unread | Das, Nanda & Das 2019 Genomics, the mapping's source cited by Gulocak 2022 | literature/eisenstein_wantlist.md:93–96 |
| binary codes with transition-aware distance | Yu, Li & Yu 2018 (narrows the novelty claim) | library papers_md/summaries/yu_2018_encoding-survey.md:24–31 |

## Speed (laila, Xeon Gold 6330 at 2.0 GHz, cycles per input k-mer, cs10, k = 31)

Input: 736,547,799 k-mers (results/paper_bench_2026-11-11.md:8). Every timing below comes from ONE job, 08653d9e6a7453a2, run by
`manuscript/run_paper_bench.sh`; the table is `manuscript/paper_bench_table.py` on its raw counters.

| number | meaning | job | source |
|---|---|---|---|
| 11.37 | batched ntHash | 08653d9e6a7453a2 | results/paper_bench_2026-11-11.md:22 |
| 23.53 (0.48×) | ntHash library roll() | 08653d9e6a7453a2 | results/paper_bench_2026-11-11.md:23 |
| 24.69 (0.46×) | Eisenstein, every k-mer, 2 lanes | 08653d9e6a7453a2 | results/paper_bench_2026-11-11.md:26 |
| 27.82 (0.41×) | Eisenstein, every k-mer, 1 lane | 08653d9e6a7453a2 | results/paper_bench_2026-11-11.md:25 |
| 47.37 (0.24×) | Gaussian 61-bit, every k-mer, 2 lanes | 08653d9e6a7453a2 | results/paper_bench_2026-11-11.md:28 |
| 26.43 (0.43×) | CM31 AVX2, every k-mer, 8 lanes | 08653d9e6a7453a2 | results/paper_bench_2026-11-11.md:42 |
| 25.86 (0.44×) | CM31 AVX-512, every k-mer, 16 lanes | 08653d9e6a7453a2 | results/paper_bench_2026-11-11.md:43 |
| 9.40 (1.21×) | Eisenstein sampler, branchy, 1/1000 | 08653d9e6a7453a2 | results/paper_bench_2026-11-11.md:30 |
| 10.01 (1.14×) | Eisenstein sampler, branchy, 1/100 | 08653d9e6a7453a2 | results/paper_bench_2026-11-11.md:31 |
| 11.74 (0.97×) | Eisenstein sampler, batched, 1/16 | 08653d9e6a7453a2 | results/paper_bench_2026-11-11.md:36 |
| 17.08 (0.67×) | Eisenstein sampler, batched, 1/4 | 08653d9e6a7453a2 | results/paper_bench_2026-11-11.md:37 |
| ≈ 28.5 | Eisenstein cycles per kept k-mer (slope 1/16 to 1/4) | 08653d9e6a7453a2 | results/paper_bench_2026-11-11.md:47 |
| ≈ 1/20 (1/20.3) | Eisenstein sampler break-even against batched ntHash | 08653d9e6a7453a2 | results/paper_bench_2026-11-11.md:47 |
| 2.07× | batched ntHash against the library roll() (23.53 / 11.37), "about twice" | 08653d9e6a7453a2 | results/paper_bench_2026-11-11.md:22–23 |
| 64.3 | CM31 v2 cycles per kept k-mer (batched), previous best | 08653d9e6a7453a2 | results/paper_bench_2026-11-11.md:49 |
| ≤ 1.3% | spread between the two repeats, every row in the paper's table | 08653d9e6a7453a2 | results/paper_bench_2026-11-11.md:12 |

Superseded (earlier drafts, other jobs): 11.41 / 23.44 / 24.59 / 27.95 / 47.51 / 25.83 / 9.41 / 9.97 / 11.87 / 17.21
from 7fc2a700073c1f4f, and the break-even "≈ 1/19", which was an arithmetic slip (intercept 9.9 for 10.09; see
results/paper_bench_2026-11-11.md, per-kept section).

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

## Spaced seeds (consolidated job 08653d9e6a7453a2)

| number | meaning | source |
|---|---|---|
| 0.75× per hash (0.76× wall-clock; single wall-clock run) | Eisenstein spaced-seed hash against ntHash2 SeedNtHash, one mask, first 60 Mb of cs10 (was 0.78× in job af1eb52f65c6d471) | results/paper_bench_2026-11-11.md:58 |
| 0.747 (0.74722) | fraction of windows unchanged at 1% SNPs, non-palindromic mask: exact double-stranded value, matched by ours | results/paper_bench_2026-11-11.md:61 |
| 0.818 (0.81794) | SeedNtHash's value on the same test | results/paper_bench_2026-11-11.md:62 |

## Transition frequency (KP-1)

| number | meaning | source |
|---|---|---|
| 1.58–1.67 (61–63% transitions); 1.85 (65%) between the two haplotypes outside repeats | Ts/Tv, cannabis NC_044378.1 alignments (strata quoted separately after round 1) | results/ideas_eisenstein_2026-11-10.md:23–26 |
| 61.3–64.9% | fraction of substitutions that are transitions | results/ideas_eisenstein_2026-11-10.md:23–26 |

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
| 1.8× | IBLT / KP-10 bytes (20.53 / 11.27): a two-round sketch against a one-round table | results/kp10_vs_shibuya_2026-11-11.md:10, :121 |
| 15/30 | SUPERSEDED: IBLT successes at n = D on partition 0 (numpy-era run); whole genome now 18/30 | results/kp10_vs_shibuya_2026-11-10_prior.md:90; results/kp10_vs_shibuya_2026-11-11.md:109 |
| 26/30 (4 of 30 failed), D = 488 | SUPERSEDED: IBLT on draft partitions 0–3 at 1.05·D; whole genome (D = 7,931) now 30/30 | results/kp10_vs_shibuya_2026-11-10_prior.md:92; results/kp10_vs_shibuya_2026-11-11.md:112 |
| 4.2 | extended-syncmer IBLT bytes per true k-mer difference | results/kp10_vs_shibuya_2026-11-10_prior.md:102 |
| 2.7× | KP-10 / extended syncmers (11.27 / 4.2) | results/kp10_vs_shibuya_2026-11-10_prior.md:102 |
| 1,986,807 spurious (61% of the true difference D, 38% of the returned set), 7 missed | extended-syncmer recovery, GigaDB | results/kp10_vs_shibuya_2026-11-10_prior.md:105 |
| 48,293,748 | differing extended syncmers, draft vs final | results/kp10_vs_shibuya_2026-11-10_prior.md:110 |
| m + c = 66 | field multiplications per element | results/kp10_vs_shibuya_2026-11-10_prior.md:8 |
| 6.9 million (6,913,034) | elements in the mergeability test | results/kp10_vs_shibuya_2026-11-10_prior.md:68 |
| 20,000 | sampled k-mers on which the numpy hashes equal the C | results/kp10_vs_shibuya_2026-11-10_prior.md:33 |
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

## Consolidated benchmark: DONE (job 08653d9e6a7453a2, 2026-11-11)

Items 1–6 of the earlier requirement (ntHash batched and library; Eisenstein scalar, 1 and 2 lanes; Gaussian 61-bit;
CM31 AVX2 and AVX-512; Eisenstein and CM31 samplers at 1/1000 to 1/4; spaced seeds with one and two masks) ran in
one job, two interleaved repeats, three load baselines per binary. Item 7, the KP-10 sketch build and recovery, is
a numpy prototype and is not timed on laila; the paper quotes its times only as prototype times.

## Round 1 additions (2026-11-11; results/round1_owed_2026-11-11.md)

| number | meaning | source |
|---|---|---|
| 7,931 / 137,216 B / 1 round | sketch of k-mers as written, k = 31, buckets from bits of S, exact | results/round1_owed_2026-11-11.md |
| 398,217,248 | differences seen by a forward-only hash in the same sketch | results/round1_owed_2026-11-11.md |
| 21,382 / 548,864 B / 1 round | k = 63, S from the sequence, exact, KMC 3.2.4 agrees element by element | results/round1_owed_2026-11-11.md |
| 0 of 736,547,799 | batched ntHash vs library mismatches, full cs10 | results (nteq), commit 9b1347e |
| 0 of 199,970 | SeedNtHash windows agreeing with their reverse complement, non-symmetric masks | results (seedrc), commit ef2fd8c |
| 45.36% | draft bases on flipped contigs | results/round1_owed_2026-11-11/draft_orientation.json |
| 1,583 gaps / 151 loci / 147 in repeats / 8–96 per seed | long gaps after reseeding x and λ independently | commit 0db2b6e |
| 16,512 | reversal orbits of 8-mers (F·R merges further, 16,362) | commit ec9ec66 |
| 38% / 61% | spurious share of the returned set / of D | results/kp10_vs_shibuya_2026-11-10_prior.md:105 |
| 1,275 | GigaDB buckets resent in round 2 (numpy run; 1,263 in the C run, mean 1,254 over 30 OS seeds) | results/kp10_vs_shibuya_2026-11-10_prior.md:60; results/kp10_vs_shibuya_2026-11-11.md:104–105 |
| 11.6% | exact 95% upper bound on IBLT failure from 30/30 | adjudication-L2.md (round 1) |
| 9.5% | one-sided 95% Clopper–Pearson bound from 30/30, both KP-10 (30 OS seeds) and the IBLT (30 seeds), whole genome | results/kp10_vs_shibuya_2026-11-11.md:116 |
| 21.70 | KP-10 B/diff sized for one round (2^17 buckets), GigaDB, one fixed-seed run | results/kp10_vs_shibuya_2026-11-11.md:107 |
| 0/30, 18/30 | IBLT successes at n = 0.95·D and n = D, GigaDB, whole genome | results/kp10_vs_shibuya_2026-11-11.md:109–110 |
| 17.30 vs 20.54 | KP-10 (1 round, 30/30) vs IBLT at 1.05·D (30/30), draft pair B/diff | results/kp10_vs_shibuya_2026-11-11.md:111–112 |
| 20.76 vs 38.49 | B/diff from D_est = D/8, GigaDB: KP-10 (5 rounds) vs IBLT resend-at-double (4 rounds, 10/10) | results/kp10_vs_shibuya_2026-11-11.md:141 |
| 326 vs about 460 CPU-s | FASTA to recovered k-mers, GigaDB, one thread, laila: KP-10 (kp10c) vs km-peeler (+ canonical sets) | results/kp10_vs_shibuya_2026-11-11.md (section 5) |
| 110.3 vs 3.5 CPU-s | decode: KP-10 vs IBLT diff + list, GigaDB | results/kp10_vs_shibuya_2026-11-11.md (section 5) |
| 1 of 37 | KP-10 runs with an S collision (shared classes, difference still exact), GigaDB | results/kp10_vs_shibuya_2026-11-11.md (section 3, hash collisions) |
| 2.5 pairs | expected colliding pairs allowed by the bound over 442M classes | adjudication-L3.md (round 1) |
| 2.2–3.0 GHz | measured core clock in job 08653d9e6a7453a2 (turbo on) | adjudication-L4.md (round 1) |
| 8 vs 3 | multiplications per k-mer, fourth-root vs Eisenstein | adjudication-L1.md F8 (round 1) |
| 35% | haplotype 31-mer difference over the union (42% over one haplotype) | adjudication-L4.md F12 (round 1) |
| 0.6% | max minimizer-density deviation after reseeding (v1gate seed 1, w = 50) | results/kp7_stats_seeds_2026-11-11.gaps.tsv |
| 62 | loci hit by all five methods after reseeding | results/round1_owed_2026-11-11.md |
