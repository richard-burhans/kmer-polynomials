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
| ≈ 2k/q ≈ 2.7 × 10⁻¹⁷ | collision bound of the Eisenstein symmetric hash per pair, k = 31, q = 2⁶¹ − 1 | results/kp7_eisenstein_2026-11-10.md:29–30 |
| ≈ 1.4 × 10⁻¹⁷ | collision bound of the CM31 (x complex) hash per pair, k = 31 | results/kp7_m31_simd_2026-11-10.md (Method paragraph) |
| 1 / 3 / 4 | squared distances in ℂ: transitions, A–C/G–T transversions, complements (Eisenstein) | results/kp7_eisenstein_2026-11-10.md:35–36 |
| 3 multiplies per k-mer | Eisenstein rolling cost | results/kp7_eisenstein_2026-11-10.md:19–20 |
| 16 KB | Eisenstein v2 table size (8-byte entries) | results/kp7_eisenstein_2026-11-10.md:22 |
| θ = π/3 | Das, Nanda & Das 2017 fix θ explicitly | literature/eisenstein_wantlist.md (Resolution section) |

## Speed (laila, Xeon Gold 6330 at 2.0 GHz, cycles per input k-mer, cs10, k = 31)

Input: 736,547,799 k-mers (results/kp7_eisenstein_2026-11-10.md:58).

| number | meaning | job | source |
|---|---|---|---|
| 11.41 | batched ntHash | 7fc2a700073c1f4f | results/kp7_eisenstein_2026-11-10.md:63 |
| 23.44 | ntHash library roll() | 7fc2a700073c1f4f | results/kp7_eisenstein_2026-11-10.md:64 |
| 24.59 | Eisenstein, every k-mer, 2 lanes (0.46×) | 7fc2a700073c1f4f | results/kp7_eisenstein_2026-11-10.md:68 |
| 27.95 | Eisenstein, every k-mer, 1 lane | 7fc2a700073c1f4f | results/kp7_eisenstein_2026-11-10.md:67 |
| 47.51 | Gaussian 61-bit, every k-mer, 2 lanes | 7fc2a700073c1f4f | results/kp7_eisenstein_2026-11-10.md:70 |
| 25.83 | CM31 AVX-512, every k-mer, 16 lanes | 7fc2a700073c1f4f | results/kp7_eisenstein_2026-11-10.md:71 |
| 9.41 (1.21×) | Eisenstein sampler, branchy, 1/1000 | 7fc2a700073c1f4f | results/kp7_eisenstein_2026-11-10.md:73 |
| 9.97 (1.14×) | Eisenstein sampler, branchy, 1/100 | 7fc2a700073c1f4f | results/kp7_eisenstein_2026-11-10.md:74 |
| 11.87 (0.96×) | Eisenstein sampler, batched, 1/16 | 7fc2a700073c1f4f | results/kp7_eisenstein_2026-11-10.md:79 |
| 17.21 (0.66×) | Eisenstein sampler, batched, 1/4 | 7fc2a700073c1f4f | results/kp7_eisenstein_2026-11-10.md:80 |
| ≈ 28.5 | Eisenstein cycles per kept k-mer (slope 1/16 to 1/4) | 7fc2a700073c1f4f | results/kp7_eisenstein_2026-11-10.md:91 |
| ≈ 1/19 | Eisenstein sampler break-even against batched ntHash | 7fc2a700073c1f4f | results/kp7_eisenstein_2026-11-10.md:96 |
| 2.09× | batched ntHash against the library roll() | 5a03e1322aef25af | results/kp7_next_steps_2026-11-10.md:18 |
| 64.4 | CM31 v2 cycles per kept k-mer (batched) | 5deca06fa147bf0e | results/kp7_kept_eval_2026-11-10.md:43 |

## Sampling quality (local, full cs10, 5 seeds per method, no job)

| number | meaning | source |
|---|---|---|
| 0.14 | geometric expectation of gaps > 200 k-mers per seed, rate 2/21 | results/kp7_next_steps_2026-11-10.md:50–51 |
| about 1,450 long gaps, all but 2 in soft-masked repeats | across 4 methods × 5 seeds | results/kp7_next_steps_2026-11-10.md:64–65 |
| 0.00088 to 0.00103 | mean max per-bit deviation from ½ (ntHash, v2 gate, v1 gate, Gaussian provable) | results/kp7_next_steps_2026-11-10.md:55–58 |
| 47.7% | soft-masked fraction of cs10 | results/kp7_next_steps_2026-11-10.md:50 |

Note: these statistics were measured with the Gaussian provable hash and the packed-code gate. The Eisenstein sampler
uses the same gate (results/kp7_eisenstein_2026-11-10.md:46–47), but the Eisenstein provable hash itself has no
multi-seed statistics yet.

## Spaced seeds (KP-14, laila job af1eb52f65c6d471)

| number | meaning | source |
|---|---|---|
| 0.78× | Eisenstein spaced-seed hash against ntHash2 SeedNtHash, one mask | results/ideas_eisenstein_2026-11-10.md:105 |
| 0.74724 | fraction of windows unchanged at 1% SNPs, non-palindromic mask: exact double-stranded value, matched by ours | results/ideas_algebra_2026-11-10.md:55 |
| 0.818 | SeedNtHash's value on the same test | results/ideas_eisenstein_2026-11-10.md (KP-14 section) |

## Transition frequency (KP-1)

| number | meaning | source |
|---|---|---|
| 1.58–1.85 | Ts/Tv, cannabis NC_044378.1 pairwise columns | results/ideas_eisenstein_2026-11-10.md:23–26 |
| 61.3–64.9% | fraction of substitutions that are transitions | results/ideas_eisenstein_2026-11-10.md:23–26 |

## Set reconciliation (KP-9; KP-10 pending)

| number | meaning | source |
|---|---|---|
| 0.42 \|A\| | KOMPa vs KOMPb chr7 difference | results/kp9_kompa_kompb_2026-11-10.md (section 1 table) |
| 2,180× / 273× | sketch smaller than the set, 100 / 1,000 SNPs | results/kp9_kompa_kompb_2026-11-10.md (section 2 table) |
| 17 bytes, 25–30 measured | bytes per difference (design load / with power-of-two rounding) | results/kp9_kompa_kompb_2026-11-10.md section 2 |
| pending | KP-10 against Shibuya's km-peeler | recon/ (running) |

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

## Consolidated benchmark required before submission

The timings above come from at least four laila jobs (7fc2a700073c1f4f, 5a03e1322aef25af, 5deca06fa147bf0e,
af1eb52f65c6d471). Before submission, ONE laila job must time, in interleaved repeats with a `__load_only__` baseline:
1. batched ntHash (`ntbatch.h`) and the ntHash library roll()
2. Eisenstein every-k-mer: scalar reference, 1 lane, 2 lanes
3. Gaussian 61-bit every-k-mer, 2 lanes, and CM31 AVX2 and AVX-512
4. Eisenstein sampler, branchy and batched, at 1/1000, 1/100, 1/16, 1/4
5. CM31 v2 sampler (batched) at the same rates, as the previous best
6. the Eisenstein spaced-seed hash and ntHash2 SeedNtHash (one and two masks)
7. the KP-10 sketch build and recovery, once that code is final

Also needed: multi-seed sampling statistics for the Eisenstein provable hash itself.
