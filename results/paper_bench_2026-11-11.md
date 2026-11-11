# The paper's consolidated benchmark (2026-11-11)

**One laila job, 08653d9e6a7453a2**, times every method the manuscript reports, so no timing in the paper is compared
across jobs. Laila: Intel Xeon Gold 6330 at 2.0 GHz, kernel 6.1.0-52, gcc 14.2.0 in the job, `-O3 -march=native`,
single-threaded, user-space counters via `cannabis-genome/scripts/laila_perf.py run` (UDT laila-perf).

- **Runner:** `manuscript/run_paper_bench.sh` (f56af99). kmer-polynomials at f56af99, ntHash at c26bd45. Input: whole
  cs10 (`cs10.fa`, 221 records, 876,147,649 bases, 736,547,799 31-mers).
- **Order:** for each of 2 repeats: a `__load_only__` run of `bench_e`, the 20 `bench_e` methods, a `__load_only__`
  run of `bench`, the 2 `bench` SIMD methods. Then one more load run per binary, then `bench_spaced_e` once. 51 runs,
  every exit 0.
- **Per k-mer:** each binary's mean load baseline is subtracted (bench_e 6.98e9 cycles, spread 0.7% over 3 runs;
  bench 7.65e9, 0.1%), then divided by 736,547,799. Values are the mean of the 2 repeats; "spread" is
  (max − min) / mean. Every spread is at most 1.3%, except the scalar reference roller (4.0%).
- **Table:** `python3 manuscript/paper_bench_table.py <results dir>`, output `paper_bench_2026-11-11/table.md`.
  Raw: `paper_bench_2026-11-11/` (counters.tsv, run.log, meta.tsv, job.json, launcher.log).

## Cycles per input k-mer

| method | binary | cycles/k-mer | spread | instructions/k-mer | branch-miss rate | vs batched ntHash |
|---|---|---|---|---|---|---|
| ntHash batched (ntbatch) | bench_e | 11.37 | 0.7% | 47.0 | 0.34% | 1.00× |
| ntHash library roll | bench_e | 23.53 | 0.0% | 110.1 | 0.23% | 0.48× |
| e scalar roller (reference) | bench_e | 77.30 | 4.0% | 172.0 | 3.83% | 0.15× |
| e sym_seq L1 | bench_e | 27.82 | 0.3% | 72.6 | 0.32% | 0.41× |
| e sym_seq L2 | bench_e | 24.69 | 0.4% | 78.5 | 0.33% | 0.46× |
| gauss61 sym_seq L1 | bench_e | 50.00 | 0.0% | 160.2 | 0.31% | 0.23× |
| gauss61 sym_seq L2 | bench_e | 47.37 | 0.2% | 179.3 | 0.33% | 0.24× |
| m31 simd (avx512 else avx2) | bench_e | 25.99 | 0.6% | 60.2 | 0.37% | 0.44× |
| e branchy 1/1000 | bench_e | 9.40 | 0.2% | 43.0 | 0.35% | 1.21× |
| e branchy 1/100 | bench_e | 10.01 | 0.4% | 43.8 | 0.42% | 1.14× |
| e branchy 1/16 | bench_e | 13.32 | 0.7% | 48.8 | 0.86% | 0.85× |
| e branchy 1/4 | bench_e | 25.39 | 0.6% | 66.6 | 2.35% | 0.45× |
| e batched 1/1000 | bench_e | 9.95 | 0.5% | 46.0 | 0.34% | 1.14× |
| e batched 1/100 | bench_e | 10.22 | 0.3% | 46.9 | 0.34% | 1.11× |
| e batched 1/16 | bench_e | 11.74 | 0.8% | 51.9 | 0.34% | 0.97× |
| e batched 1/4 | bench_e | 17.08 | 0.3% | 70.1 | 0.33% | 0.67× |
| m31 branchy 1/1000 | bench_e | 9.52 | 0.5% | 43.1 | 0.35% | 1.19× |
| m31 branchy 1/100 | bench_e | 10.44 | 1.3% | 44.8 | 0.42% | 1.09× |
| m31 batched 1/16 | bench_e | 13.67 | 1.0% | 58.0 | 0.33% | 0.83× |
| m31 batched 1/4 | bench_e | 25.72 | 0.6% | 94.3 | 0.34% | 0.44× |
| m31 avx2 8-lane | bench | 26.43 | 0.3% | 61.7 | 0.32% | 0.43× |
| m31 avx512 16-lane | bench | 25.86 | 0.6% | 60.2 | 0.34% | 0.44× |

## Per-kept cost and break-even (batched samplers, slope from 1/16 to 1/4)

- **Eisenstein:** (17.08 − 11.74) / 0.1875 = **28.5 cycles per kept k-mer**, intercept 9.96. It beats batched ntHash
  (11.37) below a sampling rate of **1/20.3**, and the 2-lane every-k-mer hash (24.69) below 1/1.9.
- **CM31 v2 (previous best):** 64.3 cycles per kept k-mer, intercept 9.66, beats batched ntHash below 1/37.6.
- **Correction to earlier drafts.** results/kp7_eisenstein_2026-11-10.md gave the Eisenstein break-even as "about
  1/19" from an intercept of 9.9. Its own numbers give 11.87 − 28.5 × 0.0625 = 10.09 and a break-even of 1/21.5;
  the 9.9 was an arithmetic slip. This job gives 1/20.3.

## Spaced seeds (run_51; bench_spaced_e's own wall-clock rates, first 60,000,000 bases of cs10, 3 rounds)

| masks | Gaussian | Eisenstein | SeedNtHash | Eisenstein / SeedNtHash |
|---|---|---|---|---|
| 1 | 11.1 M/s | 21.5 M/s | 28.8 M/s | **0.76×** |
| 2 | 10.2 M/s | 20.1 M/s | 28.5 M/s | **0.71×** |

Windows unchanged under substitutions, non-palindromic mask, at 0.001 / 0.01 / 0.05: ours 0.97146 / 0.74722 /
0.22567, equal to the exact double-stranded fraction; SeedNtHash 0.98023 / 0.81794 / 0.35839. (On
cs10_head60M.fa the earlier run gave 0.74724; the input here is the first 60 Mb of the full FASTA.)

## Against the numbers the paper had (job 7fc2a700073c1f4f and the KP-14 job)

All within the 3–6% between-job variance. Paper-relevant moves: batched ntHash 11.41 → 11.37, library 23.44 →
23.53 (0.49 → 0.48×), e L2 24.59 → 24.69, e L1 27.95 → 27.82, gauss61 L2 47.51 → 47.37, AVX-512 25.83 → 25.86 (from
`bench`; `bench_e`'s copy gave 25.99), sampler 1/1000 9.41 → 9.40, 1/100 9.97 → 10.01, 1/16 11.87 → 11.74, 1/4
17.21 → 17.08 (0.66 → 0.67×), spaced 0.78× → 0.76×, break-even 1/19 → 1/20.3. New: AVX2 8-lane 26.43.

Not in this job: the KP-10 sketch build and recovery (a numpy prototype; its times are prototype times, in
results/kp10_vs_shibuya_2026-11-11.md).
