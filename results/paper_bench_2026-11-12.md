# The paper's frozen benchmark (2026-11-12)

**One laila job, c84a395ba84b315e, at one commit, 5bee54c**, times every number the manuscript quotes. It replaces job
08653d9e6a7453a2 (results/paper_bench_2026-11-11.md) as the source of every timing in the paper.

- **Host:** laila, Intel Xeon Gold 6330 at 2.0 GHz, kernel 6.1.0-52, gcc 14.2.0 (docker.io/library/gcc:14.2.0),
  `perf_event_paranoid` 2, one of the job's 16 slots (runs are sequential and single-threaded).
- **Code:** kmer-polynomials at **5bee54c03781e73421618055a2d37c46faefe636**, bundled with `laila_perf --git` (git
  archive of that commit), ntHash at c26bd45. Every binary built in the job with `-O3 -march=native`: `bench_e`,
  `bench`, `bench_kept` (hash/Makefile) and `bench_spaced_e` (ideas/eisenstein/Makefile). Same flags on both sides.
- **Runner:** `manuscript/run_paper_bench.sh 5bee54c 5`. **Table:** `python3 manuscript/paper_bench_table.py
  results/paper_bench_2026-11-12` (output `paper_bench_2026-11-12/table.md`). Raw: counters.tsv, run.log, meta.tsv,
  job.json, launcher.log in `paper_bench_2026-11-12/`.
- **Input:** whole cs10 (`cs10.fa`), k = 31, 736,547,799 k-mers.
- **Order:** 5 repeats of 32 entries (a `__load_only__` run per binary and every method, one counted process each);
  odd repeats run the list reversed. Then one more load run per binary, `bench_kept … check`, and `bench_spaced_e`
  (5 interleaved rounds). 165 runs, every exit 0, job state ok, 1,656 s.
- **Per k-mer:** each binary's median load run (6 runs) is subtracted, then divided by 736,547,799. **Values are
  medians of 5.** Load spread: bench_e 1.2%, bench 2.7%, bench_kept 0.8%.
- **Ratios** are against batched ntHash timed by the same binary (bench_e 11.37, bench_kept 11.43; `bench` has none
  and uses bench_e's).
- **Identity:** `test_kp7kept` passed in the build (29,880 sampler configurations); `bench_kept check`: every evaluator
  identical to its reference; every repeat of every bench_kept method kept the same k-mers with the same XOR, and every
  sampler at a rate kept the same count (the table script refuses otherwise).

## Every k-mer, and the samplers as the paper first timed them (bench_e, bench)

| method | cycles/k-mer | spread | vs batched ntHash | job 08653d9e (mean of 2) |
|---|---|---|---|---|
| ntHash batched | 11.37 | 0.4% | 1.00× | 11.37 |
| ntHash library roll() | 23.92 | 1.5% | 0.48× | 23.53 (0.48×) |
| Eisenstein, every k-mer, 1 lane | 27.85 | 1.6% | 0.41× | 27.82 |
| Eisenstein, every k-mer, 2 lanes | 24.67 | 2.8% | 0.46× | 24.69 |
| Gaussian 61-bit (fourth-root), 2 lanes | 47.19 | 0.6% | 0.24× | 47.37 |
| CM31 AVX2, 8 lanes | 26.44 | 0.6% | 0.43× | 26.43 |
| CM31 AVX-512, 16 lanes | 25.94 | 1.3% | 0.44× | 25.86 |
| Eisenstein sampler, branchy, 1/1000 | 9.43 | 1.5% | 1.21× | 9.40 (1.21×) |
| Eisenstein sampler, branchy, 1/100 | 9.94 | 0.7% | 1.14× | 10.01 (1.14×) |
| Eisenstein sampler, batched, 1/16 | 11.74 | 3.8% | 0.97× | 11.74 |
| Eisenstein sampler, batched, 1/4 | 17.12 | 1.5% | 0.66× | 17.08 (0.67×) |

Per kept k-mer (batched, slope 1/16 to 1/4): **28.7 cycles**, intercept 9.95, break-even with batched ntHash at
**1/20.2** (was 28.5 and 1/20.3). Every paper number moved by at most 1.7% (the library roll, 23.53 → 23.92), and no
ratio by more than 0.01.

## Samplers and the round-1 baselines in one harness (bench_kept)

| method | 1/1000 | 1/100 | 1/16 | 1/4 |
|---|---|---|---|---|
| batched ntHash, every k-mer | 11.43 (1.00×) | | | |
| **canonical code + universal hash, every k-mer** (`uh all`) | **10.24 (1.12×)** | | | |
| Eisenstein sampler, branchy (`lib-e 8 br`) | 9.39 (1.22×) | 9.90 (1.15×)\* | — | — |
| Eisenstein sampler, batched (`lib-e 8 bat`) | 9.93 (1.15×) | 10.19 (1.12×) | 11.21 (1.02×) | 15.43 (0.74×) |
| `ec 13`, branchy | 9.15 (1.25×) | 9.65 (1.18×) | — | — |
| gate + ntHash on kept, branchy / batched (`nt 8`) | 9.72 (1.18×) | 10.23 (1.12×) | 11.08 (1.03×) | 14.70 (0.78×) |
| gate + universal hash on kept, branchy / batched (`uh 0`) | 9.15 (1.25×) | 9.52 (1.20×) | 10.25 (1.12×) | 11.36 (1.01×) |

\* One of 5 runs was slow (spread 43%, minimum 9.88); the median is unaffected.

Per kept k-mer (batched, slope 1/16 to 1/4): **Eisenstein 22.5** (break-even 1/13.8), **gate + ntHash 19.3**,
**gate + universal hash 5.9**.

Reading:
- **The universal hash on every k-mer is 1.12× batched ntHash** (was 1.13× in job 88c9cc63e42a63a9). For k ≤ 32 it is
  the cheaper provable strand-symmetric hash, which supports the paper's reframing.
- **At sketch rates the sampler's gain is the packed-code filter, not the hash.** The same filter followed by the
  universal hash is as fast at 1/1000 (9.15 against 9.39) and faster at 1/100 (9.52 against 9.90).
- **`ec 13` branchy** reproduces its gain (1.25× / 1.18× against 1.22× / 1.15×), the same speed as gate + universal hash
  at 1/1000. The kept-eval note found it to be gate-loop code generation, not evaluation. The paper does not quote it.

## Two harnesses, one per-kept cost range

The two binaries call the same library samplers (`kp7e_sample_b`, `kp7e_sample_batch`) with different harnesses. At
sketch rates they agree (1/1000: 9.43 and 9.39). At 1/4 they differ: 17.12 in bench_e and 15.43 in bench_kept, with
70.1 and 68.9 instructions per k-mer. So the per-kept slope is 28.7 in bench_e and 22.5 in bench_kept, in the same job.
The cause is not diagnosed. The paper quotes the range (about 22–29 cycles per kept k-mer), compares against the
baselines only within bench_kept, and quotes the break-even (1/20) from bench_e with its own ntHash.

## Spaced seeds (run_165, bench_spaced_e's own wall-clock medians of 5 rounds, first 60,000,000 bases of cs10)

| masks | Gaussian | Eisenstein | SeedNtHash | Eisenstein / SeedNtHash, per hash |
|---|---|---|---|---|
| 1 | 10.7 M/s | 21.8 M/s | 29.2 M/s | **0.75×** (0.76× wall-clock) |
| 2 | 9.8 M/s | 20.3 M/s | 28.3 M/s | 0.72× (0.73× wall-clock) |

Windows unchanged at 1% substitutions, non-palindromic mask: ours 0.74722 (the exact double-stranded fraction),
SeedNtHash 0.81794. Unchanged from job 08653d9e6a7453a2.

## Not in this job

The KP-10 sketch build and recovery (laila job b69c5c76857f93ac, results/kp10_vs_shibuya_2026-11-11.md), and the
Supplement's timings (earlier jobs, as the paper says).
