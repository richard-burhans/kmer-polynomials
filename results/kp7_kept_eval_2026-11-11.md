# KP-7 kept-k-mer evaluation, redone: what is left to cut (2026-11-11)

**Verdict: none of the four ideas makes a kept k-mer cheaper than the current Eisenstein 16 KB table evaluator.**

- **The current evaluator** (`kp7e_sample_batch`): 21.8 cycles per kept k-mer in the confirm job and 25.3 in the
  sweep (28.5 in the paper's job).
- **Chunk widths:** 11-, 13- and 16-bit chunks come in at 23.8–30.3, the same or worse.
- **Batch sizes:** 64 and 4,096 give 27.0, against 21.8 at the default 512.
- **The bit-plane, real-x evaluators:** 45–65 cycles.
- **The CM31 16 KB tables:** 55.7–59.9.

The per-kept cost is now about 22–28 cycles of lookups, additions and one multiply. At sketch rates it hardly
matters: 1/100 × 25 = 0.25 cycles per input k-mer. The sampler's cost there is the gate loop, about 9 cycles per
input k-mer.

**Best variant against batched ntHash** (confirm job 88c9cc63e42a63a9, median of 5, batched ntHash 11.19 cycles per
k-mer in the same job):

| rate | best variant | cycles/k-mer | vs batched ntHash | current sampler (`kp7e_sample_b`) |
|---|---|---|---|---|
| 1/1000 | `ec 13`, branchy | 8.91 | **1.26×** | 9.14 (1.22×) |
| 1/100 | `ec 13`, branchy | 9.44 | **1.18×** | 10.41 (1.08×) |

The 1/1000 gain comes from the gate loop, not the evaluator. The kp7kept branchy loop retires 40.0 instructions per
k-mer against kp7e's 43.0, with only 0.1% of k-mers kept. The rate-0 (nothing kept) runs carry it: `lib-e 8 br 0`
measured 9.09.

**The round-1 baselines are now measured** (critic gaps 2–3, reports/adversarial_log.md):
- the canonical code plus a universal hash on EVERY k-mer: **9.91 cycles, 1.13× batched ntHash**
- the same packed-code gate followed by ntHash on the kept k-mers: 20.2 cycles per kept k-mer
- the gate followed by a multiply-shift universal hash on the kept k-mers: 5.9 cycles per kept k-mer

So for k ≤ 32 the universal hash is cheaper than KP-7 on kept k-mers, and also on every k-mer. The case for KP-7
stays its algebra, as the round-1 reframing already says.

## Method

**Code** (all at kmer-polynomials 8ea3e51 for the confirm job and da4c32a for the sweep):
- `hash/kp7kept.{h,c}`: the evaluators
- `hash/test_kp7kept.cpp`: the tests
- `hash/bench_kept.cpp`: the bench
- `hash/run_kept_laila.sh`: the runner
- `hash/kept_table.py`: the table script

**Measurement:**
- Laila: Xeon Gold 6330, kernel 6.1.0-52, gcc 14.2.0 in the job.
- Every binary was built in the job by `hash/Makefile` with `OPT = -O3 -march=native`. That covers the KP-7 code,
  ntbatch and both baselines alike, and all of them were compiled into the same `bench_kept` binary.
- Counters come from cannabis-genome `scripts/laila_perf.py run`.
- The bundle is `git archive` of the commit (`--git`), so the timed code is exactly that commit.
- Input: whole cs10, k = 31, 736,547,799 k-mers.
- One counted process per method (`bench_kept $DATA 31 "only=<method>"`).
- Repeats are interleaved, and odd repeats run the list reversed.
- The job's mean `__load_only__` run is subtracted, then the result divided by the k-mer count. Each value is the
  median over repeats.
- Every job also ran:
  - `test_kp7kept`: 29,880 sampler configurations against direct references, k = 1–32
  - `make mutation_kept`: both mutants killed
  - `bench_kept … check`: every evaluator identical to its family's existing sampler on whole cs10, at 1/16 and
    1/1000, batched and branchy

The identity holds:
- `kept_table.py` refuses to print unless the check passed, every repeat of a method gives the same count and XOR,
  every sampler at a rate keeps the same k-mers, and each hash family gives one XOR per rate.
- Kept counts: 745,710 / 7,382,545 / 45,892,012 / 183,856,486 at 1/1000 / 1/100 / 1/16 / 1/4.

**Commands:**

```
hash/run_kept_laila.sh sweep   da4c32a 3     # job 3dddf48cdf81c0ca, 65 methods x 3, 200 runs
hash/run_kept_laila.sh confirm 8ea3e51 5     # job 88c9cc63e42a63a9, 32 methods x 5, 167 runs
python3 hash/kept_table.py results/kp7_kept_eval_2026-11-11/<job> [--min]
```

Raw outputs (counters.tsv, run.log, meta.tsv, job.json, launcher.log) and the tables (`table.md` with medians,
`table_min.md` with minima) are in `results/kp7_kept_eval_2026-11-11/{sweep,confirm}/`.

**What the evaluators are** (`hash/kp7kept.h`). With complement 3 − c and val(3 − c) = −val(c), true for both the
Gaussian and the Eisenstein encodings, the per-digit value is affine in the digit's two bits. The b0·b1 coefficient,
val(3) − val(2) − val(1) + val(0), is 0. So:
- **chunk tables can use any bit width:** `ec`/`mc` with B = 8, 10, 11, 12, 13 or 16 bits per lookup, ⌈64/B⌉ lookups
  per strand;
- **one table of subset sums of x^t serves both bit planes and both strands:** `ep`/`gp`/`mp`, which take 4 `pext`
  and 4·⌈32/B⌉ lookups;
- **Gaussian planes need only additions,** because (1 ∓ i) is free in ℤ[i];
- **Eisenstein planes need one multiply by ω per strand.**

`gp` is the real-x idea: with x real (KP7_XREAL) the subset sums are real 8-byte entries. CM31 cannot use a real x,
because its collision bound would become 2(k − 1)/p ≈ 3 × 10⁻⁸.

## Results

### Per kept k-mer: the four ideas

Values are the slope between 1/16 and 1/4, batched collection, batch 512 unless stated. Both jobs are shown, so the
between-job spread is visible.

| idea | evaluator | table | lookups per kept | sweep (3dddf48c) | confirm (88c9cc63) |
|---|---|---|---|---|---|
| (current) | `lib-e 8`, Eisenstein 4-base chunks | 16 KB | 16 | 25.3 | **21.8** |
| chunk size | `ec 11` | 96 KB | 12 | 26.3 | 24.0 |
| chunk size | `ec 13` | 320 KB | 10 | 26.0 | 23.8 |
| chunk size | `ec 16` | 2 MB | 8 | 30.3 | — |
| batching | `ec 8`, batch 64 / 4,096 | 16 KB | 16 | — | 27.0 / 27.0 |
| real-x planes | `gp 11` / `gp 16`, Gaussian 61-bit, x real | 48 KB / 1 MB | 12 / 8 | 49.9 / 54.1 | — |
| planes | `ep 11` / `ep 16`, Eisenstein | 48 KB / 1 MB | 12 / 8 | 44.8 / 49.5 | — |
| CM31 16 KB tables | `lib-m 8`, kp7pack_sample_m | 16 KB | 16 | 59.9 | 55.7 |
| CM31 planes | `mp 11` / `mp 16` | 48 KB / 1 MB | 12 / 8 | 57.7 / 65.4 | — |
| **baseline** | `nt 8`: ntHash from rotated-seed tables | 16 KB | 16 | 21.4 | 20.2 |
| **baseline** | `uh`: min(fwd, rc), multiply-shift | — | 0 | 7.3 | 5.9 |

Reading:
- **Fewer, larger lookups do not pay.**
  - 11- and 13-bit chunks save 4–6 loads but measure within noise of the 16 KB table (±3 cycles between jobs on the
    same code: `ec 8` = `lib-e 8` measured 29.9 against 25.3 in the sweep).
  - Tables beyond L1 pay for it. At 1/4, every 16-bit table and `ec 13` took about 1.58e9 L1-dcache load misses per
    run (sweep counters.tsv), against roughly 0.05–0.14e9 for the 16 KB tables. L2 and L3 misses were not counted.
- **Planes are worse.** Together the 4 `pext` (port 1 only), the plane combination and, for Eisenstein, the two
  ω-multiplies cost more than the 4–8 lookups they save. Real-x Gaussian adds the four-product F_{q²} multiply, 50–54
  cycles. That answers the brief's real-x item: it is the cheapest form of the Gaussian hash, but it stays about 2×
  the Eisenstein one.
- **Batching:** neither 64 nor 4,096 beats 512.
- **The CM31 16 KB tables remain about 2.5× the Eisenstein cost.** Their 31-bit sums are split into re and im, and
  the product takes 4 multiplies and 4 reductions.
- **The gap to ntHash-on-kept is small: 21.8 against 20.2 in the same job.** Both do the same 16 lookups. KP-7 adds
  one 61-bit multiply and two partial reductions, so the KP-7 evaluator is now within about 2 cycles of the
  equivalent table form of ntHash. This is the "arithmetic-bound" floor: lookups plus the adds that consume them.
  Only the universal hash (no tables) is far below it.

### Per input k-mer against batched ntHash (confirm job, median of 5)

| method | 1/1000 | 1/100 | 1/16 | 1/4 |
|---|---|---|---|---|
| batched ntHash (every k-mer) | 11.19 (1.00×) | | | |
| UH on every k-mer (provable, k ≤ 32) | 9.91 (1.13×) | | | |
| current Eisenstein sampler, branchy (`lib-e 8 br`) | 9.14 (1.22×) | 10.41 (1.08×)* | — | — |
| current Eisenstein sampler, batched (`lib-e 8 bat`) | 9.68 (1.16×) | 9.93 (1.13×)* | 11.05 (1.01×) | 15.15 (0.74×) |
| **`ec 13`, branchy** | **8.91 (1.26×)** | **9.44 (1.18×)*** | — | — |
| `ec 13`, batched | 9.76 (1.15×)* | 9.74 (1.15×)* | 10.89 (1.03×) | 15.34 (0.73×)* |
| CM31 16 KB tables, batched | — | — | 13.39 (0.84×) | 23.83 (0.47×)* |
| gate + ntHash on kept, branchy / batched | 9.51 (1.18×)* | 10.12 (1.11×) | 10.80 (1.04×)* | 14.58 (0.77×)* |
| gate + UH on kept, branchy / batched | 8.92 (1.25×)* | 9.35 (1.20×)* | 9.99 (1.12×)* | 11.09 (1.01×)* |

\* At least one of the 5 runs was much slower than the rest (spread above 15%; see Caveats). The minima
(`confirm/table_min.md`) give the same ordering:
- `ec 13` branchy 9.22 / 9.76 (1.23× / 1.17× against ntHash's minimum, 11.39)
- `lib-e 8` branchy 9.39 / 10.06 (1.21× / 1.13×)
- UH every k-mer 10.19 (1.12×)

**Break-even against batched ntHash.** From the 1/16–1/4 line of the confirm job, the batched Eisenstein sampler
breaks even at about 1/14.6 (intercept 9.69). The paper's job gave 1/20.3 against an ntHash of 11.37. The difference
is between-job variance in the intercept and in ntHash, not a code change.

## What did not work (and why)

1. **Larger chunk tables** (11, 13, 16 bits): fewer loads, same or worse cost. Beyond 16 KB the table no longer stays in L1 at dense rates.
2. **Bit-plane evaluation** (one shared subset-sum table, "real x"): `pext` and the plane combination cost more than
   the lookups saved. With Eisenstein the two ω-multiplies make it worse still. CM31 and Gaussian real-x planes do not
   beat the CM31 16 KB tables or the Eisenstein table.
3. **Batch size** 64 or 4,096 instead of 512: no gain. 64 adds per-batch overhead. 4,096 makes the arrays (96 KB)
   compete with the table for L1.
4. **SIMD (AVX-512 gathers): not built.** 16 gathers per 8 kept k-mers are bound by the same 2 load ports as scalar
   loads. The 61-bit multiply has no vector form. At sketch rates the evaluator is 0.3% of the cost anyway.
5. **A harness bug, found and fixed before the confirm job.**
   - The sweep's sink updated `s->x`/`s->n` through a pointer that may alias the hash buffer. That cost every
     every-k-mer method about 10 instructions per k-mer: batched ntHash measured 57.1 instructions and 13.76 cycles in
     the sweep, against 47.0 / 11.37 in bench_e.
   - Fixed in 8ea3e51: the confirm job measured 48.0 / 11.19.
   - **The sweep's "vs ntHash" column and its UH-every-k-mer row are therefore invalid.** Its per-kept slopes are
     unaffected, because every sampler pays the sink once per kept k-mer.

## Caveats

- **Laila was busy during both jobs.**
  - The load runs spread 1.4% in the sweep and 18% in the confirm job.
  - Single runs came in up to 2.3× slow, with identical instruction counts (e.g. `uh 0 br 1/100` 164%; batched ntHash
    99%: one run of 5).
  - Medians of 5 absorb up to two such runs, and the minima agree on every ordering above. But differences under
    about 0.3 cycles per k-mer (3%), such as `ec 13` against `lib-e 8` batched at 1/1000, are not resolved.
- **Per-kept cost at 1/1000 and 1/100 cannot be measured** from per-input-k-mer counters: noise of 0.1 cycles becomes
  100 or 10 cycles per kept. The `per kept at 1/1000 / 1/100` columns of the tables are noise and should not be quoted.
  The slope from 1/16 to 1/4 is the measurement.
- **The `ec 13` branchy gain is gate-loop code generation, not evaluation.** It is 3 fewer instructions per k-mer in a
  different compilation unit. If the paper adopts it, re-time it inside the paper's consolidated job.
- **ntbatch.h is wrong at k = 32** (found here; unchanged):
  - Its multi-bit `srol(x, d)` disagrees with NtHash::roll() on rolled windows at k = 32 (1 of 29 windows equal in a
    60-base test).
  - k = 31, the only k ever timed, is exact (also 0/736,547,799 differences by hash/nteq).
  - The new ntHash table evaluator `nt` matches the library at both k.
- UH's bound is 2⁻⁶³ per pair of distinct canonical codes (multiply-shift with a random odd 128-bit multiplier,
  Dietzfelbinger et al. 1997). It needs k ≤ 32 and has none of KP-7's algebra.

## What the paper should use

`manuscript/paper.md` was not edited.

- **Keep the paper's speed numbers from job 08653d9e6a7453a2** (1.21× at 1/1000, 1.14× at 1/100, 28.5 cycles per
  kept k-mer, break-even about 1/20). This work found no evaluator that beats them by more than between-job variance.
- **If the paper states the per-kept cost as a range,** "about 22–28 cycles" covers the three jobs (this
  confirm 21.8, sweep 25.3, paper 28.5).
- **Add the round-1 baselines, preferably timed in the consolidated job.** Here (confirm job 88c9cc63e42a63a9):
  - **UH on every k-mer: 9.91 cycles, 1.13× batched ntHash.** A provable (universal) strand-symmetric hash on every
    k-mer for k ≤ 32, faster than batched ntHash. This supports the reframing that, for k ≤ 32, the canonical code
    plus a universal hash is the cheaper provable baseline.
  - **The same gate followed by ntHash on kept k-mers:** 20.2 cycles per kept k-mer, against KP-7's 21.8. With UH:
    5.9.
- To quote "best" sampler rates of 1.26× / 1.18× (`ec 13` branchy), add `ec 13 br` to `manuscript/run_paper_bench.sh`
  first, so that they come from the consolidated job.
