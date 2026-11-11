# KP-7 hybrid v2: cheaper evaluation of kept k-mers (2026-11-10)

**Verdict: the cost per kept k-mer falls from 142 to 64 cycles (CM31, batched, 4-base chunks). That raises v2's
break-even rate against batched ntHash from about 1/71 to about 1/36, and against the CM31 SIMD every-k-mer hash
from about 1/8.5 to about 1/4. Every variant is bit-identical to its reference.**

Code: hash/kp7pack.{c,h} (439cffd, 3c4bce8). Benchmark: hash/bench_pack.cpp (85ed875 and later). Measured on laila
(Xeon Gold 6330) with hardware counters via cannabis-genome `scripts/laila_perf.py run`, in two jobs: 9303509ac7d3567d
and 5deca06fa147bf0e. Raw logs are in `laila_counters_kept_2026-11-10/`. Cycles are per INPUT k-mer, net of the same
job's `__load_only__` run, over 736,547,799 k-mers of cs10, k = 31.

## What changed

1. **Profile first.** In a local cpu-clock profile of "chunk-eval all", 27% of the time was `kp7_sym`, a separate
   out-of-line product that also tests the Karatsuba flag on every call. 52% was the chunk loop, mostly address
   arithmetic, because `nchunk` is a runtime count and the loop wasn't unrolled. After the changes below, the fast
   path is about 63% arithmetic (table load plus add, the reductions) and about 32% loads, with branches and calls
   around 1%.
2. **Fast evaluators**, `kp7pack_ftab` (61-bit, q = 2⁶¹−1) and `kp7pack_mtab` (CM31, p = 2³¹−1, x complex):
   - λ is folded into chunk 0's table, so the sums are G = F + λ and H = R + λ directly
   - all chunks are unrolled; tables for chunks past k are zero, so they add nothing
   - one partial reduction (61-bit: < 2q, then the wide product takes parts < 4q)
   - the product and fold are inline
   - CM31 entries are 8 bytes (re | im << 32), which halves the 4-base table to 16 KB
3. **Chunk width:** 4-base chunks (W = 8) beat 2-base chunks (W = 4) in both fields. At 1/16, 61-bit: 15.94 against
   16.94 cycles; CM31: 13.76 against 14.76. Kept.
4. **Collection:** "batched" collects kept k-mers without a branch (three stores per k-mer), then evaluates them in a
   tight loop. "Branchy" branches on the gate and evaluates inline. Branchy has the lowest base cost but pays a branch
   mispredict per kept k-mer. Batched pays the stores on every k-mer but evaluates fastest. In job 1 the 61-bit
   batched sampler's base cost was 12.06 cycles (the stores, plus probably register pressure from the 128-bit
   product), so for the 61-bit hash the branchy form is the one to use.
5. **The real-x structure (item 3 of the brief):** not pursued. λ is complex, so G and H are complex anyway and the
   final product is a full F_{q²} multiply. With x real, the table entries would have one zero component per base,
   but chunk entries mix A/T and C/G bases, so the table doesn't shrink. The CM31 table already halves the entry size.

## Results (job 5deca06fa147bf0e; batched ntHash 11.12 and 11.11 cycles per k-mer, first and last run)

| method | per kept k-mer | base (cycles/k-mer) | 1/1000 | 1/100 | 1/16 | 1/4 | break-even vs batched ntHash | vs CM31 SIMD (25.85) |
|---|---|---|---|---|---|---|---|---|
| old chunk eval (61-bit) | 141.7 | 9.10 | 9.25 (1.20×) | — | — | 44.53 (0.25×) | 1/71 | 1/8.5 |
| fast 61-bit, branchy | 102.4 | 9.18 | 9.28 (1.20×) | 10.27 (1.08×) | 15.81 (0.70×) | 34.78 (0.32×) | 1/53 | 1/6.1 |
| fast CM31, branchy | 96.5 | 9.08 | 9.18 (1.21×) | 10.04 (1.11×) | 15.42 (0.72×) | 33.20 (0.33×) | 1/47 | 1/5.8 |
| **fast CM31, batched** | **64.4** | 9.31 | 9.37 (1.19×) | 10.59 (1.05×)* | **13.25 (0.84×)** | **25.40 (0.44×)** | **1/36** | **1/3.9** |

Per kept and base come from the 1/1000 and 1/4 points (linear fit). The break-even rate is where base + rate × per
kept equals the competitor. *1/100 for CM31 batched is from job 9303509ac7d3567d (10.59); its 1/1000 there was 9.95,
against 9.37 here, so between-job variance is about 3–6%.

Job 9303509ac7d3567d also measured every k-mer through the evaluators: old chunk 117.2, fast 61-bit 71.1, fast CM31
71.4 cycles per k-mer. Both stay far from the CM31 SIMD rolling hash (25.85), because table evaluation of every
k-mer can't compete with rolling.

## Which hash a sketch uses

- `kp7pack_sample_f`/`_fb` produce the 61-bit KP-7 hash, bit-identical to `kp7_roll64(.., KP7_OUT_SYM, ..)`.
- `kp7pack_sample_m`/`_mb` produce the CM31 KP-7 hash, bit-identical to `kp7m31_direct` + `kp7m31_sym`. That is a
  different function, with collision bound 2(k−1)/p² ≈ 1.3e-17 per pair at k = 31 (corrected 2026-11-11 from 1.4e-17).
- A sketch must use one of them throughout.

## Recommendation

- **For a sketch sampling 1/36 or sparser,** use v2 with CM31 batched (`kp7pack_sample_m`, W = 8). It beats batched
  ntHash, and its provable hash goes on every kept k-mer. At 1/1000 and 1/100 every fast variant is within noise,
  about 1.2× and 1.05–1.1×.
- **If the 61-bit hash is required,** use `kp7pack_sample_fb` (break-even 1/53).
- **For a provable hash on every k-mer, or sampling denser than about 1/4,** use the CM31 SIMD rolling hash
  (results/kp7_m31_simd_2026-11-10.md).
