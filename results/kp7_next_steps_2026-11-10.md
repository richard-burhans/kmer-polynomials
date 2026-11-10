# KP-7 next steps: table evaluation, batched ntHash, multi-seed statistics (2026-11-10)

Laila: Xeon Gold 6330, hardware counters via cannabis-genome `scripts/laila_perf.py run`. Jobs 5a03e1322aef25af (all
methods) and 96950d85975bd8b9 (three repeats of the sampled variants, spread ≤ 3%). Each job had its own `__load_only__`
baseline, which was subtracted, then divided by 736,547,799 input k-mers. Raw logs: `laila_counters_next_2026-11-10/`.
kp7 at 8994885 for the counter runs. Roast multiz jobs were running on laila at the same time; the repeats show this
did not move the numbers by more than 3%.

**A benchmark bug, found and corrected.** `only=` matched by substring, so "pack … sample 1/100" also timed its
"… 1/1000" sibling in the same counted process. The two 1/100 rows below are corrected by subtracting the 1/1000 row
from the same job. Both corrected values fit the per-kept cost line. `only=` now prefers an exact name (2291e62).

## Cycles per input k-mer (hashing only)

| method | cycles | instructions | vs ntHash library | vs ntHash batched |
|---|---|---|---|---|
| ntHash library `NtHash::roll()` | 23.58 | 110.1 | 1.00× | 0.48× |
| **ntHash batched** (`ntbatch.h`, bit-identical) | **11.26** (median of 3) | 48.0 | 2.09× | 1.00× |
| v2 gate only, every k-mer emitted | 18.54 | 58.9 | 1.27× | 0.61× |
| v2 sample 1/1000, direct evaluation | 11.15 (median of 3) | 46.6 | 2.11× | 1.01× |
| v2 sample 1/100, direct evaluation (corrected) | 27.2 | 86.8 | 0.87× | 0.41× |
| **v2 sample 1/1000, chunk tables** | **9.59** (median of 3) | 45.6 | **2.46×** | **1.17×** |
| **v2 sample 1/100, chunk tables** (corrected) | **10.64** | 48.5 | **2.22×** | **1.06×** |
| v2 sample 1/16, chunk tables | 18.23 (median of 3) | 65.4 | 1.29× | 0.62× |
| v2 sample 1/16, position tables | 29.12 | 108.6 | 0.81× | 0.39× |
| v2 provable for EVERY k-mer, chunk tables (`kp7pack_sym_all`) | 117.12 | 361.6 | 0.20× | 0.10× |
| CM31 AVX-512, provable for every k-mer | 25.85 | 60.2 | 0.91× | 0.44× |
| v1 sample 1/1000 (rolling field hash + gate) | 33.28 | 142.4 | 0.71× | 0.34× |

From the sampled rows, a v2 k-mer costs a base of about 9.5 cycles plus, per KEPT k-mer, about 120–140 cycles with chunk
tables (≈ 32 table lookups and additions plus one F_{q²} product). Direct Horner evaluation costs about 1,780. Break-even
sampling rates for v2 with chunk tables:
- against batched ntHash: about 1/72 (sparser sampling beats it)
- against CM31 SIMD (provable for every k-mer): about 1/8
- against v1: about 1/5.5

**Verdicts:**
1. **Table evaluation works.** Chunk tables cut the per-kept cost about 13×, from 1,780 to about 130 cycles, and are
   bit-identical to kp7_roll64(SYM) (test_kp7pack: every k from 1 to 32, both variants, all thresholds, N layouts,
   block sizes). Position tables (k lookups) are about twice as slow as chunk tables.
2. **v2 with chunk tables at sketch rates (1/1000, 1/100) is faster than batched ntHash** (1.17× and 1.06×), and every
   kept k-mer carries the provable KP-7 hash.
3. **Provable for every k-mer by chunk tables (117 cycles) does NOT beat the CM31 SIMD version (25.9).** For a
   provable hash on every k-mer, use CM31 SIMD.
4. **Batched ntHash is 2.1× faster than ntHash's own per-call API.** Against it, the CM31 SIMD provable hash is 0.44×
   (not the 0.91× measured against the library), and the v2 gate-only hash is 0.61×.

## Sampling statistics, 5 seeds per method (full cs10, `stats_seeds`, results/kp7_stats_seeds_2026-11-10.tsv)

cs10: 221 records, 876,147,649 bases, 47.7% soft-masked. Sampling rate 2/21; geometric expectation for gaps > 200
k-mers: 0.14 per seed.

| method | max bit deviation (mean) | density ratio w=10/20/50 (mean) | gaps > 200 per seed | of those OUTSIDE repeats | longest gap per seed |
|---|---|---|---|---|---|
| ntHash | 0.00088 | 1.0004 / 1.0014 / 1.0036 | 8, 87, 56, 94, 96 | 0, 0, 0, 1, 1 | 729, 14275, 724, 14290, 14278 |
| v2 gate | 0.00103 | 1.0009 / 1.0013 / 1.0024 | 85, 79, 59, 86, 54 | 0, 0, 0, 0, 0 | 15354, 14295, 739, 14319, 753 |
| v1 gate | 0.00101 | 1.0006 / 1.0014 / 1.0021 | 87, 50, 90, 87, 57 | 0, 0, 0, 0, 0 | 14310, 14287, 14297, 14293, 740 |
| provable | 0.00103 | 1.0005 / 1.0008 / 1.0028 | 59, 53, 82, 49, 89 | 0, 0, 0, 0, 0 | 713, 721, 14284, 732, 14308 |

(A gap is inside a repeat when at least half the bases it spans are soft-masked. ntHash seeds 1–4 use ntHash's own
multi-hash extension, the v2 gate uses mix(fwd ^ s) + mix(rc ^ s), and v1 and the provable hash use kp7_init seeds.)

Reading:
- **Every method has far more long gaps than the geometric expectation (tens per seed, against 0.14), and all but 2 of
  about 1,450 long gaps lie in soft-masked repeats.** The 14,3xx-k-mer gap appears for every method at some seeds: one
  long repeat whose few distinct k-mers either all miss the threshold or don't, depending on the seed.
- **No method has more deserts outside repeats than another.** The v2 gate's 14,310 gap in the earlier single-seed run
  was this repeat, not a weakness of the gate.
- Density and bit uniformity are within seed-to-seed noise across all four.
- On Marçais et al. 2024: all four hashes are strand-symmetric (each folds a k-mer with its reverse complement), so
  this comparison cannot show the folding effect itself. It shows only that among symmetric hashes the long gaps come
  from repeats, not from the hash.

## Which method for which job

| job | use | cycles/k-mer (laila) |
|---|---|---|
| FracMinHash-style sketch at rate ≤ about 1/72, with provable hashes for the algebra (KP-9/KP-10) | v2 with chunk tables (`kp7pack_sample_eval`, `KP7PACK_EVAL_CHUNK`) | 9.6–10.6, faster than batched ntHash |
| sampling between about 1/72 and 1/8 | v2 with chunk tables (beats every provable alternative; slower than batched ntHash) | 10.6–25 |
| a provable hash for every k-mer, or sampling denser than about 1/8 | CM31 SIMD (AVX-512 on laila, AVX2 elsewhere) | 25.9 |
| a fast strand-symmetric hash with no proof | batched ntHash (`ntbatch.h`) | 11.3 |
