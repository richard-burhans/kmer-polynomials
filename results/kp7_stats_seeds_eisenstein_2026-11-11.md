# Multi-seed sampling statistics for the Eisenstein provable hash (2026-11-11)

**Verdict.** On all of cs10, with 5 seeds, the Eisenstein provable hash S = (F + λ)(R + λ) over F_q samples like
ntHash and like the other hashes:
- **Bit balance:** the mean max per-bit deviation from ½ is 0.00090, against 0.00088 for ntHash.
- **Minimizer density:** within 0.33% of 2/(w + 1) at w = 10, 20 and 50, against 0.40% for ntHash.
- **Long threshold-sampling gaps (> 200 k-mers, rate 2/21; 0.14 expected per seed by chance):** 367 over 5 seeds
  (49–87 per seed), against 341 for ntHash, and **none outside soft-masked repeats** (ntHash: 2).

Raw: `kp7_stats_seeds_2026-11-10.eisenstein.tsv` (the file name carries the runner's local date). Run:
`hash/run_stats_seeds.sh 31 5 eisenstein` at 27b09e3 (`stats_seeds` method `eprovable`, `kp7e_init(k, 0x5eed +
1000·seed)`, then `kp7e_sym_seq` with 2 lanes). Local and untimed; 221 records, 876,147,649 bases, 47.7% soft-masked.

S is a 61-bit field element, so it is fed to the 64-bit statistics as (S << 3) | (S >> 58). Every bit then carries a
uniform field bit, and orders (minimizers, thresholds) are those of S up to its low 3 bits.

| method | mean max bit deviation (range over seeds) | max density deviation | long gaps, 5 seeds (per seed) | outside repeats | longest gap |
|---|---|---|---|---|---|
| ntHash | 0.00088 (0.00065–0.00102) | 0.40% | 341 (8–96) | 2 | 14,290 |
| packed-code gate v2 | 0.00103 (0.00077–0.00128) | 0.31% | 363 (54–86) | 0 | 15,354 |
| packed-code gate v1 | 0.00101 (0.00075–0.00121) | 0.36% | 371 (50–90) | 0 | 14,310 |
| Gaussian provable (F_{q²}) | 0.00103 (0.00084–0.00132) | 0.33% | 332 (49–89) | 0 | 14,308 |
| **Eisenstein provable (F_q)** | **0.00090 (0.00078–0.00110)** | **0.33%** | **367 (49–87)** | **0** | 14,295 |

The first four rows reproduce results/kp7_stats_seeds_2026-11-10.tsv exactly (same seeds, deterministic). Over all
five methods, 1,774 long gaps, of which 2 lay outside soft-masked repeats. The long gaps are repeats, not a property
of any hash.
