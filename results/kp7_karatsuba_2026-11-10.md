# KP-7: Karatsuba for the F_{q²} multiply (2026-11-10)

**Verdict: no measurable gain. Karatsuba stays in the code behind `params.karatsuba` (default off) and is not used.**

Change: (a+bi)(c+di) with imaginary part (a+b)(c+d) − ac − bd, 3 widening multiplies instead of 4 (hash/kp7.c
`fmul_k`). It gives bit-identical values:
- test_kp7 checks edge values, 1 M random pairs, and every rolled window for both variants
- the benchmark checks that the XOR of all 736.5 M hashes matches
- `make mutation` kills a dropped term

Run: `hash/run_bench.sh 31 5 throughput` on cs10, single thread, k = 31, 5 interleaved rounds, `-O3 -march=native`
for kp7 and ntHash alike. Raw output: `kp7_rolling_vs_nthash_2026-11-10.throughput.raw.txt` (kp7 commit 0d5e66b).

| method | median M k-mers/s | range over 5 rounds |
|---|---|---|
| ntHash canonical | 243.7 | |
| kp7 symmetric, x real | 96.1 | 93.0–98.6 |
| kp7 symmetric, x real, Karatsuba | 97.0 | 93.5–100.4 |
| kp7 symmetric, x complex | 86.1 | 81.1–88.9 |
| kp7 symmetric, x complex, Karatsuba | 83.5 | 79.2–85.3 |

The differences (+0.9% with x real, −3% with x complex) lie inside the round-to-round spread.

**Why it does not help:** with x real, Karatsuba removes 1 of the 8 widening multiplies per k-mer. A 64×64→128
multiply is cheap on this CPU, and the extra additions and wider operands offset the saving. The multiplies are not
the bottleneck.

**Likely costs, still to try:**
- the conditional reductions in every `addq` (2 per `fadd`, 4 `fadd` per k-mer, plus 2 for λ)
- the two separate table lookups per strand: out_f[co] and in_f[ci] could merge into one 16-entry table [co][ci],
  saving one `fadd` per strand
- adding λ every k-mer: keep G = F + λ and H = R + λ as the rolling state instead. The update constant absorbs λ(1 −
  x⁻¹) and λ(1 − x), which saves 2 `fadd` per k-mer.
- lazy reduction: keep sums in [0, 2q) or [0, 4q) and reduce once before the multiply
