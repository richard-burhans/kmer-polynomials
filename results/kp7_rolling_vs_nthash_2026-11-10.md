# KP-6/KP-7 rolling hash against ntHash on cs10 (2026-11-10)

**Verdict:**
- **Speed:** the strand-symmetric hash works but is not competitive. It runs at **90.9 M k-mers/s against
  ntHash's 220.5 M (2.4× slower)**, single-threaded.
- **Minimizer quality:** its random-minimizer density matches the theoretical 2/(w+1) and ntHash to within 0.3% at
  w = 10, 20 and 50, and its output bits are uniform.
- **Deserts:** inconclusive. All hashes show repeat-driven long gaps, and one seed per method cannot rank them.

## Setup

- **Input:** cs10 (`hub-build/2bit/cs10_NCBI_RefSeq_softmasked.2bit` via `twoBitToFa`), 221 records,
  876,147,649 bases, upper-cased (case-blind). k = 31, giving 736,547,799 valid k-mers (windows with N skipped)
  for every method.
- **Command:** `hash/run_bench.sh 31 5`. It builds ntHash via `hash/build_nthash.sh` (pinned commit c26bd45,
  bcgsc/ntHash, 2026-04-01), builds and tests kp7, and runs `hash/bench`. The raw output is in
  `kp7_rolling_vs_nthash_2026-11-10.raw.txt`, with per-round times.
- **Compiler:** gcc/g++ 15.2.0, both libraries `-O3 -march=native`. Each library's roll is a non-inlined call into
  its own object file, so both pay the same call overhead.
- **Machine:** Intel Xeon (AVX2), 32 threads available. One thread used.
- **Timing:** the input is loaded once and I/O is excluded. 5 interleaved rounds over all methods; medians reported.
- **kp7:** q = 2⁶¹ − 1, F_{q²} = ℤ[i]/(q), A = 1, C = i, G = −i, T = −1. The 64-bit output is
  re ⊕ rotl(im, 32) (justified in `hash/kp7.h`). "x real" means x₀ ∈ F_q, with collision bound (k−1)/q per pair.
  "x complex" means x₀ ∈ F_{q²}, with bound (k−1)/q².

## Throughput

| method | M k-mers/s | vs ntHash |
|---|---|---|
| ntHash canonical (fwd + rev) | **220.5** | 1.00 |
| kp7 forward only, x real | 135.4 | 0.61 |
| kp7 canonical min(F, R), x real | 122.1 | 0.55 |
| kp7 symmetric (F+λ)(R+λ), x real | **90.9** | **0.41** |
| kp7 symmetric, x complex | 80.9 | 0.37 |

Per k-mer, the symmetric hash costs two rolling multiplications by a real x (2 modular products each) plus one full
F_{q²} product (4), against ntHash's rotations and XORs. Obvious speed-ups remain untried:
- Karatsuba for the F_{q²} product (3 products instead of 4)
- computing the symmetric value only for the k-mers a sketch needs
- vectorizing across sequences

## Minimizer density (random minimizers, leftmost minimum of w consecutive k-mers)

| method | w = 10 | w = 20 | w = 50 |
|---|---|---|---|
| expected 2/(w+1) | 0.18182 | 0.09524 | 0.03922 |
| ntHash canonical | 0.18197 (1.0009) | 0.09520 (0.9996) | 0.03937 (1.0039) |
| kp7 forward | 0.18221 (1.0022) | 0.09558 (1.0036) | 0.03942 (1.0053) |
| kp7 canonical min | 0.18229 (1.0026) | 0.09576 (1.0055) | 0.03945 (1.0061) |
| **kp7 symmetric, x real** | **0.18192 (1.0005)** | **0.09526 (1.0002)** | **0.03934 (1.0032)** |
| kp7 symmetric, x complex | 0.18172 (0.9995) | 0.09532 (1.0008) | 0.03932 (1.0026) |

Ratios to the expectation are in parentheses. The symmetric hash sits at or below canonical-min at every w, and level
with ntHash. The w = 20 gap means are 10.50 (ntHash) and 10.50 (symmetric). The maximum gap is w by construction for
every method.

## Bit uniformity (max per-bit deviation of the fraction of ones from 0.5)

| method | max deviation |
|---|---|
| ntHash canonical | 0.00083 |
| kp7 forward | 0.00191 |
| kp7 symmetric, x real | 0.00132 |
| kp7 symmetric, x complex | 0.00116 |
| kp7 canonical min | 0.25162 |

- **canonical min** is not uniform, by construction: the minimum of two uniform values is skewed toward small
  numbers. At a fixed threshold it therefore samples at about twice the rate (0.185 instead of 0.095). Its
  context-free gap numbers below are **not comparable** to the others.
- **the other kp7 deviations** are about 1-2·10⁻³ over 7.4·10⁸ k-mers. That is larger than pure noise (about
  2·10⁻⁵) but similar in size to ntHash's, and plausibly the same repeat-driven effect: the many identical k-mers in
  repeats all share one hash value, so they weight that value's bits.

## Context-free sampling gaps ("deserts"), threshold rate 2/21

Long gaps (> 200 k-mers) at this rate should essentially never happen in random sequence (expected 0.1):

| method | sampled | gaps > 100 (expected ≈ 3,150) | gaps > 200 | gaps > 400 | max gap |
|---|---|---|---|---|---|
| ntHash canonical | 0.0957 | 3,002 | 8 | 5 | 729 |
| kp7 forward | 0.0974 | 3,518 | 26 | 1 | 623 |
| kp7 symmetric, x real | 0.0950 | 4,473 | 59 | 7 | 713 |
| kp7 symmetric, x complex | 0.0950 | 4,227 | 84 | 32 | 14,298 |

**Every method has long gaps far beyond the geometric expectation, because of genome repeats.** A satellite with
period p has only p distinct k-mers, so it is skipped entirely with probability (1 − r)^p. The 14,298-k-mer gap under
one seed is such a region. Which hash "loses" a given repeat is luck per seed, so these counts measure the genome, not
the hash, and **one seed per method cannot rank the hashes on deserts.** Settling the KP-7 desert claim (Marçais et
al. 2024: canonical k-mers create deserts) needs three things:
- several seeds per method
- canonical-min at a matched sampling rate
- excluding (or separately reporting) masked or repeat regions

## Tests (hash/test_kp7.c, `make test`; `make mutation`)

- rolling equals direct evaluation at every window, on random sequences with about 5% N and lowercase, at k = 1, 8 and
  31, for both variants of x
- F(rc s) = R(s), and sym(s) = sym(rc s), at every window
- exhaustive at k = 7 and k = 8: F is injective, and sym, both as a field element and after the 64-bit fold, has
  exactly one value per double-stranded class (8,192 and 32,896)
- ntHash's canonical hash is also exact at k = 8 (32,896 values; `bench --nthash-classes`). Its non-polynomial
  hashes do not share the degeneracy of the polynomial F + R
- `make mutation`, which flips the sign of one rolling term, fails 149,968 checks: the mutant is killed

## What this means for KP-6/KP-7

- **Correctness and sketch quality:** the algebraic symmetric hash is exact and samples as well as ntHash.
- **Cost:** it is 2.4× slower single-threaded. Its case rests on its algebraic properties, not speed:
  - an exact strand-symmetric value with no comparison
  - a provable collision bound
  - use as a field element in KP-9/KP-10 sketches
- **Next:** Karatsuba plus a lazily computed symmetric value, to measure how much of the gap closes. A multi-seed,
  repeat-aware desert comparison.
