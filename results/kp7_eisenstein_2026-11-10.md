# KP-7 over the Eisenstein encoding: one field element per hash (2026-11-10)

**Verdict.** Two things improve:
- **The provable hash on every k-mer:** the Eisenstein encoding (A = 1, C = ω, G = −ω, T = −1, ω ∈ F_q,
  q = 2⁶¹ − 1) halves its cost, from 47.5 to **24.6 cycles/k-mer** (2 lanes, scalar C). That beats the CM31 AVX-512
  version (25.8) without any SIMD.
- **The hybrid-v2 sampler:** a kept k-mer costs about **29 cycles**, against 64 with CM31 tables. The sampler beats
  batched ntHash below about **1/19** (it was 1/36). It beats every provable every-k-mer method below about 1/2.

Everything is bit-identical to the scalar reference.

Code: `hash/kp7e.{h,c}` (7ed553d), tests `hash/test_kp7e.c`, bench `hash/bench_e.cpp` (e636d3b).

## Why it is cheaper

ω is a primitive cube root of unity. It lies in F_q because q ≡ 1 (mod 3), while i does not (q ≡ 3 mod 4), which is
why the Gaussian encoding needs F_{q²}. So F, R and the symmetric hash S = (F + λ)(R + λ) are each ONE 61-bit field
element:
- rolling with the λ state and merged 16-entry tables costs one modular multiply per strand plus one for S, so
  **3 multiplies per k-mer**;
- the Gaussian x-real version needs 4 plus a full F_{q²} product;
- the v2 tables hold 8-byte entries, so the table is **16 KB**: 8 lookups and additions per strand, one partial
  reduction, one multiply.

## Guarantees

- **Collision bound.** The four base values are distinct in F_q, so distinct k-mers give a nonzero difference
  polynomial of degree ≤ k − 1, and Pr[F(s) = F(u)] ≤ (k − 1)/q. As a polynomial in λ, S factors into (λ + F)(λ + R),
  so it determines the double-stranded class. By Schwartz–Zippel, distinct classes collide with probability
  ≤ 2(k − 1)/q ≈ 2.6 × 10⁻¹⁷ at k = 31: the difference has total degree ≤ 2k − 2 because the λ² terms cancel (the
  looser 2k/q ≈ 2.7 × 10⁻¹⁷; corrected 2026-11-11, round-1 adjudication L1 F2).
  - This is the same as the Gaussian KP7_XREAL bound and comparable to CM31 with x complex.
  - It is not the F_{q²} x-complex bound 2(k − 1)/q².
- **Exactness, checked exhaustively:** at k = 7 and 8, F is injective on all 4^k k-mers, and S gives exactly 8,192 and
  32,896 values, one per double-stranded class.
- **Distances.** In the complex embedding (ω = e^{2πi/3}) the squared distances are 1 for transitions (A–G, C–T), 3 for
  the A–C and G–T transversions and 4 for complements. test_kp7e checks this. That metric lives in ℂ only: F_q has no
  conjugation or absolute value, so the field hash keeps the algebra (complement = negation, ω² + ω + 1 = 0), not the
  metric.

## Tests (`make test_e`, `make mutation_e`)

- field: ω³ = 1 and ω² + ω + 1 = 0; distinct base values; complement = negation
- rolling == direct at every window (k = 31, 8, 1, with N); F(rc s) = R(s); S strand-symmetric
- `kp7e_sym_seq` with 1 and 2 lanes bit-identical (hashes and positions) to the reference roller, for block sizes
  1/3/5/64/4096 and every N layout
- table evaluation == direct for k = 1–32; the gate equals `kp7pack_gate`, so a sampler keeps exactly the k-mers
  kp7pack's samplers keep
- both samplers == {gate < thr}, with the reference S, for 4 thresholds, 4 buffer sizes and 3 N layouts
- mutants: ROLL (sign of the incoming base on the reverse strand), LANE (remainder offset) and TAB (λ missing from the
  table) are all killed
- the bench's identity checks also pass (L1/L2 against the roller; branchy against batched samplers; sampler hashes
  against the roller at the gate positions)

## Laila hardware counters

Laila is a Xeon Gold 6330, job 7fc2a700073c1f4f. Each method ran in its own counted process
(`bench_e cs10.fa 31 1 only=<name>`). The same-job `__load_only__` baseline (7.00e9 cycles, 20.75e9 instructions)
was subtracted, then divided by 736,547,799 input k-mers. Raw: `laila_counters_eisenstein_2026-11-10/`. The
counter cycles are core clock cycles, which ran at about 2.2 to 3.0 GHz with turbo, not at the 2.0 GHz base frequency, so the ratios are what count (round-1 adjudication L4 F3).

| method | cycles/k-mer | instructions | IPC | vs batched ntHash |
|---|---|---|---|---|
| ntHash batched (`ntbatch.h`) | 11.41 | 47.0 | 4.12 | 1.00× |
| ntHash library roll() | 23.44 | 110.1 | 4.70 | 0.49× |
| **Eisenstein, provable for every k-mer** | | | | |
| e scalar reference roller | 76.09 | 172.0 | 2.26 | 0.15× |
| e sym_seq, 1 lane | 27.95 | 72.6 | 2.60 | 0.41× |
| **e sym_seq, 2 lanes** | **24.59** | 78.5 | 3.19 | **0.46×** |
| Gaussian 61-bit sym_seq, 1 lane | 50.47 | 160.2 | 3.17 | 0.23× |
| Gaussian 61-bit sym_seq, 2 lanes | 47.51 | 179.3 | 3.77 | 0.24× |
| CM31 AVX-512, 16 lanes | 25.83 | 60.2 | 2.33 | 0.44× |
| **samplers (provable hash for kept k-mers), per input k-mer** | | | | |
| e branchy 1/1000 | 9.41 | 43.0 | 4.57 | 1.21× |
| e branchy 1/100 | 9.97 | 43.8 | 4.40 | 1.14× |
| e branchy 1/16 | 13.28 | 48.8 | 3.68 | 0.86× |
| e branchy 1/4 | 25.32 | 66.6 | 2.63 | 0.45× |
| e batched 1/1000 | 9.95 | 46.0 | 4.62 | 1.15× |
| e batched 1/100 | 10.19 | 46.9 | 4.60 | 1.12× |
| **e batched 1/16** | **11.87** | 51.9 | 4.38 | **0.96×** |
| **e batched 1/4** | **17.21** | 70.1 | 4.07 | **0.66×** |
| CM31 branchy 1/1000 | 9.43 | 43.1 | 4.57 | 1.21× |
| CM31 branchy 1/100 | 10.44 | 44.8 | 4.29 | 1.09× |
| CM31 branchy 1/16 | 15.53 | 54.9 | 3.53 | 0.73× |
| CM31 branchy 1/4 | 33.49 | 90.8 | 2.71 | 0.34× |
| CM31 batched 1/16 | 13.64 | 58.0 | 4.25 | 0.84× |
| CM31 batched 1/4 | 25.71 | 94.3 | 3.67 | 0.44× |

## Per-kept cost and break-even rates

These come from the batched samplers' slope between 1/16 and 1/4:
- **Eisenstein:** (17.21 − 11.87) / (0.25 − 0.0625) ≈ **28.5 cycles per kept k-mer**, base ≈ 9.9.
- **CM31:** (25.71 − 13.64) / 0.1875 ≈ **64 cycles**, which matches results/kp7_kept_eval_2026-11-10.md.

| comparison | Eisenstein v2 (batched) | previous best (CM31 v2) |
|---|---|---|
| beats batched ntHash (11.41) below | **about 1/19** (9.9 + 28.5r = 11.41) | 1/36 |
| beats the best provable every-k-mer hash below | **about 1/2** (against e 2-lane, 24.6) | 1/3.9 (against CM31 SIMD, 25.9) |

At very sparse rates (1/1000) the branchy sampler is best (9.41, 1.21×), because the branch is almost never taken. At
1/16 and denser, batched collection wins. At 1/4 the branchy sampler takes 218 M branch misses (189 M above the load
run's 29 M), against 30.2 M for the batched one (about 1 M above the load).

## SIMD and AVX-512 IFMA (not built)

A 61-bit modular multiply has no vector instruction. AVX-512 IFMA multiplies 52-bit limbs, so it would need a prime
p < 2⁵² with p ≡ 1 (mod 3) and Montgomery or Barrett reduction (about 3 IFMA per modmul), at 8 lanes. That buys
vector width, but costs a weaker bound (about 2k/p ≈ 10⁻¹⁴ at k = 31) and the Mersenne shortcut. The scalar
Eisenstein 2-lane loop already matches the CM31 AVX-512 code, so the next speed-up for the every-k-mer case would
come from more ILP (3-4 interleaved stretches) before SIMD. Not built.

## Which to use (updated)

- **Sketches at 1/19 or sparser:** the Eisenstein v2 sampler. It is faster than batched ntHash, with a provable hash
  on every kept k-mer. Use the branchy variant below about 1/100 and the batched variant above.
- **Sampling denser than 1/19, up to about 1/2:** the Eisenstein v2 batched sampler, the fastest provable option.
- **A provable hash on every k-mer:** Eisenstein sym_seq with 2 lanes (24.6 cycles), on par with CM31 AVX-512 and
  needing no SIMD.
- A sketch must use one hash family throughout (Eisenstein, Gaussian 61-bit or CM31). Their values differ.
