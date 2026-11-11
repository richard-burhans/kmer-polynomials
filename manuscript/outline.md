# Outline and framing

**One paper, one claim:** a rolling k-mer hash whose strand symmetry and collision behaviour are proven, not measured,
and which runs at practical speed. The case rests on guarantees and algebra. Speed is competitive, and better than
batched ntHash only in the sampling setting.

## Title (working)

Provably strand-symmetric polynomial hashing of k-mers

## Sections

1. **Abstract** (provisional, single source copied into the repository README by `make readme`). To be revised after
   KP-10.
2. **Introduction** (PENDING). Intent: why comparing k-mer content across genomes needs a strand-independent hash.
   What ntHash and ntHash2 provide (fixed seeds, fwd + rev, empirical quality). What a proof buys: a per-pair bound for
   any input, exact strand handling, and a hash that is native to algebraic set sketches. Final framing waits for
   KP-10.
3. **Methods**
   - 3.1 Encoding. The θ = π/3 trigonometric map of Das, Nanda & Das 2017, written as sixth roots of unity.
     Complement is negation, contrasted with Anastassiou's conjugation. Also Mendizabal's real negation map.
   - 3.2 The transition criterion. With complement as negation, transitions separate from A–C/G–T transversions iff
     Re(A·C̄) ≠ 0, so no fourth-root map can. This is ours: the reason for choosing π/3.
   - 3.3 Polynomial hash over F_q. q = 2⁶¹ − 1 ≡ 1 (mod 3), so ω ∈ F_q. The identities and rolling updates.
   - 3.4 The strand-symmetric product S = (F + λ)(R + λ). Unique factorization in λ, the Schwartz–Zippel bound, and why
     F + R and F·R fail.
   - 3.5 Implementation. λ state, merged tables, lazy reduction, batch loop and lanes, then the CM31 SIMD variant and
     the hybrid sampler (packed-code gate, then table evaluation of kept k-mers).
   - 3.6 Testing and measurement. Exhaustive class counts, rolling = direct, bit identity across variants, mutants.
     Laila hardware counters, load baseline, batched ntHash as the baseline.
4. **Results**
   - 4.1 The algebraic shortcuts collapse, and the product is exact.
   - 4.2 Speed against batched ntHash: every k-mer, and the sampler.
   - 4.3 Sampling quality: density, bit uniformity, long gaps in repeats.
   - 4.4 Spaced seeds stay strand-consistent for any mask.
   - 4.5 Algebraic set sketches (KP-9/KP-10) against Shibuya's km-peeler (PENDING).
5. **Discussion.** What is different from ntHash/ntHash2 and what is shared (IDEAS.md). The Marçais scope: deserts
   affect window-guaranteed schemes, not threshold sampling. The metric lives in ℂ, the algebra in F_q. Limits.
6. **Supplement S1, negative and niche results:** KP-3 (adjacency is a restatement), KP-4 (DFT bound prunes too
   little), KP-13 (FFT diagonals miss indel-broken alignments), KP-5 (palindromes, niche for long arms), KP-9 on
   haplotypes (difference 42% of the set), and KP-11 (mismatch counting by transform: cite as probably known, pending a
   literature check).

## Owed before submission

- the consolidated laila benchmark (ledger.md, last section)
- KP-10 against Shibuya, which decides sections 2 and 4.5 and the abstract's last sentence
- multi-seed statistics for the Eisenstein provable hash
- a literature check for KP-11 (Krawtchouk/MacWilliams mismatch counting)
- user decisions: co-authors, AI-use disclosure, venue, licence, when to make the repository public
