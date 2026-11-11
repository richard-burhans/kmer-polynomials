# KP-7: which strand-symmetric hash (2026-11-10)

Let F = p_s(x₀), the forward hash, and R = rc(p_s)(x₀) = −x₀^{k−1} p_s(1/x₀), the reverse-complement hash, in
F_{q²} = ℤ[i]/(q), with q ≡ 3 (mod 4). A strand-symmetric hash must give the same value for s and rc(s), and different
values for different double-stranded k-mers.

## Result: use (F + λ)(R + λ), not F + R or F·R

Exhaustive over every k-mer, q = 2³¹ − 1, random x₀ and λ (`prototype/kp7_symmetric_check.py`,
output in `results/kp7_symmetric_check_2026-11-10.txt`):

| candidate | k = 7 (8,192 classes) | k = 8 (32,896 classes) |
|---|---|---|
| e₁ = F + R | 729 values | 6,561 values |
| e₂ = F · R | 4,160 values | 16,362 values |
| (F + λ)(R + λ) | 8,192 (exact) | 32,896 (exact) |

**e₁ is degenerate.** e₁ = Σ sₜ (xᵗ − x^{k−1−t}) depends only on the mirror differences s_t − s_{k−1−t}, and the middle
base cancels for odd k. With s ∈ {±1, ±i}, a difference takes 9 values, so e₁ has at most 9^⌊k/2⌋ values (9³ = 729,
9⁴ = 6,561, as measured). Every ordinary palindrome hashes to 0.

**e₂ is the autocorrelation.** e₂ = −x^{k−1} p(x) p(1/x), the aperiodic autocorrelation without conjugation. So it is
invariant under plain reversal (it merges s with reverse(s), confirmed), and it merges homometric sequences too.

**(F + λ)(R + λ) = λ² + λe₁ + e₂ is right.** As a polynomial in λ over F(x), it factors uniquely into its roots, so it
determines the unordered pair {p_s, rc(p_s)}, i.e. the double-stranded k-mer. By Schwartz–Zippel, two distinct
double-stranded k-mers collide with probability ≤ 2(k−1) / q² over random (x₀, λ).

## Rolling, both strands in O(1)

    F_{j+1} = (F_j − s_j)·x⁻¹ + s_{j+k}·x^{k−1}
    R_{j+1} = x·(R_j + s_j·x^{k−1}) − s_{j+k}

Precompute x⁻¹ and x^{k−1}. The symmetric hash then costs 2 additions and 1 multiplication per k-mer more than the
forward hash. It needs no comparison, unlike canonical min(F, R).

**Corrected twice, 2026-11-10.** The first version said this hash avoids the "sketching deserts" of Marçais et al.
2024. The first correction swung too far, saying any strand-symmetric hash creates them. What the paper actually says
(marcais_2024_knonical, source lines 514 and 574): canonicalization, i.e. folding a k-mer with its reverse complement,
can negate the window guarantee of methods that HAVE one (syncmers and the other guaranteed schemes), producing much
larger deserts. The paper's Table 2 caption says "orders of magnitude", but its own table shows 5–8× for three of five
methods at k = 31 (Mykkeltveit 5.5×, syncmer 0 5.4×, syncmer 1 7.6×, chromosome 1). At k = 15 it holds for most
methods (Mykkeltveit 64–83×, syncmer 0 48–72×, Champarnaud and syncmer 1 > 100×) but not for syncmer k/2 (11–13×).
Fractional sampling is not reliably unchanged either: chromosome 10 at k = 31 goes from 24.68 to 57.52 Mb, though
the paper notes high variance for that method. It has "no visible effect" on FRACTIONAL (FracMinHash-style hash-threshold) sampling,
which has no guarantee to lose. So:
- in window-guaranteed schemes, a strand-symmetric hash (this one, min(F, R), or fwd + rev) cannot avoid the
  canonicalization deserts;
- in threshold sampling, which is how this repo uses the hash for sketches, canonicalization costs nothing extra.

ntHash2 already uses fwd + rev as its canonical hash. So KP-7's novelty is the exact product form (F + λ)(R + λ),
the quadratic factor over F_{q²}, which is injective on double-stranded classes with a provable collision bound.

## Next

- implement the rolling form with q = 2⁶¹ − 1 (128-bit products), and compare it against ntHash/ntHash2 on cs10 for
  speed, minimizer density and spacing
- use it in place of canonical codes in the KP-9 prototype
