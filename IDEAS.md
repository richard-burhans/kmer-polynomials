# k-mers as polynomials over the Gaussian integers: ideas to follow

**Status: OPEN, exploration. Started 2026-11-10 from the user's idea. Nothing here is implemented or checked
against the literature yet.**

Each idea has an ID (KP-n) so later notes, code and experiments can cite it. The "check" line under each idea says
what would show whether it is worth pursuing.

---

## 0. The encoding

Bases map to the fourth roots of unity:

| base | A | C | G | T |
|---|---|---|---|---|
| value | 1 | i | −i | −1 |

A k-mer s₀ s₁ … s_{k−1} is the polynomial p(x) = Σ sⱼ xʲ, degree k − 1, coefficients in {±1, ±i}.

Identities (the reason for this choice):

- complement: −p(x), since A ↔ T and C ↔ G are both sign flips
- reverse: x^{k−1} p(1/x)
- reverse complement: rc(p)(x) = −x^{k−1} p(1/x)

**Prior art (corrected 2026-11-10 after the papers were read in full):** the map is NOT Anastassiou's. His complex
representation makes complement a conjugation, and ours makes it a negation. Ours matches the "trigonometric" mapping
in Das et al. 2021 (das_2021_dna-encoding) at θ = π/4, up to a rotation; Das, Nanda & Das 2017 introduced it. Genomic
signal processing evaluates such encodings at roots of unity, i.e. takes the DFT (for example the period-3 signal of
coding exons at x = e^{2πi/3}). The encoding is not new. What might be new is the polynomial-algebra layer built on it
below.

---

## 1. Single k-mers

### KP-1. Comparison as a metric

Treat k-mers as vectors in ℂᵏ. Per base, |a − b|² is:

- 0 for the same base
- 2 for A–C, A–G, C–T, G–T
- 4 for A–T, C–G (complements)

So ‖p − q‖² = 2·(other mismatches) + 4·(complementary mismatches). It is a true metric that charges a
complementary mismatch double. Plain Hamming distance needs four 0/1 indicator vectors instead.

The inner product Re⟨p, q⟩ = Σ Re(aⱼ b̄ⱼ) scores +1 per match, −1 per complement and 0 otherwise.

*Check:* does the 2/4 weighting help or hurt any real task (seed matching, clustering) compared with Hamming?

### KP-2. Edits as polynomial operations

- substitution of a by b at position j: p + (b − a) xʲ
- slide to the next k-mer (drop s₀, append c): (p − s₀)/x + c·x^{k−1}
- reverse complement: −x^{k−1} p(1/x)

### KP-3. De Bruijn adjacency as a polynomial identity

u → v is an edge (v is u shifted by one base, with c appended) if and only if

    x·v(x) − u(x) = c·xᵏ − u₀

The difference has terms only at degree 0 and degree k. Paths and unitigs then correspond to sequence polynomials.

*Check:* is there an algorithm (compaction, path enumeration, bubble detection) that becomes simpler or faster
in this form, or is it only a restatement?

### KP-4. Rotations, the DFT and Parseval lower bounds

Mod xᵏ − 1, multiplying by x rotates the k-mer, and the DFT turns rotation into a phase change. So |DFT| is
rotation-invariant (cyclic k-mers, necklaces).

Parseval: ‖p − q‖² = (1/k) Σⱼ |P(ωʲ) − Q(ωʲ)|². Keeping only a few DFT coefficients gives a **lower bound** on the
distance, which can discard candidate pairs before an exact comparison. This is the GEMINI / DFT lower-bound
technique from time-series search (Faloutsos et al. 1994), the same field the cross-domain yama pilot mined
(bench/data/cross_domain_pilot/field_timeseries.md).

*Check:* on real k-mer pairs, what fraction of candidates do the first m coefficients prune at a given
mismatch threshold?

### KP-5. RC palindromes and inverted repeats

A reverse-complement palindrome satisfies p(x) = −x^{k−1} p(1/x): an anti-self-reciprocal polynomial. That can be
tested algebraically, possibly as a rolling test over a whole sequence.

*Check:* compare speed and simplicity against existing palindrome and inverted-repeat finders.

---

## 2. The arithmetic: hashing in F_{q²}

### KP-6. Gaussian integers mod a prime give a field that holds the encoding exactly

For a prime q ≡ 3 (mod 4), ℤ[i]/(q) is the field F_{q²}, so 1, i, −i and −1 live there natively.
**q = 2⁶¹ − 1 is ≡ 3 (mod 4)**, so we get Mersenne-prime reduction (cheap) with i built in.

Hash: h(s) = p_s(x₀) for a random x₀ ∈ F_{q²}. Two different k-mers collide at most at k − 1 values of x₀, the
standard polynomial-identity bound.

Both strands roll in O(1) per base:

- forward: KP-2's slide
- reverse complement: from rc(p)(x) = −x^{k−1} p(1/x), the RC hash updates in O(1) as well

### KP-7. Strand invariance without choosing a canonical strand

Instead of canonical = min(h(s), h(rc s)), which needs an order, use symmetric functions of the pair:

- e₁ = h(s) + h(rc s)
- e₂ = h(s) · h(rc s)

or the factor (t − h(s))(t − h(rc s)) = t² − e₁t + e₂, which is identical for both strands. That carries into the set
representations of section 3 as double-stranded sets.

**Status 2026-11-10:** the symmetric functions alone are degenerate (e₁ = F + R has ≤ 9^⌊k/2⌋ values, and e₂ = F·R merges a k-mer with its reverse). (F + λ)(R + λ) is exact. See notes/kp7_strand_symmetric_hash.md.

**Status 2026-11-10 (C, cs10, k = 31):** the symmetric hash runs at 90.9 M k-mers/s against ntHash's 220.5 M (2.4× slower), with random-minimizer density within 0.3% of 2/(w+1) and of ntHash. Deserts are inconclusive (repeat-driven, one seed). See results/kp7_rolling_vs_nthash_2026-11-10.md and hash/.

Karatsuba (3 multiplies) tried: no measurable gain (results/kp7_karatsuba_2026-11-10.md). The remaining cost is additions and reductions, not multiplies.

Addition-side changes (merged tables, λ-state, lazy reduction): +18%, to 113.9 M k-mers/s with x real (0.47× ntHash), bit-identical (results/kp7_addside_2026-11-10.md).

Batch loop (`kp7_sym_batch`): 158.7 M k-mers/s with x real, 0.66× ntHash; 127.1 with x complex (results/kp7_batch_2026-11-10.md).

Laila hardware counters: the best version (2 lanes) is 47.4 cycles/k-mer against ntHash's 23.3. The gap is 1.6× instructions and 1.25× IPC (results/kp7_laila_counters_2026-11-10.md).

CM31 SIMD (p = 2³¹ − 1, x complex): on laila, 0.87× ntHash with AVX2 and 0.91× with AVX-512, the provable symmetric hash at about 60 instructions/k-mer. The remaining gap is IPC (results/kp7_m31_simd_2026-11-10.md).

Hybrid v2 (packed 2-bit codes rolled by shifts, the KMC idea, plus provable hashes only for kept k-mers): 11.6 cycles/k-mer at 1/1000 sampling, 2.0× FASTER than ntHash (23.4); the gate alone is 1.26× faster (results/kp7_pack_hybrid_2026-11-10.md).

*Check (KP-6 + KP-7):* compare against ntHash on cs10 for throughput, uniformity of minimizer selection (density and
spacing between chosen positions), and collision rate.

---

## 3. Collections of k-mers

### KP-8. A set as the roots of a polynomial: set operations become algebra

χ_S(t) = Π_{s∈S} (t − h(s)) over F_{q²}.

| set operation | polynomial operation |
|---|---|
| membership of s | χ_S(h(s)) = 0 |
| A ∩ B | gcd(χ_A, χ_B) |
| A ∪ B | lcm(χ_A, χ_B) |
| A ∖ B | χ_A / gcd(χ_A, χ_B) |
| multiset | repeated roots |
| double-stranded set | product of KP-7's quadratic factors |

### KP-9. Set reconciliation: recover the difference at a cost set by the difference (most promising)

Minsky, Trachtenberg and Zippel (2003), characteristic-polynomial set reconciliation:

1. Evaluate χ_A and χ_B at the same m points.
2. Their ratio is χ_A/χ_B = χ_{A∖B} / χ_{B∖A}.
3. Rational-function interpolation recovers both difference sets exactly whenever |A ∖ B| + |B ∖ A| ≤ m.

The cost scales with the size of the **difference**, not the size of the sets. Two haplotypes of one plant, or two
close assemblies, share nearly all k-mers. This would list exactly the k-mers they do not share, from m field
elements per genome.

*Check:* KOMPa vs KOMPb (two haplotypes of one plant). Recover the k-mer difference, verify it against a brute-force
set difference, and measure sketch size and time against difference size.

**Status 2026-11-10: prototyped** ([results](results/kp9_kompa_kompb_2026-11-10.md)). Recovery is exact. For KOMPa vs
KOMPb chr7, D is 42% of the set, so the sketch is as large as the set and KP-9 does not help. For a near-identical
pair (chr7 with 1,000 SNPs) the sketch is 1.13 MB against 308 MB.

### KP-10. Power-sum sketches: linear, mergeable, decodable

Pⱼ = Σ_{s∈S} h(s)ʲ for j = 1…m.

- **linear:** multiset union is addition and multiset difference is subtraction
- **mergeable and streamable:** add k-mers as they arrive, and combine sketches from chunks or machines
- **decodable:** Newton's identities turn power sums into χ, then root-finding recovers A − B when its size is ≤ m

It is an algebraic counterpart to an invertible Bloom lookup table (IBLT), with no hash-table failure
probability. It is also the streaming form of KP-9.

*Check:* the same data as KP-9, and compare against an IBLT at equal sketch size.

### KP-11. Generating functions by k-mer code: concatenation becomes multiplication

F_S(y) = Σ_{s∈S} y^{code(s)}, where code(s) is the base-4 index of s.

- union and intersection are coefficient-wise (OR and AND), not multiplication
- multiplication is **concatenation**: since code(uv) = code(u)·4^{|v|} + code(v), the product
  F_U(y^{4^{|v|}}) · F_V(y) is the set of all concatenations uv, with coefficients counting the ways to form each

That is formal-language algebra on k-mers, natural for spaced seeds and for questions like "which (k₁ + k₂)-mers can
be assembled from these two sets".

*Check:* is there a concrete question (seed design, assembly-graph reachability) where this beats a direct
computation?

---

## 4. Whole sequences

### KP-12. All k-mer hashes of a sequence, both strands at once

Encode the sequence as S(x) of length n. All windows' hashes come from KP-6's rolling update in O(n), for both
strands.

### KP-13. All-offset correlation by FFT, a diagonal finder

S_A(x) · conj(S_B)(1/x) gives, at every offset, the match-minus-complement score of KP-1, in O(n log n) by one FFT.
MAFFT uses FFT correlation in this way to find homologous segments.

For our work, this could rank diagonals or offsets between two pieces before lastz or yama run.

*Check:* KOMPa chr7 against cs10 NC_044378.1. Do the correlation peaks fall where the chains (rerun
2cc2600a989b3f73) say the alignment is?

### KP-14. Spaced seeds as masks

A 0/1 mask polynomial applied coefficient by coefficient gives spaced k-mers, and KP-11 composes them.

---

## 5. Experiments, in suggested order

1. **KP-9/KP-10 reconciliation** on KOMPa vs KOMPb. If it works, comparing many assemblies costs per difference,
   not per genome.
2. **KP-6/KP-7 hash** against ntHash on cs10.
3. **KP-4 Parseval pruning** rate for mismatch search.
4. **KP-13 FFT diagonals** against the KOMPa chains.

## 6. Literature to check before claiming anything

**Acquisition (2026-11-10):** [literature/wantlist.md](literature/wantlist.md): 24 papers; IBLT reconciliation of k-mer sets is published (Shibuya et al. 2022); the characteristic-polynomial / power-sum form for k-mers was not found.

- Anastassiou 2001, the complex DNA representation (the encoding in section 0)
- Minsky, Trachtenberg and Zippel 2003, characteristic-polynomial set reconciliation (KP-8, KP-9)
- whether reconciliation or power-sum sketches have already been applied to genomic k-mer sets
- ntHash (Mohamadi et al. 2016) and ntHash2, the rolling canonical hash to beat (KP-6, KP-7)
- IBLT (Goodrich and Mitzenmacher 2011) and its genomics uses (KP-10)
- Faloutsos, Ranganathan and Manolopoulos 1994, DFT lower bounds (KP-4)
- the MAFFT FFT step (Katoh et al. 2002) (KP-13)

## 7. Open questions

- Which ideas give an algorithm that is actually faster or smaller, rather than an elegant restatement?
- Do the 2/4 distance weights (KP-1) mean anything biologically, for example transitions against transversions?
  With this map, A–G and C–T (transitions) both cost 2, and so do A–C and G–T (transversions). Only complements
  cost more, so the metric does not separate transitions from transversions. Would a different root assignment
  do that while keeping complement = negation?
- How large must m be in KP-9/KP-10 for two real haplotypes at k = 21 or 31?
- What does the k-mer difference between two haplotypes look like? Is it dominated by SNPs (k k-mers per SNP)?
