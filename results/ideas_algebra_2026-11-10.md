# Untested ideas, algebra side: KP-3, KP-14 and the root assignment (2026-11-10)

Code is in `ideas/algebra/` (`make test`, `make mutation`, `./dbg`, `./bench_spaced`, `python3 roots.py`). Raw outputs
are in `ideas_algebra_2026-11-10/`. Timings are local, single thread, `-O3 -march=native` for both sides, on a shared
host, so treat the ratios as approximate.

## KP-3: de Bruijn adjacency as a polynomial identity. Verdict: correct and provable, but slower and bigger than packed codes; a niche, not a speed-up

From u → v iff x·v(x) − u(x) = c·xᵏ − u₀, a k-mer's successor hashes follow from its hash F(u) ALONE as 16 candidates,
(F(u) − b(g))·x⁻¹ + b(c)·x^{k−1}, over the guessed first base g and the appended base c. Predecessors work the same
way. A wrong guess leaves a nonzero constant term, so it matches a stored k-mer only at a root of a nonzero
polynomial of degree ≤ k. That gives a provable bound of ≤ 12·|S|·k/q expected false successors per node, and there
are no false negatives.

`dbg` on the first 20 Mb of cs10 (13,530,525 distinct forward 31-mers, a directed de Bruijn graph):

| | packed 2-bit codes | hash-only (F per node) |
|---|---|---|
| edges | 13,608,457 | 13,608,457 true, 0 false, 0 missed |
| nodes whose out/in-degree differs | — | 0 |
| non-branching nodes | 13,391,419 | 13,391,419 |
| time | 140 ns/node | 1,039 ns/node (7.4× slower) |
| storage per node | 8 bytes | 16 bytes (one F_{q²} element) |

The provable bound on expected false successors for the whole set is 0.03.

- **Packed codes are cheaper when k-mers are stored.** For k ≤ 32 adjacency is one shift and mask on the stored code,
  so packed codes win on both time and space.
- **What the hash form adds:** navigation and unitig compaction over a set that stores only hash values, with a
  provable false-positive bound, at a constant 16 bytes per node whatever k is. Packed codes take 2k bits, so the
  hash form uses less storage only for k > 64.
- **ntHash supports the same navigation.** Its rolling is invertible and linear over GF(2), so the 16-candidate trick
  works there too, but with no proven bound.
- **So KP-3 is mostly a restatement with a niche:** hash-only, large-k graph navigation with a guarantee.

## KP-14: spaced seeds as masks. Verdict: works, and is strand-consistent for ANY mask (ntHash2 is not); half ntHash2's speed

- **Rolling** (`spaced.h`): F_m(j+1) = x⁻¹·F_m(j) − Σ_blocks s_{j+a}x^{a−1} + Σ_blocks s_{j+b}x^{b−1}, with the matching
  formula for the reverse-complement side R_m.
  - Every care block shares ONE multiply per strand.
  - Each block adds two table entries; a base times a power of x is a table lookup.
- **Strand symmetry:** S_m = fold((F_m+λ)(R_m+λ)), where R_m is F_m of the reverse complement. So S_m is
  strand-symmetric for every mask, palindromic or not.
- **Tests:**
  - rolling equals direct evaluation for 6 masks with N;
  - the all-ones mask is bit-identical to `kp7_roll64` symmetric output;
  - strand symmetry holds;
  - exhaustive k = 7: distinct hashes equal the exact masked double-stranded classes, for 4 masks (e.g. 2,080 = 2,080);
  - the mutant is killed.
- **Against ntHash2's `SeedNtHash`** (same masks, one hash per seed, 58.8 Mb of cs10):
  - **Speed:** ours runs at 27.4 against 59.3 M hashes/s for one mask (0.47×), and 0.46× for two masks.
  - **Sensitivity to substitutions:** the fraction of windows unchanged at SNP rates 0.001 / 0.01 / 0.05.
    - **Palindromic mask:** ours 0.97635 / 0.78574 / 0.29189, identical to SeedNtHash and to the exact fraction.
      Sensitivity is a property of the mask, not the hash.
    - **Non-palindromic mask:** ours 0.97150 / 0.74724 / 0.22572, exactly the exact double-stranded fraction.
      SeedNtHash gives 0.98025 / 0.81795 / 0.35843 and prints "Seed … is not symmetric, reverse-complement hashing will
      be inconsistent". Its canonical spaced hash is strand-consistent only for palindromic seeds.
  - **Window counts:** ours 41.28 M against 41.85 M. We skip any window with N in any of its k positions; SeedNtHash
    appears to keep windows whose N falls only on don't-care positions.
- **Many masks at once:** no algebraic shortcut beyond what is implemented.
  - Each mask needs its own state, and one multiply per strand per mask.
  - F_m = F − F_D (the don't-care positions) lets a sparse-don't-care mask reuse one full rolling hash, as ntHash2
    does with ignore blocks.
  - An NTT evaluation of a mask over all positions is possible over CM31 (its multiplicative group has 2-power roots
    of unity), but at O(log n) per position it doesn't beat O(blocks) rolling for k ≤ 32.

## Root assignment (IDEAS section 7). Verdict: yes. Sixth roots of unity (Eisenstein) separate transitions; fourth roots provably cannot

With complement = negation (T = −A, G = −C), both transitions always have the same distance |A + C|. The transversions
have |A − C| (A–C, G–T), 2|A| (A–T) and 2|C| (C–G). Since |A+C|² − |A−C|² = 4·Re(A·C̄):
- transitions separate from the A–C/G–T transversions **iff Re(A·C̄) ≠ 0**;
- they are the shorter **iff the angle between A and C exceeds 90°**.

The current fourth-root map (A = 1, C = i) is exactly the orthogonal case, so it cannot separate them. Every squared
distance is 2, except complements at 4.

**Sixth roots** (A = 1, C = ω = e^{2πi/3}, G = −ω, T = −1, i.e. the Eisenstein integers ℤ[ω]): squared distances 1
(transitions), 3 (A–C, G–T) and 4 (complements). Transitions become the cheapest substitution, a Kimura-like
weighting. `roots.py` checks this, and checks the 1–179° sweep against the Re(A·C̄) criterion.

**Rings.** A primitive cube root ω lies in F_q itself when q ≡ 1 (mod 3). Both Mersenne primes used here qualify
(2⁶¹ − 1 and 2³¹ − 1, each ≡ 1 mod 3; ω found and checked). i would need q ≡ 1 (mod 4), and both are ≡ 3 (mod 4).
So the Eisenstein encoding puts every hash value in F_q, ONE field element instead of the Gaussian encoding's two:
- complement is still negation, and the rc identity R = −x^{k−1}p(1/x) still holds;
- (F + λ)(R + λ) is still exact: 8,192 / 32,896 distinct values for 8,192 / 32,896 classes at k = 7 / 8;
- a base times a power of x is still a table lookup;
- the Euclidean weighting lives only in the complex embedding (F_q has no conjugation);
- the collision bound is (k−1)/q (x ∈ F_q), the same as today's x-real 61-bit hash.

**Follow-up worth measuring (not done here):** an Eisenstein-in-F_q KP-7 roller costs about 3 F_q multiplies per
k-mer (one per strand plus the product), against about 7 F_q-equivalent multiplies for the Gaussian x-real version.
It could be substantially faster at the same collision bound.
