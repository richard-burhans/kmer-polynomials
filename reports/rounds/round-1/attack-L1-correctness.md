# Round 1 attack, slot L1, lens: correctness

Attacked text: `reports/rounds/round-1/ATTACKED-paper.md` (identical to `manuscript/paper.md`, draft d049b77).
Scratch scripts: `/home/agent/.claude/jobs/0d0eac3e/tmp/L1/` (`exh.py`, `fr.py`, `sm.py`, `ovf.py`, and a copy of `hash/` for the mutation run).

Prior rounds: `scripts/round history "strand-symmetric"` and `"overflow"` return no prior mention. This is round 1 for this repository, so no earlier adjudication applies to any finding below.
Method cards: the brief carried no appendix. `kmer-polynomials` has no method-card directory. The only hit is `library/reports/method_cards_audit.md`, and none of its cards covers Schwartz–Zippel or field arithmetic. No finding below rests on a card.

**Steering note.** The claim list does not point at the weakest part of C1, which is how the code draws x and λ (F1). It also restates the bound as "~2k/q", the ledger's convention, while the paper's own derivation gives 2(k−1)/q, and the two disagree in the printed digit (F2).

## F1 -- the implementation derives λ from x, so the Schwartz–Zippel premise fails for every reported number (C1, C18)

The proof needs x and λ to be independent and uniform over F_q ("over the random choice of x and λ"). The code does not do this. `hash/kp7e.c:46-62` draws both from one 64-bit seed with two consecutive `splitmix64` outputs, `x = rand_q(&s); lam = rand_q(&s);`. The splitmix mixer is a bijection, so x fixes the seed up to the 3 masked-off high bits, and then x fixes λ up to 8 values. `sm.py` inverts the mixer for the bench seed `0x5eedE` (`hash/bench_e.cpp:128`). From x alone it recovers 8 candidate λ, and the true λ is among them. `kp7m31.c:83-84` (CM31) does the same.

So the polynomial argument does not cover the hash that was measured. In that hash, S(s)−S(u) = λ(x)·a(x) + b(x), where λ(x) is a non-polynomial mixing function, and Schwartz–Zippel does not bound that.

Every experiment also uses fixed seed constants (`0x5eedE`, `0x5eed + 1000·seed`). That weakens the contrast the abstract and C18 draw: "ntHash fixes its seeds, so its collisions are the same in every run" describes these runs too.

**Correction.**
- Draw x and λ independently, each from at least 61 bits of fresh entropy, and say the bound assumes this.
- Add one sentence: "The experiments use fixed seeds; the bound applies to a deployment that draws them at random."

## F2 -- 2.7×10⁻¹⁷ is not 2(k−1)/q, and 1.4×10⁻¹⁷ is neither bound for CM31 (C1)

The Methods text gives "probability at most 2(k−1)/q … At k = 31 that is about 2.7 × 10⁻¹⁷". But 60/(2⁶¹−1) = **2.60×10⁻¹⁷**. The figure 2.7×10⁻¹⁷ is 62/q, i.e. 2k/q, which is what the abstract uses ("about 2k/q").

The CM31 bound "about 1.4 × 10⁻¹⁷" (Implementation; source `results/kp7_m31_simd_2026-11-10.md:12-13` and `results/kp10_vs_shibuya_2026-11-11.md`) is also wrong. With p² = (2³¹−1)², 2(k−1)/p² = **1.30×10⁻¹⁷** and 2k/p² = **1.34×10⁻¹⁷**. Only 64/p² rounds to 1.4. All values were recomputed with Python integers.

**Correction.** Methods: "at most 2(k−1)/q, about 2.6 × 10⁻¹⁷ at k = 31". Abstract: "at most 2(k−1)/q ≈ 2.6 × 10⁻¹⁷". Implementation: "about 1.3 × 10⁻¹⁷". Fix the ledger rows to match.

The algebra itself holds, and I checked it:
- The difference S(s)−S(u) = λ(ΔF+ΔR) + Δ(FR) is the zero polynomial only if {p_s, rc p_s} = {p_u, rc p_u}.
- The four base values are distinct in F_q.
- The total degree is max(k, 2k−2).

## F3 -- "determines the unordered pair … which determines the double-stranded k-mer" is true of polynomials, false of values (C1)

As a single field value, S does not determine {F, R}. One equation in two unknowns has a solution for every F. And the values {F, R} do not determine the k-mer, because F collides with probability (k−1)/q. The sentence moves from polynomials to values without saying so. Separately, "The product with λ is exact at both sizes, as the factorization argument requires" misstates what the k = 7/8 check shows. The argument requires nothing of a single draw. What the check does show is useful: no pair of classes has S(s)−S(u) ≡ 0, since such a pair would collide at every (x, λ). I reproduced both counts independently for the Eisenstein encoding (`exh.py`: 8,192 and 32,896, with F injective on all 4^k k-mers), and `test_kp7e` reproduces them too.

**Correction.**
- Replace the first sentence with: "As polynomials in (x, λ), S(s) and S(u) are equal only if {p_s, rc p_s} = {p_u, rc p_u}, i.e. s and u are the same double-stranded k-mer."
- Replace the second with: "…is exact at both sizes, so no pair of classes has identically equal S."

## F4 -- overflow detection is not deterministic: a constructed bucket with 140 differences decodes as zero (C16)

The paper says "a bucket holding more than m differences is detected rather than resolved wrongly". The Discussion says the sketch "detected its own overflow deterministically". Both are false. The c extra power sums only make a wrong decode unlikely for generic hash values.

Counterexample (`ovf.py`, which calls the paper's own decoder `recon/recon61.py:decode` with m = 64, c = 2):
1. Side A holds the 70 values a·ζ^i and side B the 70 values b·ζ^i, where ζ is a primitive 70th root of unity (70 divides q−1). This is a bucket with 140 differences.
2. All power sums S_1…S_66 of the difference are zero, and δ = 0.
3. `decode` returns `(0, 0)` flagged elements, i.e. "no differences". It does not return `None`, which would mean overflow.

This exact configuration is very unlikely for random hash values, and that is the point: detection holds only with high probability over the hash seed. With fixed seeds (F1) that probability is not even a random event.

**Correction.**
- Methods: "…is detected except with negligible probability over the random seed".
- Discussion: delete "deterministically".
- `results/kp10_vs_shibuya_2026-11-11.md` also says "KP-10 detects an overflow deterministically", so fix it there too.

## F5 -- the "spaced seed" responds to 29 positions, not its mask's 20; the 0.747 vs 0.818 sentence hides this (C14)

The symmetry claim is right only because `ideas/eisenstein/spaced_e.h:4` defines R_m = F_m(rc(window)), which reads the mirrored mask. The cost is that S_m depends on mask ∪ mirror(mask). For the non-palindromic mask `1101100111010111011100101110101` (`bench_spaced_e.cpp:37`), the weight is 20 and the union has **29** positions. The numbers match exactly:

| | measured | prediction |
|---|---|---|
| ours | 0.74722 | 0.99²⁹ = 0.74717 |
| SeedNtHash | 0.81794 | 0.99²⁰ = 0.81791 |
| ours, rate 0.001 | 0.97146 | 0.999²⁹ |
| SeedNtHash, rate 0.001 | 0.98023 | 0.999²⁰ |

So the hash has lost 9 of the mask's don't-care positions. Substitutions there change the hash, and a spaced seed exists to be blind to them. The paper calls 0.747 "the exact double-stranded value" and sets it against 0.818 without saying which is better or why.

The test never checks SeedNtHash's strand consistency directly. It never compares h(w) with h(rc w). The "inconsistent" reading rests only on ntHash2's warning, plus the inference that a hash depending on 20 forward positions cannot be strand-symmetric.

**Correction.** "Strand symmetry forces the hash to depend on the mask and its mirror, 29 positions for this weight-20 mask, so it kept 0.747 = 0.99²⁹ of windows unchanged. SeedNtHash kept 0.818 = 0.99²⁰ because it reads only the forward mask's 20 positions, which is why it is not strand-consistent."

## F6 -- 1/20 is above the measured break-even, and 28.5 cycles per kept k-mer holds only in collected mode (C11)

The source puts the break-even at 1/20.3 (`results/paper_bench_2026-11-11.md:47`), and 1/20 = 0.05 is a higher rate than 1/20.3 = 0.0493. At 1/20 the paper's own line gives 9.96 + 28.48 × 0.05 = **11.38 cycles**, against 11.37 for batched ntHash. Interpolating instead between the measured batched rows at 1/100 (10.22) and 1/16 (11.74) also gives 11.38. So at 1/20 the sampler ties ntHash or is marginally slower. It is not faster.

The batched intercept is consistent with the measured 1/100 and 1/1000 rows, which the paper does not table. So the extrapolation is sound, but its endpoint is 1/20.3, not 1/20.

The inline mode used for the sparse table rows costs (10.01 − 9.40) / 0.009 ≈ **68 cycles per kept k-mer**, so "a kept k-mer costs about 28.5 cycles" holds only for collected evaluation. The margin at 1/20 (0.1%) is inside the 1.3% repeat spread.

**Correction.**
- Abstract: "at sampling rates below about 1/20". The Results text already says "about 1/20".
- Results: "with collected evaluation, a kept k-mer costs about 28.5 cycles".

## F7 -- "61 percent of it spurious" is 61 percent of D; of the superset it is 38 percent (C15)

The paper says km-peeler's extended-syncmer method "returned a superset, 61 percent of it spurious". The source says "1,986,807 spurious k-mers (61% of D)" (`results/kp10_vs_shibuya_2026-11-11.md:105`). The superset is 3,238,094 true differences plus 1,986,807 spurious ones, 5,224,901 in total, so the spurious share is **38.0%**. The figure 1,986,807 / 3,238,101 = 61.4% is relative to the true difference. The sentence installs a larger number than the data support.

Also, "equal reliability" sets KP-10's 11.3 B/diff, which already includes a second round (1,275 of 65,536 buckets split and resent, table line "rounds 2"), against a one-shot IBLT at 30/30. The byte counts are honest (36.5 MB = 65,536 × 67 × 8 plus the resent buckets), but the protocols differ in rounds, and the sentence does not say so.

**Correction.** "returned a superset in which 38 percent of the k-mers, 1,986,807, were spurious". Add "(KP-10 in two rounds, the IBLT in one)" to the 0.55× / 1.8× comparison.

## F8 -- the fourth-root hash needs 7 to 8 multiplications against 3, not "about twice as many" (C4)

The congruences hold: 2⁶¹−1 ≡ 1 (mod 3) and ≡ 3 (mod 4). But the measured fourth-root hash is KP7_XREAL (`bench_e.cpp:130`, `kp7.c:49 fmul_real`). Multiplying by a real x costs 2 F_q multiplications per strand, so 4 for both strands, plus a full F_{q²} product for S (3 with Karatsuba, 4 schoolbook). That is **7–8 against the Eisenstein 3**, as `results/kp7_eisenstein_2026-11-10.md:21` itself says ("4 plus a full F_{q²} product"). The cycle ratio is 47.37 / 24.69 = 1.92, and the instruction ratio is 179.3 / 78.5 = 2.28.

**Correction.** "…so that hash lives in F_{q²}, needs 7 or 8 multiplications per k-mer against 3, and costs about twice the cycles."

## F9 -- F·R collapses further than reversal explains, and the table caption names the wrong field (C2)

Autocorrelation identifies a k-mer with its reverse, complement and reverse complement. Burnside's lemma gives 4,160 orbits at k = 7, which matches, but **16,512** at k = 8, while the table prints 16,362. `fr.py` finds **146 further collision groups**, the same at three independent draws of x, so they are polynomial identities and not chance. One example is CAAATCCC against AAACCCCT, which are homometric sequences. Under the Eisenstein encoding F·R gives 16,500 at k = 8 (`exh.py`), so the number of collapses depends on the encoding.

The caption says "fourth-root encoding, q = 2³¹ − 1". But i ∉ F_q for that q (q ≡ 3 mod 4), and `prototype/kp7_symmetric_check.py:1` computes over F_{q²}. Everywhere else the paper uses q = 2⁶¹ − 1.

The F+R formula and the 9^⌊k/2⌋ bound are correct. I reproduced 729 and 6,561 for both encodings.

**Correction.**
- Text: "F·R is an autocorrelation, so it cannot tell a k-mer from its reverse or from other homometric k-mers".
- Caption: "fourth-root encoding over F_{p²}, p = 2³¹ − 1".

## F10 -- ntHash's canonical hash is not GF(2)-linear (C18)

The Discussion says the shared properties "hold for ntHash too, which is linear over GF(2)". The strand hashes are GF(2)-linear (XOR of bit-permuted seeds). But the canonical value that ntHash outputs, and that the paper benchmarks, is `f + r` in integer arithmetic mod 2⁶⁴ (`hash/ntbatch.h:75,80`), and that is not GF(2)-linear. The multi-hash extension, multiply then xorshift (`stats_seeds.cpp:111`), is non-linear too. The properties still hold, but only by keeping the (f, r) pair, the same way this hash needs (F, R).

**Correction.** "…which hold for ntHash's strand hashes, each linear over GF(2), kept as a pair".

## F11 -- the bit-for-bit check is not in the paper's job, and it is a checksum plus a 2-million-k-mer prefix (C12)

The paper's job runs `bench_e`. Its identity pairs (`bench_e.cpp:164-165`) do not include batched ntHash against the library. The check lives in `bench_pack.cpp:276-296`, a different binary. It compares the XOR of all hashes plus the count genome-wide, and the element-wise (position, hash) pairs only for the first 2,000,000 k-mers of the first record, about 0.27% of the input. "About twice" holds (23.53 / 11.37 = 2.07).

**Correction.** "…reproduces the library's hashes (identical XOR checksum and count on all of cs10, element for element on the first 2 million k-mers)".

## F12 -- the upper end of the Ts/Tv range comes from a repeat-masked subset (C9)

The source table (`results/ideas_eisenstein_2026-11-10.md:21-24`) shows three full alignments at Ts/Tv 1.585, 1.584 and 1.669 (61.3–62.5% transitions). The value 1.847 (64.9%) is "uppercase only (outside soft-masked repeats)" for the haplotype pair alone. The conversion r/(1+r) itself is correct: 61.2% and 64.9%.

**Correction.** "61 to 63 percent of substitutions are transitions (Ts/Tv 1.58 to 1.67), and 65 percent outside soft-masked repeats".

## Claims checked and held (no separate finding)

- **C3.** Derived both recurrences. R_{j+1} = xR_j + s_j x^k − s_{j+k} holds. Single-strand bound (k−1)/q is correct.
- **C5.** The identity and the "exactly when" are correct. A fourth-root map with complement as negation forces C = ±iA, so Re(AC̄) = 0. I did not search for priority.
- **C8.** |1+ω|² = 1, |1−ω|² = 3, |2|² = 4. The sentence disclaims distance in F_q.
- **C10.** Recomputed every ratio from `results/paper_bench_2026-11-11.md`: 0.46, 0.48, 0.41, 0.24, 0.43, 0.44, 1.21, 1.14, 0.97, 0.67. All match, and every spread is ≤1.3%.
- **C13.** N·r·(1−r)²⁰⁰ with r = 2/21 and N = 736,547,799 gives 0.142. The long gaps sum to 341+363+371+332+367 = 1,774. The Eisenstein hash has 0 outside repeats.
- **C17.** 2⁸ × 67 × 8 = 137,216 B, and 17.3 B/diff, one round, as stated.
- **C20.** Rebuilt in scratch: `make test_e` passes. The ROLL, LANE and TAB mutants are killed (24, 12 and 64,072 failed checks).

| claim | C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 | C9 | C10 | C11 | C12 | C13 | C14 | C15 | C16 | C17 | C18 | C19 | C20 | C21 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| verdict | NARROWED | NARROWED | SURVIVED | NARROWED | SURVIVED | UNCHALLENGED | UNCHALLENGED | SURVIVED | NARROWED | SURVIVED | NARROWED | NARROWED | SURVIVED | NARROWED | NARROWED | REFUTED | SURVIVED | NARROWED | UNCHALLENGED | SURVIVED | UNCHALLENGED |
