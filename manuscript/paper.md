---
title: "Provably strand-symmetric polynomial hashing of k-mers"
author:
  - Richard Burhans (ORCID 0009-0006-5949-0958)
date: "Draft, 2026-11-11"
link-citations: true
---

<!-- TODO(user): co-authors, affiliation, AI-use disclosure, target venue, licence, and when the repository goes public. -->

# Abstract

<!-- ABSTRACT:BEGIN -->
**Provisional abstract, pending the KP-10 sketch comparison.** Comparing k-mer content across genomes needs a hash that gives a k-mer and its reverse complement one value. The widely used rolling hashes ntHash and ntHash2 do this with fixed seeds, so their collision behaviour is measured, not proven. We read a k-mer as a polynomial over the prime field $\mathbb{F}_q$, $q = 2^{61} - 1$, with bases mapped to sixth roots of unity—$A = 1$, $C = \omega$, $G = -\omega$, $T = -1$, the $\theta = \pi/3$ member of the trigonometric encoding of Das, Nanda and Das. The cube root of unity $\omega$ lies in $\mathbb{F}_q$ because $q \equiv 1 \pmod 3$. Complement becomes negation and reversal a reflection of the polynomial, so the product $(F + \lambda)(R + \lambda)$ of the forward and reverse-complement evaluations is the same for both strands. Its unique factorization in $\lambda$ separates every pair of distinct double-stranded k-mers: two collide with probability at most about $2k/q$, $2.7 \times 10^{-17}$ at $k = 31$, for any input. The obvious alternatives fail, and on the 8,192 double-stranded 7-mers the sum $F + R$ takes 729 values and the product $F \cdot R$ takes 4,160. The same encoding makes a transition the cheapest substitution, at squared distance 1 against 3 for the A–C and G–T transversions. On a 2.0 GHz Xeon, a rolling implementation hashes every k-mer of the Cannabis sativa cs10 reference at 24.6 cycles per k-mer, 0.46 times the speed of a batched ntHash. A sampler that computes the provable hash only for the k-mers a packed-code filter keeps runs faster than batched ntHash at sampling rates of 1/19 and below, 1.21 times as fast at 1/1000.
<!-- ABSTRACT:END -->

# Introduction

<!-- PENDING: KP-10 vs Shibuya. Intent: why genome comparison needs a strand-independent k-mer hash, what ntHash and ntHash2 guarantee (fixed seeds, fwd + rev, empirical quality), and what a proof adds: a per-pair bound for any input, exact strand handling, and a hash native to algebraic set sketches. The framing depends on whether the sketch comparison shows an advantage. -->

# Methods

## A complex encoding in which complement is negation

We map the four bases to the sixth roots of unity $A = 1$, $C = \omega$, $G = -\omega$, $T = -1$, where $\omega = e^{2\pi i/3}$. Complementary bases are negatives of each other. Up to a rotation of the complex plane, this is the trigonometric mapping of Das, Nanda and Das [@das_2017_trigonometric-mapping], who fix its angle at $\theta = \pi/3$ and motivate it by the shared imaginary parts of the purines and pyrimidines and by complements summing to zero. Our earlier fourth-root map ($A = 1$, $C = i$, $G = -i$, $T = -1$) is the $\theta = \pi/4$ member of the same family [@das_2021_dna-encoding]. Anastassiou's complex representation makes complement a conjugation, so it differs [@anastassiou_2001_genomic-signal-processing]. The real-valued mapping in Mendizabal and colleagues also makes complement a negation, and their Table 6 shows that its spectra cannot tell a strand from its complement [@mendizabal_2017_dna-numerical].

A k-mer $s_0 s_1 \ldots s_{k-1}$ is the polynomial $p(x) = \sum_t s_t x^t$, and because complement is negation, its complement is $-p(x)$, its reverse is $x^{k-1} p(1/x)$ and its reverse complement is $-x^{k-1} p(1/x)$.

## Which encodings separate transitions from transversions

With complement as negation, $G = -C$ and $T = -A$, so the two transitions sit at the same distance, $|A - G| = |C - T| = |A + C|$, and both A–C/G–T transversions at $|A - C|$. For unit-modulus values, $|A + C|^2 - |A - C|^2 = 4\,\mathrm{Re}(A \bar C)$, so transitions separate from those transversions exactly when $\mathrm{Re}(A \bar C) \neq 0$. A fourth-root map has $A \perp C$, so both distances equal 2 and it cannot separate them. With $A = 1$ and $C = \omega$, $\mathrm{Re}(A \bar C) = -1/2$, so the squared distances are 1 for transitions, 3 for the A–C and G–T transversions, and 4 for complements. A transition, the commoner substitution, is the cheapest. In pairwise alignments of cannabis assemblies on chromosome NC_044378.1, 61 to 65 percent of substitutions are transitions (Ts/Tv 1.58 to 1.85). This metric lives in the complex embedding, and the field hash below keeps the algebra, not the distances.

## A polynomial hash over a prime field

We evaluate $p$ at a random point $x \in \mathbb{F}_q$ with $q = 2^{61} - 1$. A Mersenne prime makes reduction cheap. Because $q \equiv 1 \pmod 3$, a primitive cube root of unity exists in $\mathbb{F}_q$, so all four base values are field elements and each hash is one 61-bit element. The fourth-root encoding needs $i$, which $\mathbb{F}_q$ lacks since $q \equiv 3 \pmod 4$, so that hash lives in $\mathbb{F}_{q^2}$ and costs about twice as many multiplications. Write $F = p(x)$ for the forward hash and $R = -x^{k-1} p(1/x)$ for the reverse-complement hash. Both roll in constant time per base:

$$F_{j+1} = (F_j - s_j)\,x^{-1} + s_{j+k}\,x^{k-1}, \qquad R_{j+1} = x\,(R_j + s_j\,x^{k-1}) - s_{j+k}.$$

Distinct k-mers give distinct coefficient vectors, so $F(s) = F(u)$ for $s \ne u$ holds for at most $k - 1$ values of $x$, a probability of at most $(k-1)/q$.

## A strand-symmetric product with a proof

A strand-symmetric hash must give $s$ and its reverse complement one value. Reverse complementation swaps $F$ and $R$, so any symmetric function of the pair qualifies, but the two elementary ones fail. The sum $F + R = \sum_t s_t (x^t - x^{k-1-t})$ depends only on the mirror differences $s_t - s_{k-1-t}$, and the middle base cancels when $k$ is odd. The product $F \cdot R = -x^{k-1} p(x) p(1/x)$ is an autocorrelation, so it cannot tell a k-mer from its plain reverse. Under the fourth-root encoding, exhaustive enumeration confirms both collapses: on all 7-mers, $F + R$ takes 729 values and $F \cdot R$ takes 4,160, against 8,192 double-stranded classes.

We use their combination with a random $\lambda \in \mathbb{F}_q$,

$$S = (F + \lambda)(R + \lambda) = \lambda^2 + \lambda (F + R) + F R.$$

As a polynomial in $\lambda$, $S$ factors uniquely into $(\lambda + F)(\lambda + R)$. It therefore determines the unordered pair $\{F, R\}$, which determines the double-stranded k-mer. For two distinct classes the difference $S(s) - S(u)$ is a nonzero polynomial in $(x, \lambda)$ of total degree at most $2(k-1)$, so by the Schwartz–Zippel lemma [@schwartz_1980_polynomial-identities; @zippel_1979_sparse-polynomials] they collide with probability at most $2(k-1)/q$ over the random choice of $x$ and $\lambda$. At $k = 31$ that is about $2.7 \times 10^{-17}$ per pair. The bound holds for every input, repetitive or adversarial, as long as $x$ and $\lambda$ are drawn independently of it. Exhaustive enumeration agrees: at $k = 7$ and $k = 8$ the Eisenstein $S$ takes exactly 8,192 and 32,896 values, one per double-stranded class.

## Implementation

The reference roller follows the two update rules above. The fast rollers keep $G = F + \lambda$ and $H = R + \lambda$ as the rolling state, so the constant $\lambda$ never needs adding per k-mer. They merge the two per-strand table lookups into one 16-entry table indexed by the outgoing and incoming base, and leave additions partly reduced until the next multiplication. A batch loop hashes a run of k-mers with the state in registers, and an interleaved version rolls two contiguous stretches of a block together so their dependency chains overlap. The Eisenstein roller needs three modular multiplications per k-mer, one per strand and one for $S$.

The same construction over $\mathbb{F}_{p^2}$ with $p = 2^{31} - 1$ fits 32-bit vector lanes. An AVX2 version runs 8 lanes and an AVX-512 version 16, using Mersenne-31 vector reduction and register shuffles for table lookups. With $x$ drawn from $\mathbb{F}_{p^2}$, its collision bound is about $1.4 \times 10^{-17}$ per pair at $k = 31$.

Sketches need hashes only for a sample of k-mers. The sampler rolls the 2-bit packed codes of both strands with shifts and masks, which is exact for $k \le 32$. A k-mer is kept when $\mathrm{mix}(f) + \mathrm{mix}(r) < t$, where $f$ and $r$ are the forward and reverse-complement codes and mix is a one-multiply 64-bit scrambler. This filter is strand-symmetric, and its collisions are heuristic. Only kept k-mers get the provable hash, evaluated without multiplications from a 16 KB table of precomputed four-base chunk values, followed by one product for $S$. Kept k-mers are either evaluated inline after a branch or collected and evaluated in a tight loop. The branch is better at sparse sampling, and collection at dense sampling.

## Testing and measurement

Every variant is checked against the reference roller, for hashes and k-mer positions, across windows containing N, block sizes from 1 to 4,096 and both strands. Rolling must equal direct evaluation at every window, and $S(s)$ must equal $S(\mathrm{rc}(s))$. Deliberately broken versions (a wrong sign on the incoming base, a shifted lane remainder, $\lambda$ missing from the table) must fail these tests, and they do.

Timings come from hardware performance counters on an Intel Xeon Gold 6330 at 2.0 GHz, single-threaded, counting user-space cycles and instructions through `perf_event_open`. Each method ran in its own counted process over the 736,547,799 31-mers of the cs10 reference. The cost of loading the sequence was measured in the same job and subtracted. The baseline is a batched reimplementation of ntHash's canonical roll that reproduces the library's hashes bit for bit [@mohamadi_2016_nthash; @kazemi_2022_nthash2]. It runs about twice as fast as the library's per-k-mer interface, so comparing against the library would overstate our speed. Both were compiled with `-O3 -march=native`.

# Results

## The product is exact where the obvious symmetric functions collapse

| k | double-stranded classes | $F + R$ | $F \cdot R$ | $(F + \lambda)(R + \lambda)$ |
|---|---|---|---|---|
| 7 | 8,192 | 729 | 4,160 | 8,192 |
| 8 | 32,896 | 6,561 | 16,362 | 32,896 |

: Distinct values of three strand-symmetric functions over all k-mers, fourth-root encoding, $q = 2^{31} - 1$. The sum has at most $9^{\lfloor k/2 \rfloor}$ values, as predicted from the mirror differences.

The product with $\lambda$ is exact at both sizes, as the factorization argument requires. With the Eisenstein encoding it is also exact, 8,192 and 32,896 values.

## Speed against batched ntHash

| method | cycles per k-mer | relative to batched ntHash |
|---|---|---|
| ntHash, batched | 11.41 | 1.00 |
| ntHash, library interface | 23.44 | 0.49 |
| Eisenstein, every k-mer, two lanes | 24.59 | 0.46 |
| Eisenstein, every k-mer, one lane | 27.95 | 0.41 |
| fourth-root over $\mathbb{F}_{q^2}$, every k-mer, two lanes | 47.51 | 0.24 |
| $\mathbb{F}_{p^2}$ with AVX-512, every k-mer, 16 lanes | 25.83 | 0.44 |
| Eisenstein sampler at 1/1000 | 9.41 | 1.21 |
| Eisenstein sampler at 1/100 | 9.97 | 1.14 |
| Eisenstein sampler at 1/16 | 11.87 | 0.96 |
| Eisenstein sampler at 1/4 | 17.21 | 0.66 |

: Hashing cost on cs10, $k = 31$, one laila job (7fc2a700073c1f4f). The two sparse sampler rows use inline evaluation and the two dense rows use collected evaluation.

The provable hash for every k-mer costs 24.6 cycles with the Eisenstein encoding, half the cost of the fourth-root version, and matches the 16-lane AVX-512 code without any vector instructions. A kept k-mer costs about 28.5 cycles in the sampler, so the sampler is faster than batched ntHash at sampling rates of about 1/19 and below. FracMinHash sketches commonly sample far more sparsely than that, and at those rates every kept k-mer carries the provable hash.

## Sampling quality matches ntHash

Across the whole cs10 reference and five seeds per method, minimizer density stayed within 0.4 percent of the expected $2/(w+1)$ at $w = 10$, 20 and 50 for ntHash, the packed-code filter and the provable hash, and the largest per-bit deviation from one half averaged 0.0009 to 0.0010. Every method left tens of threshold-sampling gaps longer than 200 k-mers per seed, against 0.14 expected by chance. All but 2 of about 1,450 such gaps lay in soft-masked repeats, so no method produced more long gaps outside repeats than another. These statistics were measured with the fourth-root provable hash and the shared filter, and the Eisenstein provable hash has not yet had its own multi-seed run.

## Spaced seeds stay strand-consistent for any mask

Applying a 0/1 mask to the coefficients gives a spaced-seed version of the same hash, which rolls with one multiplication per strand per mask. Because $(F_m + \lambda)(R_m + \lambda)$ is strand-symmetric for any mask, the hash of a masked window equals that of its reverse complement. ntHash2's spaced-seed hash warns that reverse-complement hashing becomes inconsistent for a non-palindromic mask. On a test with 1 percent substitutions and a non-palindromic mask, our hash kept 0.747 of windows unchanged, the exact double-stranded value, against 0.818 for ntHash2's. The Eisenstein spaced-seed hash ran at 0.78 times the speed of ntHash2's with one mask.

## Algebraic set sketches

<!-- PENDING: KP-10 vs Shibuya. Intent: power-sum sketches of the strand-symmetric hash, linear and mergeable, decoded exactly by Newton's identities, compared on real near-identical cannabis pairs with Shibuya, Belazzougui and Kucherov's km-peeler: bytes per difference, decode failure, full against sampled recovery, and strand handling. -->

# Discussion

<!-- DRAFT: the comparison below is settled. The framing paragraph waits for KP-10. -->

What separates this hash from ntHash and ntHash2 is what can be proven about it. Its collision bound holds for every pair of inputs once the seeds are random, while ntHash fixes its seeds, so its collisions are the same in every run and its quality is measured empirically [@mohamadi_2016_nthash; @kazemi_2022_nthash2]. Its strand symmetry is exact by factorization, where ntHash2's sum of the two strand hashes was exact on all double-stranded 8-mers only in our test. Its values live in a field, so set operations and reconciliation sketches can run on the hash values with one error analysis [@minsky_2003_set-reconciliation].

Several properties are shared and are not advantages. Constant-time rolling, constant-time single-base substitution and hashing of concatenations from their parts all hold for ntHash too, which is linear over GF(2). Strand symmetry by combining the two strand hashes is ntHash2's idea. Long sampling gaps are not a weakness of either. Marçais and colleagues show that folding a k-mer with its reverse complement negates the window guarantee of syncmer-style schemes, but it has no visible effect on fractional sampling [@marcais_2024_knonical], which is how a sketch would use this hash.

The distances that favour transitions live in the complex embedding. The field hash keeps the algebra of the encoding—complement as negation and $\omega^2 + \omega + 1 = 0$—but $\mathbb{F}_q$ has no absolute value, so hash values do not measure substitution distance.

# Supplement S1: ideas that did not pay off

The same algebra suggested several other uses. We tested each on cs10 and record the outcome.

- **De Bruijn adjacency from hashes.** An edge $u \to v$ satisfies $x\,v(x) - u(x) = c\,x^k - u_0$, which gives a constant-time adjacency test between hashes. It was exact on 13.6 million edges of a cs10 chromosome's graph but navigated four to six times slower than packed 2-bit codes, so it pays only for hash-only graphs with $k > 32$.
- **A DFT lower bound for mismatch search.** Windows of k-mer length behave like white noise [@voss_1992_fractal-dna], so eight low-frequency coefficients carry their chance share of the distance. On candidates sharing a seed, the bound pruned 21 percent of pairs at three mismatches where 77 percent could be pruned [@faloutsos_1994_subsequence-matching].
- **FFT correlation as a diagonal finder.** On KOMPa chromosome 7 against cs10, correlation peaks reached 97 to 99 percent precision but recovered under two thirds of aligned bases, because short insertions and deletions break alignments into short diagonals.
- **Reverse-complement palindromes.** A window is a palindrome when $F = R$. The test agreed exactly with packed codes on NC_044378.1 (7,760 palindromic 20-mers) and is faster only for arms longer than 32 bases.
- **Reconciling haplotypes.** Two haplotypes of one plant differ in 42 percent of their chromosome 7 31-mers, so a reconciliation sketch is as large as the set.
- **Counting k-mers within d mismatches.** A transform over $\mathbb{Z}_4^k$ counts, for every k-mer at once, how many genome k-mers lie within $d$ mismatches, 1,370 times faster than direct enumeration at $k = 12$, $d = 3$. The technique is probably known from the Krawtchouk structure of the Hamming scheme, so we report it without a claim of novelty.

# References
