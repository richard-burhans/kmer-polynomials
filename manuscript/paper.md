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
Comparing k-mer content across genomes needs a hash that gives a k-mer and its reverse complement one value. A k-mer can be canonicalized by comparing it with its reverse complement and then hashed, but that costs a comparison of k bases per k-mer, fails for spaced seeds, and gives sketches values with no algebra. The rolling hashes ntHash and ntHash2 provide it with fixed seeds and no proven collision bound. We read a k-mer as a polynomial over $\mathbb{F}_q$, $q = 2^{61} - 1$, with bases at sixth roots of unity, $A = 1$, $C = \omega$, $G = -\omega$, $T = -1$, a published encoding that fits $\mathbb{F}_q$ because $q \equiv 1 \pmod 3$. Complement becomes negation, so the product $(F + \lambda)(R + \lambda)$ of the two strand evaluations is the same for both strands. Distinct double-stranded k-mers collide with probability at most $2(k-1)/q$, $2.6 \times 10^{-17}$ at $k = 31$, for every input, while the sum $F + R$ takes only 729 values on 8,192 classes of 7-mers. Power sums of the hashes form a linear set sketch. Between two cannabis assemblies it needed 11.3 bytes per difference, against 20.5 for an invertible Bloom lookup table. On k-mers as written, it recovered every difference between an assembly and its differently oriented draft, at $k = 31$ and $k = 63$. The hash costs 24.7 cycles per k-mer, 0.46 times the speed of a batched ntHash.
<!-- ABSTRACT:END -->

# Introduction

Comparing genomes by their k-mer content, whether by sketching, containment or set reconciliation, needs a hash that gives a k-mer and its reverse complement one value. A sequence can be read from either strand, and assemblies write their contigs in arbitrary orientations. For $k \le 32$ the usual answer is exact and cheap: pack the k-mer into a 64-bit word, keep the smaller of it and its reverse complement, and hash that. Longer k-mers need multi-word codes or string comparisons per k-mer, a spaced seed has no canonical form that its mirror shares, and a sketch built from rolling hash values never sees the k-mer at all. These cases are better served by symmetry inside the hash, which rolls in constant time per base for any k. ntHash rolls such a hash in constant time per base, first by taking the smaller of the two strand hashes and in ntHash2 by adding them, and ntHash2 extends it to spaced seeds [@mohamadi_2016_nthash; @kazemi_2022_nthash2]. Both fix their seeds, so a collision in one run occurs in every run, and their quality is established by measurement.

We ask what can be proven instead. Read as a polynomial over a prime field and evaluated at a random point, a k-mer collides with any other with probability at most $(k-1)/q$, the classical Karp–Rabin argument [@benbassat_2014_incremental-hashing; @chen_2026_hash-functions]. If the bases are encoded so that complement is negation, reverse complementation becomes an algebraic operation on that polynomial. The product $(F + \lambda)(R + \lambda)$ of the two strand hashes is then a characteristic-polynomial fingerprint of the pair, a construction used in set reconciliation [@minsky_2003_set-reconciliation], and here it is strand-symmetric with a collision bound that holds for every input. With bases at sixth roots of unity, every value lies in $\mathbb{F}_q$ for $q = 2^{61} - 1$ and each hash is one 61-bit field element.

Field values also add. Power sums of the hashes form a linear sketch from which a set difference can be recovered exactly [@dodis_2008_fuzzy-extractors; @minsky_2003_set-reconciliation]. When the hash is strand-symmetric, so is the sketch, without canonicalizing anything first. Between a cannabis assembly and its draft, the sketch of k-mers as written recovered the 7,931 differing k-mers exactly, where a forward-only polynomial hash saw 398 million.

The encoding itself is published [@cristea_2002_genomic-signals; @das_2017_trigonometric-mapping], and so is the idea of an encoding whose distances favour transitions [@yu_2018_encoding-survey]. We contribute four things. First, a criterion for complex maps in which complement is negation: by the identity $|A + C|^2 - |A - C|^2 = 4\,\mathrm{Re}(A \bar C)$, transitions separate from the A–C and G–T transversions exactly when $\mathrm{Re}(A \bar C) \neq 0$, which no fourth-root map satisfies. Second, the use of these values as polynomial hash coefficients over $\mathbb{F}_q$, with the fingerprint above as a strand-symmetric hash. OpenAlex and web searches on 10 and 11 November 2026, and a recent review of nucleotide hashing [@chen_2026_hash-functions], found neither the criterion nor a strand-symmetric use of the fingerprint. Third, rolling and sampling implementations measured against a batched ntHash. Fourth, set sketches of these hashes compared with invertible Bloom lookup tables on real cannabis assemblies.

# Methods

## A complex encoding in which complement is negation

We map the four bases to the sixth roots of unity $A = 1$, $C = \omega$, $G = -\omega$, $T = -1$, where $\omega = e^{2\pi i/3}$. Complementary bases are negatives of each other. Cristea's complex representations already made complement a negation: one is the fourth-root map $A = 1$, $C = i$, $G = -i$, $T = -1$ up to a global sign, and another is the $\theta = \pi/4$ point of what was later called the trigonometric family [@cristea_2002_genomic-signals; @das_2021_dna-encoding]. Das, Nanda and Das fixed the family's angle at $\theta = \pi/3$ [@das_2017_trigonometric-mapping], which is our map up to a rotation of the complex plane. Their rationale is that purines share imaginary parts, as do pyrimidines, and that complements sum to zero. It does not mention distances, and the angle is asserted rather than compared with others. Their exon-prediction results need care: on the main test gene the printed sensitivity is the value of their specificity formula, and the text and two tables give three different predicted exon ranges. A later comparison credits the mapping to the same authors [@gulocak_2022_mapping-effect]. Anastassiou's complex representation makes complement a conjugation, so it differs [@anastassiou_2001_genomic-signal-processing]. The real-valued mapping in Mendizabal and colleagues also makes complement a negation, and their Table 6 shows that its spectra cannot tell a strand from its complement [@mendizabal_2017_dna-numerical].

A k-mer $s_0 s_1 \ldots s_{k-1}$ is the polynomial $p(x) = \sum_t s_t x^t$, and because complement is negation, its complement is $-p(x)$, its reverse is $x^{k-1} p(1/x)$ and its reverse complement is $-x^{k-1} p(1/x)$.

## Which encodings separate transitions from transversions

With complement as negation, $G = -C$ and $T = -A$, so the two transitions sit at the same distance, $|A - G| = |C - T| = |A + C|$, and both A–C/G–T transversions at $|A - C|$. For any complex $A$ and $C$, $|A + C|^2 - |A - C|^2 = 4\,\mathrm{Re}(A \bar C)$, so transitions separate from those transversions exactly when $\mathrm{Re}(A \bar C) \neq 0$. A fourth-root map has $A \perp C$, so both distances equal 2 and it cannot separate them. With $A = 1$ and $C = \omega$, $\mathrm{Re}(A \bar C) = -1/2$, so the squared distances are 1 for transitions, 3 for the A–C and G–T transversions, and 4 for complements. A transition, the commoner substitution, is the cheapest. Encodings whose distances reflect transitions are not new: binary codes were designed so that Hamming or Euclidean distance separates transitions from transversions [@yu_2018_encoding-survey]. Our contribution here is narrower, the criterion for complex maps in which complement is negation. In alignments of cannabis assemblies on chromosome NC_044378.1, 61 to 63 percent of substitutions are transitions (Ts/Tv 1.58 to 1.67), and 65 percent between two haplotypes of one plant outside repeats. This metric lives in the complex embedding, and the field hash below keeps the algebra, not the distances.

## A polynomial hash over a prime field

We evaluate $p$ at a random point $x \in \mathbb{F}_q$ with $q = 2^{61} - 1$. A Mersenne prime makes reduction cheap. Because $q \equiv 1 \pmod 3$, a primitive cube root of unity exists in $\mathbb{F}_q$, so all four base values are field elements and each hash is one 61-bit element. The fourth-root encoding needs $i$, which $\mathbb{F}_q$ lacks since $q \equiv 3 \pmod 4$, so that hash lives in $\mathbb{F}_{q^2}$ and costs 8 multiplications per k-mer against 3. Write $F = p(x)$ for the forward hash and $R = -x^{k-1} p(1/x)$ for the reverse-complement hash. Both roll in constant time per base:

$$F_{j+1} = (F_j - s_j)\,x^{-1} + s_{j+k}\,x^{k-1}, \qquad R_{j+1} = x\,(R_j + s_j\,x^{k-1}) - s_{j+k}.$$

Distinct k-mers give distinct coefficient vectors, so $F(s) = F(u)$ for $s \ne u$ holds for at most $k - 1$ values of $x$, a probability of at most $(k-1)/q$.

## A strand-symmetric product with a proof

A strand-symmetric hash must give $s$ and its reverse complement one value. Reverse complementation swaps $F$ and $R$, so any symmetric function of the pair qualifies, but the two elementary ones fail. The sum $F + R = \sum_t s_t (x^t - x^{k-1-t})$ depends only on the mirror differences $s_t - s_{k-1-t}$, and the middle base cancels when $k$ is odd. The product $F \cdot R = -x^{k-1} p(x) p(1/x)$ is an autocorrelation, so it cannot tell a k-mer from its plain reverse, and k-mers that share an autocorrelation merge too. Under the fourth-root encoding, exhaustive enumeration confirms both collapses: on all 7-mers, $F + R$ takes 729 values and $F \cdot R$ takes 4,160, against 8,192 double-stranded classes.

We use their combination with a random $\lambda \in \mathbb{F}_q$,

$$S = (F + \lambda)(R + \lambda) = \lambda^2 + \lambda (F + R) + F R.$$

As polynomials, $S$ factors uniquely into $(\lambda + F)(\lambda + R)$, so it determines the unordered pair of strand polynomials and hence the double-stranded k-mer. For two distinct classes the difference $S(s) - S(u)$ is therefore a nonzero polynomial in $(x, \lambda)$ of total degree at most $2(k-1)$, so by the Schwartz–Zippel lemma [@schwartz_1980_polynomial-identities; @zippel_1979_sparse-polynomials] they collide with probability at most $2(k-1)/q$, $2.6 \times 10^{-17}$ at $k = 31$. The bound holds for every input, repetitive or adversarial, as long as $x$ and $\lambda$ are drawn independently of each other and of the input. It is a per-pair bound: over the 442 million distinct 31-mers of one assembly it allows about 2.5 colliding pairs in expectation. Exhaustive enumeration agrees: at $k = 7$ and $k = 8$ the Eisenstein $S$ takes exactly 8,192 and 32,896 values, one per double-stranded class.

## Set sketches from power sums

A sketch of a k-mer set assigns each strand-symmetric hash $h$ to one of $2^b$ buckets and stores, per bucket, the count and the power sums $\sum h^j$ for $j = 1, \ldots, m + c$ [@dodis_2008_fuzzy-extractors]. The sketch is linear. The sketch of a union of disjoint sets is the sum of their sketches, and subtracting the sketches of two sets leaves the signed power sums of their symmetric difference. Newton's identities and a Padé step recover up to $m$ differences per bucket and say which set each came from [@minsky_2003_set-reconciliation]. The $c$ extra sums check the result, so a bucket holding more than $m$ differences is detected, except with small probability over the random hash. It is then split and resent, at the cost of one more round [@minsky_2002_practical-reconciliation; @naumenko_2019_erlay]. Each entry is one field element of 8 bytes. We implemented the sketch over both fields, the Eisenstein $\mathbb{F}_q$ and the Gaussian $\mathbb{F}_{p^2}$ with $p = 2^{31} - 1$, as a numpy prototype whose hashes match the C code on 20,000 sampled k-mers. The decoder held both sets in one process and used them to name the elements behind recovered differences. The GigaDB runs canonicalized 2-bit codes before hashing and bucketed by a mix of them. The draft-pair runs hashed k-mers as written, at $k = 31$ and $k = 63$, and bucketed by bits of $S$, so only the hash supplied strand symmetry. In both pairs the sketch configuration was chosen knowing the true difference.

## Implementation

The reference roller follows the two update rules above. The fast rollers keep $G = F + \lambda$ and $H = R + \lambda$ as the rolling state, so the constant $\lambda$ never needs adding per k-mer. They merge the two per-strand table lookups into one 16-entry table indexed by the outgoing and incoming base, and leave additions partly reduced until the next multiplication. A batch loop hashes a run of k-mers with the state in registers, and an interleaved version rolls two contiguous stretches of a block together so their dependency chains overlap. The Eisenstein roller needs three modular multiplications per k-mer, one per strand and one for $S$. Seeds for $x$ and $\lambda$ are drawn independently, from the operating system by default. Our experiments used fixed seeds for reproducibility, so they test the construction rather than the bound, and the GigaDB sketch ran before $x$ and $\lambda$ were separated.

The same construction over $\mathbb{F}_{p^2}$ with $p = 2^{31} - 1$ fits 32-bit vector lanes. An AVX2 version runs 8 lanes and an AVX-512 version 16, using Mersenne-31 vector reduction and register shuffles for table lookups. With $x$ drawn from $\mathbb{F}_{p^2}$, its collision bound is $1.3 \times 10^{-17}$ per pair at $k = 31$.

Sketches need hashes only for a sample of k-mers. The sampler rolls the 2-bit packed codes of both strands with shifts and masks, which is exact for $k \le 32$. A k-mer is kept when $\mathrm{mix}(f) + \mathrm{mix}(r) < t$, where $f$ and $r$ are the forward and reverse-complement codes and mix is a one-multiply 64-bit scrambler. This filter is strand-symmetric but heuristic, and the proof does not cover which k-mers it keeps. Only kept k-mers get the provable hash, evaluated without multiplications from a 16 KB table of precomputed four-base chunk values, followed by one product for $S$. Kept k-mers are either evaluated inline after a branch or collected and evaluated in a tight loop. The branch is better at sparse sampling, and collection at dense sampling.

## Testing and measurement

Every variant is checked against the reference roller, for hashes and k-mer positions, across windows containing N, block sizes of 1, 3, 5, 64 and 4,096 and both strands. Rolling must equal direct evaluation at every window, and $S(s)$ must equal $S(\mathrm{rc}(s))$. Deliberately broken versions (a wrong sign on the incoming base, a shifted lane remainder, $\lambda$ missing from the table) must fail these tests, and they do. Every sketch recovery was compared element by element with the exact set difference. An early version of the recovery read the Padé matrix at negative indices whenever a bucket held more differences on the second set's side, and Python wraps such indices silently to the end of the array. It still recovered every difference, by splitting nearly every bucket, so the bug showed only as extra rounds, 3 instead of 1 on one pair. We fixed it and added one-sided tests, which fail without the fix.

Timings come from hardware performance counters on an Intel Xeon Gold 6330, nominally 2.0 GHz with turbo enabled (2.2 to 3.0 GHz measured), single-threaded, counting user-space core cycles and instructions through `perf_event_open`. Every timing in the speed table comes from one job, in which each method ran in its own counted process over the 736,547,799 31-mers of the cs10 reference, twice, in the same order. The cost of loading the sequence was measured three times in the same job and subtracted. The spaced-seed timing is a single wall-clock run in that job, and the Supplement's timings come from earlier jobs. The baseline is a batched reimplementation of ntHash's canonical roll that reproduces the library's hash for all 736,547,799 k-mers [@mohamadi_2016_nthash; @kazemi_2022_nthash2]. It runs about twice as fast as the library's per-k-mer interface, so comparing against the library would overstate our speed. Both were compiled with `-O3 -march=native`.

# Results

## The product is exact where the obvious symmetric functions collapse

| k | double-stranded classes | $F + R$ | $F \cdot R$ | $(F + \lambda)(R + \lambda)$ |
|---|---|---|---|---|
| 7 | 8,192 | 729 | 4,160 | 8,192 |
| 8 | 32,896 | 6,561 | 16,362 | 32,896 |

: Distinct values of three strand-symmetric functions over all k-mers, fourth-root encoding over $\mathbb{F}_{p^2}$, $p = 2^{31} - 1$. The sum has at most $9^{\lfloor k/2 \rfloor}$ values, as predicted from the mirror differences. Merging only reversals would leave 16,512 values for $F \cdot R$ at $k = 8$.

The product with $\lambda$ is exact at both sizes, as the factorization argument requires, and with the Eisenstein encoding it is also exact, 8,192 and 32,896 values.

## Speed against batched ntHash

| method | cycles per k-mer | relative to batched ntHash |
|---|---|---|
| ntHash, batched | 11.37 | 1.00 |
| ntHash, library interface | 23.53 | 0.48 |
| Eisenstein, every k-mer, two lanes | 24.69 | 0.46 |
| Eisenstein, every k-mer, one lane | 27.82 | 0.41 |
| fourth-root over $\mathbb{F}_{q^2}$, every k-mer, two lanes | 47.37 | 0.24 |
| $\mathbb{F}_{p^2}$ with AVX2, every k-mer, 8 lanes | 26.43 | 0.43 |
| $\mathbb{F}_{p^2}$ with AVX-512, every k-mer, 16 lanes | 25.86 | 0.44 |
| Eisenstein sampler at 1/1000 | 9.40 | 1.21 |
| Eisenstein sampler at 1/100 | 10.01 | 1.14 |
| Eisenstein sampler at 1/16 | 11.74 | 0.97 |
| Eisenstein sampler at 1/4 | 17.08 | 0.67 |

: Hashing cost on cs10, $k = 31$, all from one laila job (08653d9e6a7453a2), the mean of two repeats, which differed by at most 1.3 percent. The two sparse sampler rows use inline evaluation and the two dense rows use collected evaluation.

The provable hash for every k-mer costs 24.7 cycles with the Eisenstein encoding, about half the cycles of the fourth-root version, and matches the 16-lane AVX-512 code without any vector instructions. With collected evaluation a kept k-mer costs about 28.5 cycles, which puts the break-even with batched ntHash at a sampling rate of 1/20.3. The sampler's advantage comes from its cheap filter, and the provable hash adds cost per kept k-mer. sourmash, a widely used FracMinHash tool, keeps 1 in 1,000 k-mers by default [@brown_2016_sourmash; @irber_2022_fracminhash], a rate at which every kept k-mer carries the provable hash and the sampler is 1.21 times as fast as batched ntHash.

## Sampling quality matches ntHash

Across the whole cs10 reference and five seeds per method, minimizer density stayed within 0.6 percent of the expected $2/(w+1)$ at $w = 10$, 20 and 50 for ntHash, two packed-code filters, the fourth-root provable hash and the Eisenstein provable hash. The largest per-bit deviation from one half averaged 0.0009 to 0.0010 for every method, a sanity check that cannot separate hashes below about $10^{-3}$. Every method left 8 to 96 threshold-sampling gaps longer than 200 k-mers per seed, against 0.14 expected by chance. The 1,583 such gaps fell at 151 distinct loci, 62 of them hit by all five methods, and 147 lay in soft-masked repeats. The other 4 loci held one gap each, against about 2 expected by chance in unmasked sequence over 25 runs, so the long gaps reflect repeats rather than any hash.

## Spaced seeds stay strand-consistent for any mask

Applying a 0/1 mask to the coefficients gives a spaced-seed version of the same hash, which rolls with one multiplication per strand per mask. Because $(F_m + \lambda)(R_m + \lambda)$ is strand-symmetric for any mask, the hash of a masked window equals that of its reverse complement. ntHash2's spaced-seed code warns that reverse-complement hashing becomes inconsistent for a non-symmetric mask, and for our non-symmetric mask it disagreed with its own reverse complement at every one of 199,970 windows. Symmetry has a cost. The product depends on the mask and its mirror image together, 29 positions for our 20-position mask, so under 1 percent substitutions it kept 0.747 of windows unchanged ($0.99^{29}$), against 0.818 ($0.99^{20}$) for ntHash2's. The Eisenstein spaced-seed hash ran at 0.75 times the speed of ntHash2's with one mask, over the first 60 million bases of cs10.

## Set sketches recover k-mer differences exactly

We compared the power-sum sketch with km-peeler, the invertible Bloom lookup tables of Shibuya, Belazzougui and Kucherov [@shibuya_2022_genomic-reconciliation; @goodrich_2011_iblt; @eppstein_2011_set-difference], on two pairs of cannabis assemblies. Cannbio-2 and its GigaDB release differ in 3,238,101 canonical 31-mers, almost all from unplaced scaffolds that the release drops. Cannbio-2 and its draft differ in 7,931.

| method | bytes per difference | exact recovery | rounds |
|---|---|---|---|
| power-sum sketch, either field | 11.3 | yes, in every run | 2 |
| km-peeler, k-mer table at 1.05 times the difference | 20.5 | in 30 of 30 trials | 1 |
| km-peeler, extended syncmers | 4.2 | no, 1,986,807 spurious and 7 missed | 1 |

: Recovering the GigaDB difference, $k = 31$. The sketch used $2^{16}$ buckets with $m = 64$ and $c = 2$, and 1,275 buckets were resent in its second round. The k-mer table was sized to 1.05 times the true difference and tested on 4 of 64 partitions, and at exactly the difference it succeeded in 15 of 30 trials.

Both fields gave sketches of the same size, 36.5 MB, because a $\mathbb{F}_{p^2}$ element with $p = 2^{31} - 1$ also takes 8 bytes. Every recovered difference matched the exact set. The k-mer table needed 1.8 times as many bytes, and its failures were probabilistic. Its 30 of 30 successes bound the failure rate only below 11.6 percent, and on the first four partitions of the draft pair, 488 differences, it failed 4 of 30 trials at 1.05 times the difference.

The extended-syncmer method of km-peeler, which its authors propose for recovering every difference, was 2.7 times smaller than our sketch but returned spurious syncmers numbering 61 percent of the true difference on this pair, against about 3 percent in the authors' tests, and missed 7 differences, which we did not trace. It reads each assembly as written. The draft's contigs are written in other orientations, 45 percent of its bases on flipped contigs, so 48,293,748 extended syncmers differed between draft and final against 7,931 canonical k-mers, and the run ran out of memory while listing. Hashing k-mers as written, with buckets keyed by bits of $S$, our sketch recovered the 7,931 exactly in 137 KB and one round. A forward-only polynomial hash in the same sketch saw 398,217,248 differences and overflowed every bucket. At $k = 63$ the hash, computed from the sequence with no packed code, recovered all 21,382 differences in 549 KB and one round, matching the canonical difference computed by KMC [@kokot_2017_counting] element by element.

The sketch is slower to build. It costs about $m + c = 66$ field multiplications per element against 3 cell updates for an invertible Bloom lookup table, and our numpy prototype took about an hour per pair, against seconds to minutes for km-peeler's C code.

# Discussion

For k-mers that fit a machine word, a canonical packed code and a universal hash already give strand symmetry with a proven bound, faster than our hash. The strand-symmetric product earns its place where that recipe is costly or unavailable: for long k-mers, where canonical codes span several words, for spaced seeds and for values that must be added. It gives a collision bound for every pair of inputs, strand symmetry by factorization and a field in which the hashes themselves form a linear set sketch. The sketch recovered every k-mer difference exactly and used 0.55 times the bytes of an invertible Bloom lookup table, at the cost of a second round and about 66 field multiplications per element. Without canonicalization it recovered every difference between an assembly and its draft at $k = 31$ and $k = 63$, where a strand-specific hash saw 398 million at $k = 31$. A table of canonical k-mers is orientation-blind too, so this shows that the hash supplies the symmetry, not that it beats canonicalization.

Against ntHash and ntHash2, the difference lies in what can be proven, not in speed. This hash's collision bound holds for every pair of inputs once the seeds are random, while ntHash fixes its seeds, so its collisions are the same in every run and its quality is measured empirically [@mohamadi_2016_nthash; @kazemi_2022_nthash2]. Its strand symmetry holds for any spaced-seed mask, where ntHash2's does not.

Several properties are shared and are not advantages. Constant-time rolling, constant-time single-base substitution and hashing of concatenations from their parts all hold for ntHash's per-strand hashes, which are linear over GF(2). Combining the two strand hashes is ntHash's idea. Long sampling gaps are not a weakness of either. Marçais and colleagues show that folding a k-mer with its reverse complement can negate the window guarantee of syncmer-style schemes, while fractional sampling, which has no such guarantee, shows long sketching deserts before and after folding [@marcais_2024_knonical]. That is how a sketch would use this hash.

The distances that favour transitions live in the complex embedding. The field hash keeps the algebra of the encoding—complement as negation and $\omega^2 + \omega + 1 = 0$—but $\mathbb{F}_q$ has no absolute value, so hash values do not measure substitution distance.

# Supplement S1: ideas that did not pay off

The same algebra suggested several other uses. We tested each on cs10 and record the outcome.

- **De Bruijn adjacency from hashes.** An edge $u \to v$ satisfies $x\,v(x) - u(x) = c\,x^k - u_0$, which gives a constant-time adjacency test between hashes. It was exact on 13.6 million edges of the graph of the first 20 Mb of cs10 but navigated four to six times slower than packed 2-bit codes, so it is useful mainly for hash-only graphs with $k > 32$.
- **A DFT lower bound for mismatch search.** Apart from the period-3 signal of coding sequence, windows of k-mer length behave like white noise [@voss_1992_fractal-dna], so eight low-frequency coefficients carry their chance share of the distance. On candidates sharing a seed, the bound pruned 21 percent of pairs at three mismatches where 77 percent could be pruned [@faloutsos_1994_subsequence-matching].
- **FFT correlation as a diagonal finder.** This is the FFT step of MAFFT [@katoh_2002_mafft]. On KOMPa chromosome 7 against cs10, correlation peaks reached 97 to 99 percent precision but recovered under two thirds of aligned bases, because short insertions and deletions break alignments into short diagonals.
- **Reverse-complement palindromes.** A window is a palindrome when $F = R$. The test agreed exactly with packed codes on NC_044378.1 (7,760 palindromic 20-mers) and is faster only for arms longer than 32 bases.
- **Reconciling haplotypes.** Two haplotypes of one plant differ in 35 percent of the 31-mers of their chromosome 7, so a reconciliation sketch is as large as the set.
- **Counting k-mers within d mismatches.** A transform over $\mathbb{Z}_4^k$ counts, for every k-mer at once, how many genome k-mers lie within $d$ mismatches, 1,370 times faster than direct enumeration at $k = 12$, $d = 3$ in one run. The technique is probably known from the Krawtchouk structure of the Hamming scheme, so we report it without a claim of novelty.

# References
