---
title: "Provably strand-symmetric polynomial hashing of k-mers"
author:
  - Richard Burhans (ORCID 0009-0006-5949-0958)
date: "Draft, 2026-11-12"
link-citations: true
---

<!-- TODO(user): co-authors, affiliation, AI-use disclosure, target venue, licence, and when the repository goes public. -->

# Abstract

<!-- ABSTRACT:BEGIN -->
Comparing k-mer content across genomes needs a hash that gives a k-mer and its reverse complement one value. For k-mers that fit a machine word, a canonical packed code and a universal hash give one with a proof, but longer k-mers and spaced seeds lack that route, and the rolling hashes ntHash and ntHash2 carry no proven collision bound. We read a k-mer as a polynomial over $\mathbb{F}_q$, $q = 2^{61} - 1$, with bases at sixth roots of unity, a published encoding in which complement is negation. The product $(F + \lambda)(R + \lambda)$ of the two strand evaluations is then the same for both strands, and distinct double-stranded k-mers collide with probability at most $2(k-1)/q$, $2.6 \times 10^{-17}$ at $k = 31$, for every input. Power sums of these hashes form a linear set sketch that needs no canonicalization. Between two cannabis assemblies differing in 3.2 million 31-mers, it recovered every difference in 30 of 30 runs with 11.3 bytes per difference over two rounds, against 20.5 for a one-round invertible Bloom lookup table. Because an overflowing bucket is resent alone, an eightfold underestimate of the difference raised this to 20.8, against 38.5 for the table. Hashing every k-mer takes 2.2 times as long as a batched ntHash, and a sampler that hashes 1 in 1,000 k-mers runs 1.21 times as fast, about as fast as the same sampler with a universal hash.
<!-- ABSTRACT:END -->

# Introduction

Comparing genomes by their k-mer content, whether by sketching, containment or set reconciliation, needs a hash that gives a k-mer and its reverse complement one value. A sequence can be read from either strand, and assemblies write their contigs in arbitrary orientations. For $k \le 32$ the usual answer is exact and cheap: pack the k-mer into a 64-bit word, keep the smaller of it and its reverse complement, and hash that. Longer k-mers need multi-word codes or string comparisons per k-mer, a spaced seed has no canonical form that its mirror shares, and a sketch built from rolling hash values never sees the k-mer at all. These cases are better served by symmetry inside the hash, which rolls in constant time per base for any k. ntHash rolls such a hash in constant time per base, first by taking the smaller of the two strand hashes and in ntHash2 by adding them, and ntHash2 extends it to spaced seeds [@mohamadi_2016_nthash; @kazemi_2022_nthash2]. Both fix their seeds, so a collision in one run occurs in every run, and their quality is established by measurement.

We ask what can be proven instead. Read as a polynomial over a prime field and evaluated at a random point, a k-mer collides with any other with probability at most $(k-1)/q$, the classical Karp–Rabin argument [@benbassat_2014_incremental-hashing; @chen_2026_hash-functions]. If the bases are encoded so that complement is negation, reverse complementation becomes an algebraic operation on that polynomial. The product $(F + \lambda)(R + \lambda)$ of the two strand hashes is then a characteristic-polynomial fingerprint of the pair, a construction used in set reconciliation [@minsky_2003_set-reconciliation], and here it is strand-symmetric with a collision bound that holds for every input. With bases at sixth roots of unity, every value lies in $\mathbb{F}_q$ for $q = 2^{61} - 1$ and each hash is one 61-bit field element.

Field values also add. Power sums of the hashes form a linear sketch from which a set difference can be recovered exactly [@dodis_2008_fuzzy-extractors; @minsky_2003_set-reconciliation]. The sketch spreads the difference over buckets and checks each one, so a bucket that holds too many differences is split and resent on its own in a second round, whereas an invertible Bloom lookup table that fails must be resent whole. When the hash is strand-symmetric, so is the sketch, without canonicalizing anything first. Between a cannabis assembly and its draft, the sketch of k-mers as written recovered the 7,931 differing k-mers exactly, where a forward-only polynomial hash saw 398 million.

The encoding itself is published [@cristea_2002_genomic-signals; @das_2017_trigonometric-mapping], and so is the idea of an encoding whose distances favour transitions [@yu_2018_encoding-survey]. We contribute four things. First, a criterion for complex maps in which complement is negation: by the identity $|A + C|^2 - |A - C|^2 = 4\,\mathrm{Re}(A \bar C)$, transitions separate from the A–C and G–T transversions exactly when $\mathrm{Re}(A \bar C) \neq 0$, which no fourth-root map satisfies. Second, the use of these values as polynomial hash coefficients over $\mathbb{F}_q$, with the fingerprint above as a strand-symmetric hash. OpenAlex and web searches on 10 and 11 November 2026, and a recent review of nucleotide hashing [@chen_2026_hash-functions], found neither the criterion nor a strand-symmetric use of the fingerprint. Third, rolling and sampling implementations measured against a batched ntHash and against a canonical code with a universal hash. Fourth, a two-round set sketch of these hashes compared with invertible Bloom lookup tables on real cannabis assemblies.

# Methods

## A complex encoding in which complement is negation

We map the four bases to the sixth roots of unity $A = 1$, $C = \omega$, $G = -\omega$, $T = -1$, where $\omega = e^{2\pi i/3}$. Complementary bases are negatives of each other (Figure 1A). Cristea's complex representations already made complement a negation: his equation (5) is the fourth-root map $A = 1$, $C = i$, $G = -i$, $T = -1$ up to a global sign, and his equation (4), scaled, is the $\theta = \pi/4$ point of what was later called the trigonometric family [@cristea_2002_genomic-signals]. Das, Nanda and Das fixed the family's angle at $\theta = \pi/3$ [@das_2017_trigonometric-mapping], which is our map up to a rotation of the complex plane. Their rationale restates Cristea's: purines share imaginary parts, as do pyrimidines, and complements sum to zero. It does not mention distances, and the angle is asserted, not compared with others. Their exon-prediction results need care: on the main test gene the printed sensitivity is the value of their specificity formula, and the text and two tables give three different predicted exon ranges. Their 2019 paper keeps $\theta = \pi/3$ and credits the 2017 one, but its printed map swaps A and G, so complement becomes conjugation [@das_2019_adaptive-kaiser-window]. A later comparison cites that version as the mapping's source [@gulocak_2022_mapping-effect]. Anastassiou's complex representation makes complement a conjugation, so it differs [@anastassiou_2001_genomic-signal-processing]. The real-valued mapping in Mendizabal and colleagues also makes complement a negation, and their Table 6 shows that its spectra cannot tell a strand from its complement [@mendizabal_2017_dna-numerical].

A k-mer $s_0 s_1 \ldots s_{k-1}$ is the polynomial $p(x) = \sum_t s_t x^t$, and because complement is negation, its complement is $-p(x)$, its reverse is $x^{k-1} p(1/x)$ and its reverse complement is $-x^{k-1} p(1/x)$.

![Complement is negation, so the product of the two strand hashes is strand-symmetric. (A) The four bases at sixth roots of unity, $\omega = e^{2\pi i/3}$. Complementary bases are reflections through the origin (dashed). Transitions (thick lines) are the shortest chords, squared length 1, against 3 for the A–C and G–T transversions and 4 for complements. (B) Reverse complementation swaps the coefficients of $F$ and $R$, so $(F + \lambda)(R + \lambda)$ takes one value on both strands.](figures/fig1_encoding.pdf)

## Which encodings separate transitions from transversions

With complement as negation, $G = -C$ and $T = -A$, so the two transitions sit at the same distance, $|A - G| = |C - T| = |A + C|$, and both A–C/G–T transversions at $|A - C|$. For any complex $A$ and $C$, $|A + C|^2 - |A - C|^2 = 4\,\mathrm{Re}(A \bar C)$, so transitions separate from those transversions exactly when $\mathrm{Re}(A \bar C) \neq 0$. A fourth-root map has $A \perp C$, so both distances equal 2 and it cannot separate them. With $A = 1$ and $C = \omega$, $\mathrm{Re}(A \bar C) = -1/2$, so the squared distances are 1 for transitions, 3 for the A–C and G–T transversions, and 4 for complements. A transition, the commoner substitution, is the cheapest. Encodings whose distances reflect transitions are not new: binary codes were designed so that Hamming or Euclidean distance separates transitions from transversions [@yu_2018_encoding-survey]. Our contribution here is narrower, the criterion for complex maps in which complement is negation. In alignments of cannabis assemblies on chromosome NC_044378.1, 61 to 63 percent of substitutions are transitions (Ts/Tv 1.58 to 1.67), and 65 percent between two haplotypes of one plant outside repeats. This metric lives in the complex embedding, and the field hash below keeps the algebra, not the distances.

## A polynomial hash over a prime field

We evaluate $p$ at a random point $x \in \mathbb{F}_q$ with $q = 2^{61} - 1$. A Mersenne prime makes reduction cheap. Because $q \equiv 1 \pmod 3$, a primitive cube root of unity exists in $\mathbb{F}_q$, so all four base values are field elements and each hash is one 61-bit element. The fourth-root encoding needs $i$, which $\mathbb{F}_q$ lacks since $q \equiv 3 \pmod 4$, so that hash lives in $\mathbb{F}_{q^2}$ and costs 8 multiplications per k-mer against 3. Write $F = p(x)$ for the forward hash and $R = -x^{k-1} p(1/x)$ for the reverse-complement hash. Both roll in constant time per base:

$$F_{j+1} = (F_j - s_j)\,x^{-1} + s_{j+k}\,x^{k-1}, \qquad R_{j+1} = x\,(R_j + s_j\,x^{k-1}) - s_{j+k}.$$

Distinct k-mers give distinct coefficient vectors, so $F(s) = F(u)$ for $s \ne u$ holds for at most $k - 1$ values of $x$, a probability of at most $(k-1)/q$.

## A strand-symmetric product with a proof

A strand-symmetric hash must give $s$ and its reverse complement one value. Reverse complementation swaps $F$ and $R$ (Figure 1B), so any symmetric function of the pair qualifies, but the two elementary ones fail. The sum $F + R = \sum_t s_t (x^t - x^{k-1-t})$ depends only on the mirror differences $s_t - s_{k-1-t}$, and the middle base cancels when $k$ is odd. The product $F \cdot R = -x^{k-1} p(x) p(1/x)$ is an autocorrelation, so it cannot tell a k-mer from its plain reverse, and k-mers that share an autocorrelation merge too. Under the fourth-root encoding, exhaustive enumeration confirms both collapses: on all 7-mers, $F + R$ takes 729 values and $F \cdot R$ takes 4,160, against 8,192 double-stranded classes (Figure 2).

We use their combination with a random $\lambda \in \mathbb{F}_q$,

$$S = (F + \lambda)(R + \lambda) = \lambda^2 + \lambda (F + R) + F R.$$

As polynomials, $S$ factors uniquely into $(\lambda + F)(\lambda + R)$, so it determines the unordered pair of strand polynomials and hence the double-stranded k-mer. For two distinct classes the difference $S(s) - S(u)$ is therefore a nonzero polynomial in $(x, \lambda)$ of total degree at most $2(k-1)$, so by the Schwartz–Zippel lemma [@schwartz_1980_polynomial-identities; @zippel_1979_sparse-polynomials] they collide with probability at most $2(k-1)/q$, $2.6 \times 10^{-17}$ at $k = 31$. The bound holds for every input, repetitive or adversarial, as long as $x$ and $\lambda$ are drawn independently of each other and of the input. It is a per-pair bound: over the 442 million distinct 31-mers of one assembly it allows about 2.5 colliding pairs in expectation. Exhaustive enumeration agrees: at $k = 7$ and $k = 8$ the Eisenstein $S$ takes exactly 8,192 and 32,896 values, one per double-stranded class.

![Of three symmetric functions of $F$ and $R$, only the product with $\lambda$ separates every double-stranded class. Distinct values over all 16,384 7-mers, fourth-root encoding over $\mathbb{F}_{p^2}$, $p = 2^{31} - 1$. The dashed line marks the 8,192 double-stranded classes.](figures/fig2_collapse.pdf)

## Set sketches from power sums

A sketch of a k-mer set assigns each strand-symmetric hash $h$ to one of $2^b$ buckets and stores, per bucket, the count and the power sums $\sum h^j$ for $j = 1, \ldots, m + c$ [@dodis_2008_fuzzy-extractors]. The sketch is linear. The sketch of a union of disjoint sets is the sum of their sketches, and subtracting the sketches of two sets leaves the signed power sums of their symmetric difference. Newton's identities and a Padé step recover up to $m$ differences per bucket and say which set each came from [@minsky_2003_set-reconciliation]. The $c$ extra sums check the result, so a bucket holding more than $m$ differences is detected, except with small probability over the random hash. It is then split and resent, at the cost of one more round [@minsky_2002_practical-reconciliation; @naumenko_2019_erlay]. Each entry is one field element of 8 bytes, in either field. The whole-genome comparison used a C implementation over the Eisenstein field with $m = 64$ and $c = 2$. It hashes every 31-mer as written, with no canonical code, so only $S$ supplies strand symmetry, and a bucket is a range of $S$. Each party checks its buckets with its own data only and finds its own differing k-mers by rescanning its genome for the recovered values. The first round uses the smallest power of two of buckets that holds at most $0.8m$ estimated differences each. We ran it from the true difference and from estimates between one eighth and four times it. The $k = 63$ run used a numpy prototype, whose hashes match the C code on 20,000 sampled k-mers, with buckets sized from the true difference.

## Implementation

The reference roller follows the two update rules above. The fast rollers keep $G = F + \lambda$ and $H = R + \lambda$ as the rolling state, so the constant $\lambda$ never needs adding per k-mer. They merge the two per-strand table lookups into one 16-entry table indexed by the outgoing and incoming base, and leave additions partly reduced until the next multiplication. A batch loop hashes a run of k-mers with the state in registers, and an interleaved version rolls two contiguous stretches of a block together so their dependency chains overlap. The Eisenstein roller needs three modular multiplications per k-mer, one per strand and one for $S$. Seeds for $x$ and $\lambda$ are drawn independently, from the operating system by default. Our experiments used fixed seeds for reproducibility, so they test the construction rather than the bound, and the GigaDB sketch ran before $x$ and $\lambda$ were separated.

The same construction over $\mathbb{F}_{p^2}$ with $p = 2^{31} - 1$ fits 32-bit vector lanes. An AVX2 version runs 8 lanes and an AVX-512 version 16, using Mersenne-31 vector reduction and register shuffles for table lookups. With $x$ drawn from $\mathbb{F}_{p^2}$, its collision bound is $1.3 \times 10^{-17}$ per pair at $k = 31$.

Sketches need hashes only for a sample of k-mers. The sampler rolls the 2-bit packed codes of both strands with shifts and masks, which is exact for $k \le 32$. A k-mer is kept when $\mathrm{mix}(f) + \mathrm{mix}(r) < t$, where $f$ and $r$ are the forward and reverse-complement codes and mix is a one-multiply 64-bit scrambler. This filter is strand-symmetric but heuristic, and the proof does not cover which k-mers it keeps. Only kept k-mers get the provable hash, evaluated without multiplications from a 16 KB table of precomputed four-base chunk values, followed by one product for $S$. Kept k-mers are either evaluated inline after a branch or collected and evaluated in a tight loop. The branch is better at sparse sampling, and collection at dense sampling.

## Testing and measurement

Every variant is checked against the reference roller, for hashes and k-mer positions, across windows containing N, block sizes of 1, 3, 5, 64 and 4,096 and both strands. Rolling must equal direct evaluation at every window, and $S(s)$ must equal $S(\mathrm{rc}(s))$. Deliberately broken versions (a wrong sign on the incoming base, a shifted lane remainder, $\lambda$ missing from the table) must fail these tests, and they do. Every sketch recovery was compared element by element with the exact set difference. An early version of the recovery read the Padé matrix at negative indices whenever a bucket held more differences on the second set's side, and Python wraps such indices silently to the end of the array. It still recovered every difference, by splitting nearly every bucket, so the bug showed only as extra rounds, 3 instead of 1 on one pair. We fixed it and added one-sided tests, which fail without the fix.

Timings come from hardware performance counters on an Intel Xeon Gold 6330, nominally 2.0 GHz with turbo enabled (2.1 to 3.1 GHz measured), single-threaded, counting user-space core cycles and instructions through `perf_event_open`. Every timing in the main text comes from one job at one commit of our code, in which each method ran in its own counted process over the 736,547,799 31-mers of the cs10 reference, five times, with every other repeat in reverse order. We report medians. The cost of loading the sequence was measured six times per program in the same job and subtracted. The spaced-seed timing is the median of five wall-clock rounds in that job, and the Supplement's timings come from earlier jobs. The baseline is a batched reimplementation of ntHash's canonical roll that reproduces the library's hash for all 736,547,799 k-mers [@mohamadi_2016_nthash; @kazemi_2022_nthash2]. It runs about twice as fast as the library's per-k-mer interface, so comparing against the library would overstate our speed. Every program was compiled in the job with `-O3 -march=native`.

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
| ntHash, library interface | 23.92 | 0.48 |
| Eisenstein, every k-mer, two lanes | 24.67 | 0.46 |
| Eisenstein, every k-mer, one lane | 27.85 | 0.41 |
| fourth-root over $\mathbb{F}_{q^2}$, every k-mer, two lanes | 47.19 | 0.24 |
| $\mathbb{F}_{p^2}$ with AVX2, every k-mer, 8 lanes | 26.44 | 0.43 |
| $\mathbb{F}_{p^2}$ with AVX-512, every k-mer, 16 lanes | 25.94 | 0.44 |
| Eisenstein sampler at 1/1000 | 9.43 | 1.21 |
| Eisenstein sampler at 1/100 | 9.94 | 1.14 |
| Eisenstein sampler at 1/16 | 11.74 | 0.97 |
| Eisenstein sampler at 1/4 | 17.12 | 0.66 |

: Hashing cost on cs10, $k = 31$, from one laila job (c84a395ba84b315e), the median of five repeats. The two sparse sampler rows use inline evaluation and the two dense rows use collected evaluation.

The provable hash for every k-mer costs 24.7 cycles with the Eisenstein encoding, about half the cycles of the fourth-root version, and matches the 16-lane AVX-512 code without any vector instructions. With collected evaluation a kept k-mer costs 28.7 cycles, which puts the break-even with batched ntHash at a sampling rate of 1/20.2 (Figure 3A). The widely used FracMinHash tool sourmash keeps 1 in 1,000 k-mers by default [@brown_2016_sourmash; @irber_2022_fracminhash]. At that rate the sampler is 1.21 times as fast as batched ntHash.

| method | 1/1000 | 1/100 | per kept k-mer |
|---|---|---|---|
| Eisenstein sampler | 1.22 | 1.15 | 22.5 |
| same filter, then ntHash | 1.18 | 1.12 | 19.3 |
| same filter, then a universal hash | 1.25 | 1.20 | 5.9 |
| canonical code and a universal hash, every k-mer | 1.12 | 1.12 | — |

: Speed relative to batched ntHash (11.43 cycles per k-mer in this benchmark program) at two sampling rates, and cycles per kept k-mer from the slope between 1/16 and 1/4, in the same job. The universal hash is multiply-shift on the smaller packed code, exact for $k \le 32$.

The sampler's advantage comes from its filter. Given the same filter, a universal hash on the smaller packed code is as fast at 1/1000 and faster at 1/100 (Figure 3B), and the same universal hash on every k-mer is 1.12 times as fast as batched ntHash. The two benchmark programs call the same sampler code but give 22.5 and 28.7 cycles per kept k-mer, a difference we have not explained.

![At sketch rates the sampler is faster than batched ntHash, and its gain comes from the filter. Cycles per input k-mer on cs10, $k = 31$, medians of five repeats in one laila job. Sparse rates use inline evaluation, dense rates collected evaluation. (A) The program of Table 2. The circle marks the break-even with batched ntHash, 1/20.2, from the sampler's slope between 1/16 and 1/4. (B) The program of Table 3, with the universal-hash baselines. The panels are not compared with each other, because the two programs give the same sampler code 28.7 and 22.5 cycles per kept k-mer.](figures/fig3_speed.pdf)

## Sampling quality matches ntHash

Across the whole cs10 reference and five seeds per method, minimizer density stayed within 0.6 percent of the expected $2/(w+1)$ at $w = 10$, 20 and 50 for ntHash, two packed-code filters, the fourth-root provable hash and the Eisenstein provable hash. The largest per-bit deviation from one half averaged 0.0009 to 0.0010 for every method, a sanity check that cannot separate hashes below about $10^{-3}$. Every method left 8 to 96 threshold-sampling gaps longer than 200 k-mers per seed, against 0.14 expected by chance. The 1,583 such gaps fell at 151 distinct loci, 62 of them hit by all five methods, and 147 lay in soft-masked repeats. The other 4 loci held one gap each, against about 2 expected by chance in unmasked sequence over 25 runs, so the long gaps reflect repeats rather than any hash.

## Spaced seeds stay strand-consistent for any mask

Applying a 0/1 mask to the coefficients gives a spaced-seed version of the same hash, which rolls with one multiplication per strand per mask. Because $(F_m + \lambda)(R_m + \lambda)$ is strand-symmetric for any mask, the hash of a masked window equals that of its reverse complement. ntHash2's spaced-seed code warns that reverse-complement hashing becomes inconsistent for a non-symmetric mask, and for our non-symmetric mask it disagreed with its own reverse complement at every one of 199,970 windows. Symmetry has a cost. The product depends on the mask and its mirror image together, 29 positions for our 20-position mask, so under 1 percent substitutions it kept 0.747 of windows unchanged ($0.99^{29}$), against 0.818 ($0.99^{20}$) for ntHash2's. The Eisenstein spaced-seed hash ran at 0.75 times the speed of ntHash2's with one mask, over the first 60 million bases of cs10.

## Set sketches recover k-mer differences exactly

We compared the power-sum sketch with km-peeler, the invertible Bloom lookup tables of Shibuya, Belazzougui and Kucherov [@shibuya_2022_genomic-reconciliation; @goodrich_2011_iblt; @eppstein_2011_set-difference], on two pairs of cannabis assemblies. Cannbio-2 and its GigaDB release differ in 3,238,101 canonical 31-mers, 0.73 percent of the set, almost all from unplaced scaffolds that the release drops. Cannbio-2 and its draft differ in 7,931. KMC [@kokot_2017_counting] confirmed both differences element by element. km-peeler ran unmodified with its default three hash functions on the canonical k-mer sets of both assemblies, and each method ran 30 times with independent random seeds on one host.

| method | bytes per difference | exact recovery | rounds |
|---|---|---|---|
| power-sum sketch | 11.3 | in 30 of 30 runs | 2 |
| power-sum sketch sized for one round | 21.7 | yes, in one run | 1 |
| km-peeler, k-mer table at 1.05 times the difference | 20.5 | in 30 of 30 trials | 1 |
| km-peeler, extended syncmers | 4.2 | no, 1,986,807 spurious and 7 missed | 1 |

: Recovering the GigaDB difference, $k = 31$. The sketch started from $2^{16}$ buckets and resent about 1,250 of them in its second round. Sized for one round, it used $2^{17}$. The k-mer table was sized to succeed in 30 of 30 trials.

Between the two assemblies, the power-sum sketch recovered every difference exactly in 30 of 30 runs with 11.3 bytes per difference over two rounds, against 20.5 for a one-round invertible Bloom lookup table sized to succeed in 30 of 30 trials. Sized for a single round, the sketch needed 21.7, so its smaller size comes from the second round. Each 30 of 30 bounds the failure rate only below 9.5 percent (one-sided, 95 percent), so the runs do not compare the two methods' reliability. Neither method returned a wrong difference, but the sketch detects an overfull bucket only with high probability over the random hash. Because an overflowing bucket is split and resent on its own (Figure 4A), the sketch stayed exact when the difference was underestimated eightfold, at 20.8 bytes per difference, while the table had to be resent at double size, at 38.5 (Figure 4B). For a difference of 7,931 between an assembly and its draft, it used 17.3 bytes per difference in one round, against 20.5. The two methods also return different things. Each party of the sketch recovers only its own differing k-mers, while one listing of the table gives both sides'. These results come from one genome and two pairs, each dominated by one structural difference.

![Resending only overfull buckets keeps the sketch small when the difference is underestimated. (A) Each bucket, a range of $S$, holds the count and the power sums $\sum h^j$, $j = 1, \ldots, m + c$. Subtracting two sketches leaves the power sums of the difference. A bucket with more than $m$ differences fails its check and is split and resent alone in the next round. (B) Bytes per difference between Cannbio-2 and its GigaDB release, with buckets or table sized from an estimate of the difference. The table is resent at double size after each failed listing. Labels give the rounds. Both methods ended exact at every point, the sketch in one fixed-seed run per point and the table in all 10 repetitions, whose mean is plotted.](figures/fig4_sketch.pdf)

The sketch hashed k-mers as written. The draft writes 45 percent of its bases on contigs in the other orientation, so a forward-only polynomial hash in the same sketch saw 398,217,248 differences and overflowed every bucket. At $k = 63$, the hash, computed from the sequence with no packed code, recovered all 21,382 differences in 549 KB and one round, matching the canonical difference computed by KMC element by element. The extended-syncmer method of km-peeler, which its authors propose for recovering every difference, also reads each assembly as written. On the GigaDB pair it was 2.7 times smaller than our sketch but returned spurious syncmers numbering 61 percent of the true difference, against about 3 percent in the authors' tests, and missed 7 differences, which we did not trace. Between draft and final, 48,293,748 extended syncmers differed, and the run ran out of memory while listing.

Recovery costs more computation than a table listing. The sketch needs about $m + c = 66$ field multiplications per element, against 3 cell updates for an invertible Bloom lookup table.

# Discussion

For k-mers that fit a machine word, a canonical packed code and a universal hash already give strand symmetry with a proven bound, and in our measurements they ran 1.12 times as fast as batched ntHash and 2.4 times as fast as our hash. The strand-symmetric product earns its place where that recipe is costly or unavailable: for long k-mers, where canonical codes span several words, for spaced seeds and for values that must be added. It gives a collision bound for every pair of inputs, strand symmetry by factorization and a field in which the hashes themselves form a linear set sketch. Because an overfull bucket is resent on its own, the sketch spent 11.3 bytes per difference over two rounds where a one-round table spent 20.5, and it stayed exact at 20.8 when the difference was underestimated eightfold. Sized for one round, it needed 21.7, and it costs about 66 field multiplications per element. Without canonicalization it recovered every difference between an assembly and its draft at $k = 31$ and $k = 63$, where a strand-specific hash saw 398 million at $k = 31$. A table of canonical k-mers is orientation-blind too, so this shows that the hash supplies the symmetry, not that it beats canonicalization.

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
