# KP-9 prototype: KOMPa vs KOMPb chr7, and a near-identical pair (2026-11-10)

**Verdict.** Partitioned characteristic-polynomial reconciliation recovers k-mer set differences exactly. It is no
use for two haplotypes of one plant: their difference is 42% of the set, so the sketch is as large as the set.
It pays off for near-identical sets. A copy of KOMPa chr7 with 1,000 SNPs is reconciled exactly from a 1.13 MB
sketch, against 308 MB for the set (273× smaller). With 100 SNPs the sketch is 141 KB (2,180× smaller).

Code: `prototype/` (field.py, kmers.py, recon.py, run_kompa.py). Tests: `tests/test_recon.py`.
Raw JSON: `bench/data/kmer_poly_work/{sets_summary,sample_*,sim_*}.json` and `runs.log`.

## Setup

- **k-mers:** canonical 31-mers (minimum of the forward and reverse-complement 2-bit codes). N-containing windows
  are skipped, and soft-masking is case-blind.
- **Hash:** h(s) = p_s(x0) over F_{q²} = Z[i]/(q), with q = 2³¹ − 1 (≡ 3 mod 4), A = 1, C = i, G = −i, T = −1, and
  one random x0. The canonical k-mer's polynomial is hashed. KP-7's strand-symmetric pair (e₁, e₂) is a
  follow-up, not done.
- **Buckets:** splitmix64 of the k-mer code, independent of h, with 2^bits buckets. bits is the smallest value
  giving an expected difference of at most m/2 = 32 per bucket. This used the measured D as an oracle; a real
  protocol must estimate D first (for example with a strata estimator or MinHash).
- **Per bucket:** m = 64 solve points, c = 4 check points and one size, so 69 field elements. Each element is 62 bits,
  costed here at 8 bytes.
- **Decode:**
  - Rational interpolation with the degree split fixed by the size difference, by Gauss-Jordan over F_{q²}, with
    free variables set to 0.
  - Each side evaluates the recovered numerator or denominator at its own bucket's hashes. No root-finding.
  - A bucket is accepted only if the flagged sets alone reproduce χ_A/χ_B at the 4 independent check points.
    Otherwise it splits on the next bucket-hash bit.
- **KP-10 power sums:** S₁..S₆₈ per bucket. Newton's identities give the series of χ^rev_{A∖B}/χ^rev_{B∖A},
  Padé recovers the pair, and S₆₅..S₆₈ are the check.

## 1. The haplotype pair: the difference is too large

| chr7 canonical 31-mers | count |
|---|---|
| \|A\| KOMPa | 38,535,935 |
| \|B\| KOMPb | 38,626,335 |
| A∖B | 8,113,246 |
| B∖A | 8,203,646 |
| A∩B | 30,422,689 |
| **D = \|A∖B\| + \|B∖A\|** | **16,316,892 (0.42 \|A\|)** |

Building both exact sets took 82 s.

At 32 expected differences per bucket, this needs 2¹⁹ = 524,288 buckets, i.e. 36,175,872 field elements: **0.94
\|A\|**. That is about the cost of sending the whole set, so reconciliation buys nothing for a haplotype pair.

Correctness on this pair was checked on a random sample of 64 of the 2¹⁹ buckets (4,664 and 4,665 elements, 1,031
and 1,032 differences). Both sketch types gave **exact** recovery, element for element against numpy `setdiff1d`, in
one round with no splits, in 0.4 s (evaluation) and 0.3 s (power sums). The full pair would take about 55 minutes
single-threaded at that rate. It was not run, because it answers nothing the size arithmetic does not.

## 2. Near-identical pairs: where it pays

B is KOMPa chr7 with N random SNPs (non-N positions, a random different base, case kept).

| N SNPs | \|A∖B\| | \|B∖A\| | D | D/\|A\| | buckets | sketch | set | ratio | exact | hash+sort | reconcile |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 1,717 | 2,925 | 4,642 | 1.2×10⁻⁴ | 256 | 141 KB | 308 MB | 2,180× | yes | 130 s | 304 s |
| 1,000 | 16,382 | 28,398 | 44,780 | 1.2×10⁻³ | 2,048 | 1.13 MB | 308 MB | 273× | yes | 126 s | 182 s |

KP-10 power sums on the 100-SNP pair: the same D = 4,642 recovered **exactly**, the same 141 KB sketch, in 157 s
(evaluation sketches: 304 s). The two sketch types decoded the same differences, and power sums are the linear,
mergeable form.

Both runs used one round and no splits. A SNP removes at most 31 canonical k-mers and adds at most 31. The removed
side is smaller (1,717 against 3,100 for 100 SNPs) because many of the removed k-mers occur elsewhere on chr7,
and a k-mer present at another position stays in the set.

The sketch costs 17 bytes per difference at the design load (69 elements × 8 bytes per 32 expected differences),
and 25-30 bytes measured, because bucket counts are rounded up to a power of two. That cost is independent of the
set size. Reconciliation therefore beats sending the set (8 bytes per k-mer) whenever D < roughly 0.3 \|A\|. Its
strength is the regime D ≪ \|A\|.

## 3. Hash collisions

The bound for two distinct k-mers is (k − 1)/q² = 6.5×10⁻¹⁸. Over the 46,739,581-element union of KOMPa and KOMPb
chr7 that gives at most 0.007 expected colliding pairs. Measured: 46,739,581 distinct hashes for 46,739,581
distinct k-mers, so none. The production field q = 2⁶¹ − 1 (IDEAS.md KP-6) would push the bound to about 10⁻³⁵ per
pair, but needs 128-bit products.

## 4. Tests and their calibration

`tests/test_recon.py`, 16 tests:

- random sets with a known difference (0, 5, 7, 24 and 43 elements), both sketch types, exact recovery
- a difference of 60 in one bucket with m = 16: decoding refuses it, splitting recovers it, both sketch types
- direct overflow refusal
- Newton's identities against direct expansion of Π(1 − a t)
- strand invariance of the canonical k-mer set
- field inverses

Mutation check (sed mutants on recon.py, each change asserted to apply):

- **Flipping one decoded flag:** caught. It first ran for minutes, because every bucket failed and split toward 62
  levels. A `max_rounds = 12` cap now turns that into a failure in 15 s.
- **Always mbar = m:** survives, and is equivalent. When m − δ is odd, `(m+δ)//2 + (m−δ)//2` is already m − 1.
- **Skipping the first check-point comparison:** survives, and is redundant. The final check that the flagged
  sets alone reproduce χ_A/χ_B at the same points is strictly stronger, so the first one only exits early.

## 5. What it implies

- **Useful:**
  - assembly versions (a polished vs an unpolished assembly)
  - re-sequenced samples of one line
  - a rerun against its original (KOMPa's 12-of-12 rerun against the 11-of-12 run would be a natural real test)
  - checking that two copies of a large k-mer set agree, or pinpointing where they differ
  - streaming: KP-10 sketches are linear, so chunks and machines merge by addition
- **Not useful:** haplotype pairs or different accessions, where D is a large fraction of the set.
- **Engineering left before any real use:**
  - estimate D before choosing the bucket count
  - q = 2⁶¹ − 1 with 128-bit arithmetic
  - a compiled inner loop: hashing 38.5 M k-mers took 2 minutes in numpy, and evaluating 2 × 68 points over every
    element is most of the reconcile time
  - KP-7's strand-symmetric hashing in place of canonical codes
