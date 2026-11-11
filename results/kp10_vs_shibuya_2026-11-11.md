# KP-10 power-sum sketches vs Shibuya et al.'s IBLTs on real near-identical cannabis assemblies (2026-11-11)

**Verdict.** On real near-identical assemblies, both KP-10 variants (Gaussian CM31 and Eisenstein F_q) recover the
k-mer difference EXACTLY. They are about **1.8× smaller than Shibuya et al.'s k-mer IBLT** at equal reliability. They
are **unaffected by contig orientation**, which wrecks km-peeler's proposed extended-syncmer method on a draft-vs-final
pair. km-peeler's extended syncmers are about 3× smaller still, but they return a superset (61% spurious here) and
only work when both assemblies are written in the same orientation. KP-10 loses on compute. Its power sums cost about
m multiplications per element (here m + c = 66), against r = 3 cell updates for an IBLT, and this prototype is numpy.
The two KP-10 fields give the SAME size, because a CM31 element is already 8 bytes.

Code: `recon/` in this repo. Raw JSON for every number: `results/kp10_vs_shibuya_2026-11-11/`.

## A. The real pairs, measured before designing (recon/kmerset.py)

Exact canonical 31-mer sets, whole genome, partitioned into 64 parts by mix64:

| pair | \|A\| | \|B\| | \|A∖B\| | \|B∖A\| | D | D/\|A\| |
|---|---|---|---|---|---|---|
| ANC-1 vs Salk_ANCa | 395,215,003 | 395,215,003 | 0 | 0 | 0 | 0 |
| **Cannbio-2 vs Cannbio-2_draft** | 442,151,894 | 442,159,825 | 0 | 7,931 | **7,931** | 1.8 × 10⁻⁵ |
| **Cannbio-2 vs Cannbio2_GigaDB** | 442,151,894 | 438,915,343 | 3,237,326 | 775 | **3,238,101** | 0.73% |
| Cannbio-2_draft vs Cannbio2_GigaDB | 442,159,825 | 438,915,343 | 3,244,482 | 0 | 3,244,482 | 0.73% |

- ANC-1 and Salk_ANCa have identical k-mer sets (hub_config.yml already excludes Salk_ANCb as a duplicate of ANC-2).
- The cs10 RefSeq vs ASM341772v2 diff ran out of memory and is not used.
- **Primary pair:** Cannbio-2 vs Cannbio2_GigaDB. GigaDB keeps only the 10 chromosomes, so the difference is the
  unplaced scaffolds, almost entirely on A's side.
- **Secondary pair:** Cannbio-2 vs its draft (8,477 contigs), with 7,931 extra draft k-mers, all on B's side.

## B. KP-10 end to end (recon/kp10.py; recon/kp7hash.py, recon/kp7e_hash.py, recon/recon61.py)

**Hashes.** Both are the provable KP-7 S = (F + λ)(R + λ) and both are strand-symmetric. The numpy ports are checked
against the C at every k-mer of a 20,000-k-mer random sample (20,000/20,000 each):
- **Gaussian:** CM31 = F_{p²}, p = 2³¹ − 1, x complex (C hash/kp7m31.c).
- **Eisenstein:** F_q, q = 2⁶¹ − 1, A = 1, C = ω, G = −ω, T = −1 (C hash/kp7e.c). Its 128-bit products are done in
  uint64 by 31-bit limbs (recon/fq61.py, checked against Python integers).

**Sketch.** Per bucket (bucket = top bits of mix64(code), nested inside the 64 partitions): power sums S_1..S_{m+c}
and the size, each counted as 8 bytes. The decode uses Newton's identities and then a Padé step, and each side flags
its own roots. The c = 2 extra sums check the decode, and an overflowing bucket is split and resent (one extra
round). Every run was verified element for element against the exact numpy differences.

**Collision bounds:** Gaussian 2(k−1)/p² = 1.3 × 10⁻¹⁷ per pair of classes; Eisenstein 2(k−1)/q = 2.6 × 10⁻¹⁷ (corrected 2026-11-11 from 1.4 and 2.7, round-1 adjudication L1 F2).
Both entries fit in 8 bytes. The "halving" expected for Eisenstein holds only against a Gaussian sketch over F_{q²}
with q = 2⁶¹ − 1 (16-byte entries), not against CM31.

**Partition 0 calibration, primary pair (D = 50,644):**

| field | buckets (2^bits) | m | B/diff | rounds | splits | exact |
|---|---|---|---|---|---|---|
| both | 2¹⁹ | 16 | 24.6 | 2 | 2 | yes |
| both | 2¹⁸ | 24 | 17.5 | 2 | 4 | yes |
| both | 2¹⁷ | 40 | 14.0 | 2 | 3 | yes |
| both | 2¹⁶ | 64 | **11.3** | 2 | 23 | yes |

**Whole genome** (2¹⁶ buckets, m = 64, c = 2):

| pair | field | sketch | B/diff | rounds | splits | exact | hash time | sketch + decode time |
|---|---|---|---|---|---|---|---|---|
| GigaDB (D = 3,238,101) | Eisenstein | 36.5 MB | **11.27** | 2 | 1,275 of 65,536 | **yes** | 37.0 min | 25.0 min |
| GigaDB | Gaussian | 36.5 MB | **11.27** | 2 | 1,275 | **yes** | 30.2 min | 28.5 min |
| draft (D = 7,931; 2⁸ buckets) | Eisenstein | 137 KB | 17.3 | 1 | 0 | **yes** | 36.3 min | 69.2 min |
| draft | Gaussian | 137 KB | 17.3 | 1 | 0 | **yes** | 28.5 min | 40.8 min |

- **Times** are numpy prototype time for both sides together (one core). The C hashes run about 25 cycles/k-mer, so
  hashing a genome in C takes seconds. Sketch time scales with |A|·(m + c), not with D.
- **The draft pair's size** is set by the 64-partition floor: at least 2⁶ buckets × 67 entries.
- **Mergeability on real data** (recon/kp10.py merge-test): the sketch of Cannbio-2 partition 0 (6,913,034 elements)
  equals the sum of the sketches of two random halves, for both fields.

**A bug found and fixed (771092b).** The Padé matrix read r at negative indices, which wrap to the end of the array in
Python and numpy, whenever a bucket had dB > dA + 1 (more differences on B's side).
- The draft pair is all on B's side, so its first Gaussian run decoded only by splitting nearly every bucket: 3
  rounds, 257 splits, 52 B/diff. It was still exact (`full_kp10_gau_draft.prefix.json`).
- The fix sets r_i = 0 for i < 0, in both this decoder and prototype/recon.py. New one-sided tests fail without it.
- The GigaDB runs (dA ≈ 56, dB ≈ 7 per bucket) never hit it.
- Results that verified exact before the fix stay valid; they only paid extra rounds.

## C. Shibuya et al. (km-peeler, github.com/yhhshb/km-peeler @ b8ec39b; recon/build_kmpeeler.sh, recon/kmpeeler_run.py)

The pinned code builds with gcc 15 only with `-fcommon`; their source is unchanged. r = 3 hash functions and
ε = 0 (their defaults); the IBLT threshold is n = ratio × the true difference; the sketch size is the IBLT file.

**k-mer IBLT (their Sections 3.1 and 3.3)** on the SAME canonical k-mer sets (written as text), on partitions 0–3:

| pair, partitions 0–3 | D | n/D | decoded exactly | B/diff | build (both, C) | diff + list |
|---|---|---|---|---|---|---|
| GigaDB, partition 0 | 50,644 | 0.90 | 0/30 | 17.6 | 1.4 s | |
| GigaDB, partition 0 | 50,644 | 1.00 | 15/30 | 19.6 | 1.4 s | |
| GigaDB | 202,557 | 1.00 | 15/30 | 19.6 | 6.7 s | |
| GigaDB | 202,557 | **1.05** | **30/30** | **20.5** | 8.1 s | 0.3 s |
| draft | 488 | 1.05 | 26/30 | 20.7 | 5.0 s | 0.0 s |

At equal reliability (30/30 at 1.05) the k-mer IBLT needs 20.5 B/diff, against KP-10's 11.27, so 1.8× larger
(66 MB vs 36.5 MB extrapolated to the whole GigaDB pair). Its failures are probabilistic: 50% at n = D, and 13% on the
small draft difference even at 1.05·D. That matches the paper's note that small IBLTs fail more often. KP-10 detects
an overflow deterministically (check sums) and resends only that bucket.

**Extended syncmers (their Section 3.5, the proposed recover-all method)**, k = 31, z = 11 (2k − z = 51), run on the
FASTA AS WRITTEN (km-peeler is strand-specific), with n = 1.05 × the true extended-syncmer difference, 5 trials:
- **GigaDB, whole genome:**
  - **Size:** D_ext = 507,221 differing extended syncmers, IBLT 13.67 MB, which is **4.2 B per true k-mer
    difference**.
  - **Accuracy:** the listing returned 507,219 of 507,221 every trial, and the recovered k-mer set **missed 7** true
    differences and had **1,986,807 spurious** k-mers (61% of D; the paper saw about 3%). The spurious ones are
    k-mers of differing extended syncmers that also occur in shared ones; the unplaced scaffolds share much sequence
    with the chromosomes.
  - **Build:** 25–65 s.
- **Draft vs final, whole genome:**
  - **The difference:** **48,293,748 differing extended syncmers** (of about 55 M per assembly), against a canonical
    k-mer difference of **7,931**.
  - **Canonicalizing the extended-syncmer strings** (recon/canon_lines.py) still leaves **33,446,895**. Which
    positions are syncmers, and the rightward extension, both depend on orientation.
  - **Cost:** the IBLT is 1.3 GB per side, and the run died (out of memory) while listing.
  - **Cause:** the draft's contigs are written in other orientations. KP-7's symmetric hash makes KP-10 independent of
    orientation by construction: the same 7,931 differences, recovered exactly in 137 KB.

## D. Where each wins

| | KP-10 (either field) | km-peeler k-mer IBLT | km-peeler extended syncmers |
|---|---|---|---|
| bytes per k-mer difference (GigaDB) | **11.3** | 20.5 at 30/30 | **4.2**, approximate |
| exact recovery | **yes**, verified | yes, when peeling succeeds | no: a superset (61% spurious) plus 7 missed |
| failure mode | deterministic overflow, then an extra round for that bucket | probabilistic peeling failure (50% at n = D) | as the IBLT, plus the superset |
| contig orientation | **invariant** (symmetric hash) | invariant only if fed canonical k-mers (we did) | **not** invariant: 48 M vs 7.9 k on draft vs final |
| linear / mergeable | yes (sums), verified | yes (sums and XOR) | yes |
| compute per element | about m + c field multiplications (66) | r = 3 cell updates | r per extended syncmer (about 2/21 of k-mers) |
| implementation here | numpy prototype, about 1 h per pair | C, seconds to minutes | C |
| collision bound | Gaussian 1.3 × 10⁻¹⁷, Eisenstein 2.6 × 10⁻¹⁷ per class pair | hash checksums, probabilistic listing | same |

**Defensible advantages of the algebraic sketch:**
1. **Size at equal reliability:** about 1.8× smaller than the exact-recovery IBLT.
2. **Exact recovery with deterministic failure detection,** instead of probabilistic peeling.
3. **Orientation invariance from the hash itself,** which km-peeler's proposed compact method lacks.
4. **A provable per-pair collision bound** on the elements.

Linearity is NOT an advantage: IBLTs subtract too. The extended-syncmer IBLT is smaller when a superset is acceptable
and both assemblies share orientation. A C implementation of the power sums is the obvious next step.
