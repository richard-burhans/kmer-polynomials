# Round 1 attack, slot L5, lens: inference

Attacker L5 (inference). Attacked text: `reports/rounds/round-1/ATTACKED-paper.md` (draft d049b77 per round.json).
Question asked of every claim: does the sentence the paper writes follow from the evidence it holds?

**Prior rounds.** `scripts/round history "strand-symmetric"` reports no prior mention: round 1 is still open and is the
only round directory. No finding below re-files an adjudicated point.

**Method cards.** They exist at `bench/repos/library/methods/`. I apply `pseudo_replication.md` (F7) and `power.md`
(F7, F9). **Literature:** I read the library summaries for minsky_2003_set-reconciliation, marcais_2024_knonical,
mohamadi_2016_nthash and kazemi_2022_nthash2 (`library/papers_md/summaries/`). Scratch scripts are in
`/home/agent/.claude/jobs/0d0eac3e/tmp/L5/`.

**Steering note.** The claim list frames C14 as "is 0.818 misread". That understates it. The test that produced
0.818 does not measure strand consistency at all (F4). C15 and C16 are framed around the bytes accounting. The larger
gap is that the decoder that produced "exact" uses both parties' full sets (F5).

---

## F1 -- C14: the spaced-seed comparison measures seed WEIGHT, not strand consistency, and ntHash2's number is the correct one

**Verdict: REFUTED** (the inference the paragraph invites).

The sentence chain in the paper runs: ntHash2 "warns that reverse-complement hashing becomes inconsistent for a
non-palindromic mask". Then "our hash kept 0.747 of windows unchanged, the exact double-stranded value, against 0.818
for ntHash2's". A reader takes 0.818 as a symptom of the inconsistency. It is not.

- `ideas/eisenstein/bench_spaced_e.cpp:37`: mask 1 is `1101100111010111011100101110101`. It has weight 20. Its
  union with its mirror covers 29 positions. Mask 0 is palindromic, and on mask 0 all three methods agree exactly
  (`results/paper_bench_2026-11-11/run.log:351,355`).
- `bench_spaced_e.cpp:149-157`: the "exact" fraction is the fraction of windows whose care bases of s AND of rc(s)
  are unchanged. That means unchanged at all 29 positions of mask ∪ mirror(mask).
- Arithmetic (independent substitutions at rate r): (1−r)^29 = 0.9714 / 0.74717 / 0.22594 at r = 0.001 / 0.01 /
  0.05. (1−r)^20 = 0.98019 / 0.81791 / 0.35849. Measured (`run.log:354,358` and the 0.05 line): ours 0.97146 /
  0.74722 / 0.22567. SeedNtHash 0.98023 / 0.81794 / 0.35839. **ntHash2 behaves exactly as a weight-20 spaced seed
  should. Ours behaves as a weight-29 seed.**
- No step of the test compares a window's hash with the hash of its reverse complement. It cannot detect ntHash2's
  warned inconsistency, and it does not.

**What follows.** A strand-symmetric hash whose value depends only on the care positions is impossible for a
non-palindromic mask. If h(s) depends only on s|m and h(s) = h(rc s), then h depends only on the positions in
m ∩ mirror(m). (F_m+λ)(R_m+λ) gets its symmetry by keying on m ∪ mirror(m). It silently turns the weight-20 seed into a
weight-29 one, and it loses seed sensitivity: 0.747 vs 0.818 at 1%, and 0.226 vs 0.358 at 5%. The paper presents
this cost as exactness.

**Correction.** Replace the 0.747/0.818 sentence with: "For a non-palindromic mask, strand symmetry makes our hash
depend on the mask's palindromic closure (29 of 31 positions for a weight-20 mask). At 1% substitutions it kept 0.747
of windows unchanged, against 0.818 for ntHash2, which hashes the 20 care positions of one strand. This test does not
measure ntHash2's strand inconsistency." The headline "strand-consistent for any mask" can stay only with that cost
stated.

## F2 -- C11: by the paper's own model, the sampler is SLOWER than batched ntHash at 1/20

**Verdict: REFUTED** (as worded in the abstract).

- `results/paper_bench_2026-11-11.md`, "Per-kept cost and break-even": slope 28.5, intercept 9.96, "beats batched
  ntHash (11.37) below a sampling rate of **1/20.3**". At 1/20 the same line gives 9.96 + 28.5/20 = 11.385 cycles,
  which exceeds 11.37. The abstract says "faster than batched ntHash at sampling rates of 1/20 and below". That is
  false at its own endpoint, by its own fit. The Results section hedges ("about 1/20") and the abstract drops the
  hedge.
- The break-even rests on two collected-mode rows (1/16 and 1/4). The 0.015-cycle margin at 1/20 is 0.13%. The
  spread between the two repeats is up to 0.8% on those rows, and the same benchmark reports 3–6% between jobs. In
  that same file, break-even moved 1/21.5 → 1/19 → 1/20.3 across drafts. A 3% shift in the 11.37 baseline moves
  break-even from about 1/16 to about 1/27.
- The only measured point near the claimed boundary is 1/16. There the batched sampler is SLOWER (11.74, 0.97×), and
  the branchy one is slower still (13.32, 0.85×). Speed is measured "below" break-even only at 1/100 and 1/1000.
- Premise check of the claim list: the brief says 1/20 is interpolated "using DIFFERENT evaluation modes". That is
  wrong. The 28.5 slope uses the batched mode only. The paper's Table, however, mixes modes (inline rows at 1/1000 and
  1/100), so the "28.5 per kept k-mer" figure does not describe the inline rows it sits beside: the inline slope is
  (10.01 − 9.40)/0.009 = 68 cycles per kept k-mer.

**Correction.** Abstract: "...is faster than batched ntHash at sampling rates of 1/100 and below (1.14× at 1/100), and
breaks even near 1/20". Results: give break-even as "between 1/16 and about 1/25, within run-to-run variation".

## F3 -- C2: "F·R cannot tell a k-mer from its plain reverse" does not explain the counts it is cited for

**Verdict: NARROWED.**

- F·R = −x^{k−1}p(x)p(1/x) is invariant under p → −p as well, so it also merges every k-mer with its plain
  COMPLEMENT. The orbits of {id, reverse, complement, reverse-complement} number 4,160 at k = 7 and 16,512 at k = 8
  (Burnside, and enumerated in `tmp/L5/fr.py`).
- Exact (characteristic-0) enumeration of p(x)p(1/x) over Gaussian integers gives **4,160 at k = 7 and 16,362 at
  k = 8** (`tmp/L5/fr.py`). The table's k = 8 value is therefore 150 values below what reversal plus complement explain.
  These are homometric sequences, distinct polynomials with the same autocorrelation. They are a structural collapse,
  not finite-field chance: about 16,512²/(2·4.6×10¹⁸) ≈ 3×10⁻¹¹ chance collisions would be expected.
- `prototype/kp7_symmetric_check.py` tests the stated mechanism on ONE k-mer (`e2(s)==e2(reverse(s))`). It never tests
  complement merging or homometry.
- Caption: the check runs over F_{q²} = Z[i]/(q) with q = 2³¹−1 (script line 1). The caption's "q = 2^{31} − 1" reads
  as the prime field, which cannot hold i (q ≡ 3 mod 4), as the paper itself says for 2⁶¹−1.

**Correction.** Methods: "The product F·R is an autocorrelation. It cannot tell a k-mer from its reverse or from its
complement, and at k = 8 it also merges 150 further homometric classes." Caption: "over F_{q²}, q = 2³¹ − 1".

## F4 -- C1: the bound is right, but its arithmetic is off and three inferences around it overreach

**Verdict: NARROWED.** The core step survives. In F_q[x][λ], S_s = (λ+F_s)(λ+R_s) with monic linear factors, so
S_s = S_u as polynomials forces {F_s, R_s} = {F_u, R_u} and hence {s, rc s} = {u, rc u}. The difference polynomial
has total degree max(k, 2k−2) = 2(k−1) for k ≥ 2. The defects are in what is drawn from it:

1. **Arithmetic.** 2(k−1)/q at k = 31 is 60/(2⁶¹−1) = **2.60×10⁻¹⁷**, not 2.7. The 2.7 is 2k/q = 2.69×10⁻¹⁷. The
   Methods sentence states 2(k−1)/q and then quotes the other quantity. Likewise the F_{p²} bound ("about 1.4×10⁻¹⁷",
   Implementation) is 60/(2³¹−1)² = 1.30×10⁻¹⁷, or 1.34×10⁻¹⁷ for 2k. It is 1.3 either way.
2. **"Exhaustive enumeration agrees."** One random (x, λ) gives 8,192 and 32,896 distinct values. A uniformly random
   61-bit function would collide among 32,896 inputs with probability about 32,896²/(2q) ≈ 2×10⁻¹⁰. The check
   therefore cannot "agree with" a bound or with injectivity. It verifies the code, not the theorem. A check that
   could fail would test the polynomial identity S_s ≠ S_u symbolically, or show F+R and F·R failing while S does not
   at the same tiny field size.
3. **"Proven, not measured" (abstract line 1) against what the proof delivers at genome scale.** Cannbio-2 has
   442,151,894 canonical 31-mers (`results/kp10_vs_shibuya_2026-11-11.md` §A). Summed over pairs, the bound gives
   9.8×10¹⁶ × 2.6×10⁻¹⁷ ≈ **2.5 expected colliding pairs in one genome**. The proven guarantee therefore does not
   exclude collisions within a single assembly. It is 60 times weaker per pair (2(k−1)/q against 1/q) than the
   ideal-hash rate the measured hashes are judged against. The abstract's contrast implies a proof buys more than it
   does here.

**Correction.** "2.6 × 10⁻¹⁷" (or state the bound as 2k/q throughout), and "1.3 × 10⁻¹⁷" in Implementation. Change
"Exhaustive enumeration agrees" to "Exhaustive enumeration found no collision, as expected for any good hash at this
size". Add one sentence: "Summed over the 4.4 × 10⁸ canonical 31-mers of an assembly, the bound allows about 2.5
colliding pairs."

## F5 -- C15: "11.3 against 20.5 of equal reliability" compares an oracle decoder with a one-message IBLT, and the paper's own numbers deny equal reliability

**Verdict: REFUTED.**

1. **The decoder is not a party.** `recon/recon61.py:100-136`, `decode(SA, SB, nA, nB, m, c, ownA, ownB)`, flags roots
   among BOTH sides' element arrays. The overflow test checks the c extra sums against `ownA[flagsA]` and
   `ownB[flagsB]` together, so it needs both parties' full sets in one process. `recon/kp10.py:80-81` then scores
   "exact" on the k-mer codes found that way. A receiver holding one sketch and its own set cannot run this check. It
   also cannot turn the other side's recovered hash values into k-mers: S is not invertible, and Minsky's elements ARE
   the field elements, whereas here they are hashes (minsky summary, "Limits"). km-peeler's k-mer IBLT lists the
   actual k-mers of both sides from one table. The two "bytes per difference" therefore buy different outputs, and
   only the IBLT's was produced under a realistic protocol.
2. **The reliabilities are not equal on the paper's own data.** At the 1.05× setting the k-mer IBLT succeeded 30/30
   on four GigaDB partitions. The same setting failed 4/30 on the draft pair (`kp10_vs_shibuya` §C table: 26/30),
   which the paper reports two sentences later. By the rule of three (`library/methods/power.md`: a zero is
   uninformative unless P(zero) is small), 30/30 bounds failure only below about 10%. On the other side, KP-10's
   "every run" is two whole-genome runs on this pair (one per field) plus partition-0 calibrations.
3. **Rounds differ.** KP-10 on GigaDB took 2 rounds with 1,275 bucket splits (`full_kp10_eis_giga.json`:
   `rounds_max` 2, `splits` 1275; the bytes do include round 2, since 4,561,762 − 65,536·67 = 1,275·2·67). The IBLT
   is one message. The paper says "one round" only for the draft pair, which invites the reader to assume the same
   for GigaDB.

**Correction.** Abstract: "...it needed 11.3 bytes per difference over two rounds, against 20.5 for an invertible
Bloom lookup table sized to succeed in 30 of 30 trials on a subset. Each side recovers its own differing k-mers and
the other side's hash values." Drop "of equal reliability" in the abstract, the Results ("for the same reliability")
and the Discussion ("0.55 times the bytes ... of equal reliability"). Also state that the measured decoder had
access to both sets.

## F6 -- C16: overflow detection is probabilistic, and the code says so

**Verdict: REFUTED.**

- An overflowing bucket passes the check if the signed multiset (true difference minus flagged elements) has
  vanishing power sums 1..m+c. A single hash collision between an A-only and a B-only k-mer (h_a = h_b) makes them
  vanish identically, for every j. Detection is therefore conditioned on the hash, and its error is governed by the
  same per-pair bound as C1. It is not deterministic.
- `recon/recon61.py:157`, the authors' own decoder: "identical sketches: no difference (with probability 1 − O(m/q))".
- Minsky 2003, the cited source: rational-function equality "is tested probabilistically, with error per check at most
  (D − 1)/|E|" (library summary, Method, "Unknown bound").
- The Discussion's "detected its own overflow deterministically" inherits this. In the measured runs, detection also
  used root membership in both parties' sets (F5.1), a stronger check than the c sums the Methods describe.

**Correction.** Methods: "...so a bucket holding more than m differences is detected, except with probability
bounded by the hash's collision bound, rather than resolved wrongly." Discussion: delete "deterministically".

## F7 -- C13: the bit-balance "match" is saturated by repeated k-mers and cannot discriminate hashes; "none outside" is a powerless zero

**Verdict: NARROWED.**

- `hash/stats_seeds.cpp:61` counts bits over every k-mer occurrence (736,547,799). For an ideal hash the SD of one
  bit's mean is 0.5/√n = 1.8×10⁻⁵, and the expected maximum over 64 bits is about 4.5×10⁻⁵. The reported 0.00088
  (ntHash) and 0.00090 (Eisenstein) are about 20 times that. Eisenstein's S is a uniform field element whose bits have
  bias near 2⁻⁶¹, so the deviation is the genome's repeated k-mers. Those occurrences are not independent draws
  (`library/methods/pseudo_replication.md`: n is the number of independent decisions, not rows). The statistic
  measures cs10, not the hash, and "0.00090 against 0.00088" is no evidence that the two hashes match in quality. Only
  a defect larger than about 0.001 would show.
- "None of the Eisenstein hash's 367 lay outside them" (`kp7_stats_seeds_eisenstein_2026-11-11.md`). The pooled rate
  outside repeats is 2/1,774, so 0.41 such gaps are expected among 367, and P(zero) = e^−0.41 = 0.66
  (`library/methods/power.md`). The zero is the modal outcome and says nothing about the Eisenstein hash.
- Minor: 1,774 = 341+363+371+332+367 pools FIVE methods (two packed-code gates). The paragraph lists four.
- "So the long gaps reflect repeats" itself follows: a gap counts as inside when at least half its bases are
  soft-masked (`stats_seeds.cpp:16-17,84-86`), and 1,772/1,774 is far above the 47.7% masked share.

**Correction.** Report bit balance on DISTINCT k-mers, or drop it as evidence of hash quality. Delete "none of the
Eisenstein hash's 367 lay outside them" or give P(zero) beside it. Say "five methods".

## F8 -- C17: the orientation cause of the 48 million is asserted, never measured

**Verdict: NARROWED.**

- No artifact measures the draft's contig orientations. `grep -i orientation` over `results/*.md` and the ledger finds
  only the assertion in `kp10_vs_shibuya_2026-11-11.md` §C ("**Cause:** the draft's contigs are written in other
  orientations"). Canonicalizing the extended-syncmer strings still leaves 33,446,895 differences, and the residual is
  attributed to orientation-dependent selection without a test. The direct test was not run: orient the draft's
  contigs to the final assembly and rerun km-peeler.
- The abstract sets "where a strand-specific method saw 48 million" beside the IBLT comparison. The Discussion concedes
  that the orientation advantage "holds against the compact extended-syncmer method only". The abstract and
  Introduction do not say so, and the orientation-blind k-mer IBLT that the paper did run on this pair (26/30 at 1.05×)
  is absent from both.
- 137 KB is a floor set by the 64-partition layout (2⁸ buckets × 67 entries × 8 bytes = 137,216; `kp10_vs_shibuya` §B: "size is set by the 64-partition
  floor"), at 17.3 B/difference. It is not a property of D = 7,931.

**Correction.** Abstract: "...where km-peeler's extended-syncmer method, which reads contigs as written, saw 48
million. A canonical k-mer table is orientation-blind too." Methods/Results: "we attribute this to contig orientation;
we did not test it directly."

## F9 -- C9: the upper end of the range comes from a repeat-filtered subset, not from another assembly pair

**Verdict: NARROWED.**

`results/ideas_eisenstein_2026-11-10.md:21-26` lists three pairs on NC_044378.1, with Ts/Tv 1.585, 1.584 and 1.669
(61.3–62.5% transitions). It then lists "same, uppercase only (outside soft-masked repeats)" for KOMPa vs KOMPb at
1.847 (64.9%). The paper's "61 to 65 percent (Ts/Tv 1.58 to 1.85)" in pairwise alignments of cannabis assemblies
mixes three pairs with one filtered subset of one of them. The arithmetic (1.58/2.58 = 61.2%, 1.85/2.85 = 64.9%) is
right. Two of the three pairs share cs10, and the third is two haplotypes of one plant.

**Correction.** "61 to 63 percent (Ts/Tv 1.58 to 1.67), and 65 percent outside soft-masked repeats".

## F10 -- C18: "fixed seeds, so measured, not proven" is a false dichotomy for a GF(2)-linear hash

**Verdict: NARROWED.**

The paper itself says ntHash "is linear over GF(2)". Per the mohamadi_2016_nthash summary, the per-strand hash is the
XOR of rotated seeds. With fixed seeds it is a fixed linear map from 2k = 62 bits to 64, so its collisions on all
31-mers can be decided exactly by one rank computation. That is a proof, deterministic and stronger in kind than a
per-pair probabilistic bound, not a measurement. What fixed seeds actually cost is adversarial predictability and the
same collisions in every run, which the Introduction does say. Separately, "linear over GF(2)" does not hold for the
canonical value. ntHash1 takes min(f, r) and ntHash2 takes f + r mod 2⁶⁴ (kazemi_2022_nthash2 summary, lines 95-97).
The Discussion's "constant-time substitution and concatenation hold for ntHash too, which is linear over GF(2)" is
therefore true only per strand.

**Correction.** Abstract: "...do this with fixed seeds, so an input that collides does so in every run; ours carries a
collision bound over random seeds." Discussion: "...hold for ntHash's per-strand hashes, which are linear over GF(2)".

## F11 -- C19: "no visible effect on fractional sampling" is contradicted by the table the summary rebuilt

**Verdict: NARROWED.**

`library/papers_md/summaries/marcais_2024_knonical.md:25-29`: in Table 2, rebuilt from the PDF, fractional sampling's
deserts shrink on chromosome 1 at k = 31 (63.20 → 38.14 Mb) and grow 2.3× on chromosome 10 (24.68 → 57.52 Mb). The
paper notes high variance (SD 34.7 Mb). "No visible effect" overreads the caption. The paper's "show that folding ...
negates the window guarantee" also upgrades the source's "can negate" / "may not preserve" (summary, line 12).
Fractional sampling has no window guarantee to negate, so the sentence's conclusion, that long gaps are no weakness
for a sketch, does not follow from the cited result.

**Correction.** "...can negate the window guarantee of syncmer-style schemes. Fractional sampling has no such
guarantee, and its desert sizes vary widely in either direction after canonicalization."

---

## Claims attacked and held

- **C3 SURVIVED.** Re-derived: F_{j+1} = (F_j − s_j)x⁻¹ + s_{j+k}x^{k−1}. With R_j = −Σ s_{j+t}x^{k−1−t}, I get
  xR_j + s_j x^k = −Σ_{t'=0}^{k−2} s_{j+1+t'}x^{k−1−t'}, so R_{j+1} = x(R_j + s_j x^{k−1}) − s_{j+k}. Correct. The
  single-strand bound needs x ≠ 0 (x⁻¹ is used), so it is (k−1)/(q−1), which is negligible.
- **C4 SURVIVED.** 2⁶¹ ≡ 2 (mod 3), so q ≡ 1 (mod 3); q ≡ 3 (mod 4). Cycles 47.37/24.69 = 1.92 and instructions
  179.3/78.5 = 2.28 (`paper_bench_2026-11-11.md`). "About twice as many multiplications" is inferred from cycles, not
  counted, but ~2 holds.
- **C5 SURVIVED** (identity and "exactly when"). |A+C|² − |A−C|² = 2Re(AC̄) − (−2Re(AC̄)) = 4Re(AC̄). A fourth-root map
  with complement as negation needs C = ±iA, so Re(AC̄) = 0. Inference note: "≠ 0" is the criterion for SEPARATION
  only. Re(AC̄) > 0 (e.g. C = e^{iπ/3}) makes transitions the costliest, so the motivating property, transitions
  cheapest, needs Re(AC̄) < 0. I did not run the priority search (prior-art lens).
- **C8 SURVIVED.** |1+ω|² = 2 + 2Re ω = 1, |1−ω|² = 3, |2|² = 4. The sentence that follows says the hash keeps the
  algebra, not the distances, so it does not imply the hash measures distance.
- **C10 SURVIVED.** Every table value matches `paper_bench_2026-11-11.md`, and the spreads are ≤1.3% except the
  scalar reference roller (4.0%, not in the paper). Note that the library interface (23.53) is faster than the
  Eisenstein hash (24.69).
- **C12 SURVIVED** on "about twice": 23.53/11.37 = 2.07. I did not inspect the bit-for-bit check (reproduction lens).

## Verdicts

| C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 | C9 | C10 | C11 | C12 | C13 | C14 | C15 | C16 | C17 | C18 | C19 | C20 | C21 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| NARROWED | NARROWED | SURVIVED | SURVIVED | SURVIVED | UNCHALLENGED | UNCHALLENGED | SURVIVED | NARROWED | SURVIVED | REFUTED | SURVIVED | NARROWED | REFUTED | REFUTED | REFUTED | NARROWED | NARROWED | NARROWED | UNCHALLENGED | UNCHALLENGED |
