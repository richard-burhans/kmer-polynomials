# Round 1 adjudication of L5 (inference)

Adjudicator: independent skeptic. Attack read first, then the claims (fixed order).
**Seal check:** `sha256(attack-L5-inference.md)` begins `b1430fe37d83`, which matches the sealed digest. Adjudication proceeds.
**History:** `scripts/round history "strand-symmetric" --kind adjudication` returns nothing. Round 1 is the only round
and is still open, so no finding re-files a settled point.
**Method cards:** `refs/methods/` does not exist in this repo. The cards live at `library/methods/`. `power.md` governs
F7's zero, and `pseudo_replication.md` governs F7's bit balance. No card governs F1–F6 or F8–F11. Those are algebra,
code-reading or source-reading disputes, and I decided them by re-running the work.
**Re-runs:** scratch is `/home/agent/.claude/jobs/0d0eac3e/tmp/ADJ-L5/` (`fr.py`), plus inline checks quoted below.

---

## F1 (C14): spaced-seed test measures seed weight, not strand consistency. UPHELD

**Premise re-run.** Mask 1 `1101100111010111011100101110101` has weight 20. Its union with its mirror covers 29
positions, and its intersection covers 11. Mask 0 is palindromic with weight 24, and (0.99)^24 = 0.78568, against
0.78573 measured. Expected unchanged fractions at r = 0.001/0.01/0.05: (1−r)^29 = 0.97140/0.74717/0.22594, and
(1−r)^20 = 0.98019/0.81791/0.35849. `run.log` gives ours 0.97146/0.74722/0.22567 and SeedNtHash 0.98023/0.81794/0.35839.
Both columns fit their binomial curves to the 4th decimal. `bench_spaced_e.cpp:149-160` never compares h(s) with
h(rc s). It only compares h(s) with h(mutated s), and its "exact" column counts unchanged care bases of both strands.
The test therefore cannot see ntHash2's warned inconsistency.

The impossibility step is also sound. If h depends only on s|m and h(s) = h(rc s), then h depends only on
m ∩ mirror(m). So the 0.747 vs 0.818 gap is the price of symmetry: a weight-20 seed is effectively turned into a
weight-29 one. Placing it right after ntHash2's warning invites the reader to read the gap as inconsistency.
The sentence "0.747 ... the exact double-stranded value" is true but misleads.

**Repair.** I endorse it, with one tightening. "Which hashes the 20 care positions of one strand" is inferred from the
fit, not read from ntHash2's code. Write "whose unchanged fraction matches a 20-position seed". With that change the fix
is free (wording only, no re-measurement).

## F2 (C11): break-even at 1/20. NARROWED (diagnosis upheld, "REFUTED" overstated)

**Re-run.** Slope (17.08 − 11.74)/0.1875 = 28.48. Intercept 9.96. At 1/20 the line gives 11.384 against 11.37, so
break-even is 1/20.2. `paper_bench_2026-11-11.md:67` confirms the 3–6% between-job variance. A ±3% move in 11.37 shifts
break-even to about 1/16.3 or 1/26.7, as the attacker says. The abstract's "1/20 and below" therefore claims a boundary
the data cannot place. The defect is that the claim is unsupported, not that it is shown false. The model's own error
at 1/20 is 0.014 cycles (0.13%), far inside noise, so REFUTED is the wrong verdict. Measured speed-ups exist at 1/100 and
1/1000 in both modes: inline 1.14×/1.21×, and batched 1.11×/1.14× (the batched rows are not in the paper's table). At
1/16 both modes are slower.

The attacker's premise correction of the claim list is right: the 28.5 slope uses batched rows only. The point about
mixed modes is also right. The table's sparse rows are inline (about 63–68 cycles per kept k-mer), while "28.5 per kept
k-mer" describes the batched mode.

**Repair.** I endorse the abstract change "faster than batched ntHash at 1/100 and below (1.14×), breaking even near
1/20" (free). For Results I endorse "break-even near 1/20, between the measured 1/16 (slower) and 1/100 (faster)" in
place of the attacker's "between 1/16 and about 1/25", which states a bound nobody measured (free). Results should also
say that 28.5 is the batched sampler's cost (free).

## F3 (C2): F·R count at k = 8 is not explained by reversal. NARROWED, as filed

**Re-run** (`ADJ-L5/fr.py`). Exact enumeration of p(x)p(1/x) over Z[i] gives 4,160 values at k = 7 and 16,362 at k = 8.
Orbits of {id, rev, comp, rc} number 4,160 and 16,512. The table's 16,362 is confirmed, and 150 classes merge beyond the
reversal symmetry. Those merges are homometric pairs, structural rather than chance. `kp7_symmetric_check.py` line 1
confirms the field is F_{q²} = Z[i]/(q), q = 2³¹−1, so the caption's bare "q = 2^{31}−1" misdescribes it.

**Correction on the attacker.** The "also merges with its plain COMPLEMENT" leg is not an independent mechanism. Within
double-stranded classes {s, rc s}, merging with rev s already implies merging with comp s = rc(rev s), so the paper's
"plain reverse" accounts for the whole orbit. The paper never claimed the count equals the orbit count. It says the
enumeration "confirms" a collapse, and it does.

**Repair.** Endorse "at k = 8 it also merges 150 further homometric classes" (free, verified). Endorse the caption fix
"over F_{q²}, q = 2³¹−1" (free). "Or from its complement" is true but redundant, so leave it optional.

## F4 (C1): bound arithmetic and surrounding inferences. NARROWED

1. **Arithmetic: UPHELD.** 60/(2⁶¹−1) = 2.602×10⁻¹⁷ and 62/(2⁶¹−1) = 2.689×10⁻¹⁷. The Methods sentence names 2(k−1)/q
   and quotes 2.7, so it is internally inconsistent. The abstract's "about 2k/q, 2.7×10⁻¹⁷" is consistent. For the
   F_{p²} bound, 60/(2³¹−1)² = 1.301×10⁻¹⁷ and 62/(2³¹−1)² = 1.344×10⁻¹⁷, so "about 1.4×10⁻¹⁷" is
   wrong under either convention. The same 1.4 appears in `results/kp10_vs_shibuya_2026-11-11.md:43` and its table D,
   and those must change together. I checked the whole set: the abstract (2k/q, 2.7) and the Introduction ((k−1)/q,
   single strand) are unaffected. **Fix: "2.6 × 10⁻¹⁷" in Methods and "1.3 × 10⁻¹⁷" in Implementation and the kp10 note.
   The fix is free.**
2. **"Exhaustive enumeration agrees": NARROWED, and the repair is OVERTURNED.** The attacker is right that finding no
   collision among 32,896 values cannot confirm a 2.6×10⁻¹⁷ bound. The attacker is wrong that the check "verifies the
   code, not the theorem" and could not fail. In the same table at the same field, F+R and F·R do fail structurally,
   so the enumeration has a working negative control for structural non-injectivity, which is exactly what unique
   factorization rules out. The attacker's own suggested test ("show F+R and F·R failing while S does not") is
   already in the paper. The proposed wording "as expected for any good hash at this size" throws away that control.
   Endorse instead: "Exhaustive enumeration finds no collision, where F+R and F·R collapse" (free).
3. **Genome-scale reading: NARROWED.** The arithmetic is right: 442,151,894 canonical 31-mers give about 9.8×10¹⁶
   pairs, times 2.6×10⁻¹⁷, which is 2.54 expected colliding pairs as an upper bound, against 0.042 for an ideal 61-bit
   hash. No sentence in the paper claims zero collisions per genome, so this is missing context, not a false statement.
   Adding "summed over the 4.4×10⁸ canonical 31-mers of an assembly, the bound allows at most about 2.5 colliding
   pairs" is accurate and free. I recommend it, but it is not owed as a correction.

## F5 (C15): "11.3 against 20.5 of equal reliability". NARROWED (numbers stand; framing falls)

**Re-run.** `full_kp10_eis_giga.json` gives sketch_elements 4,561,762 = 65,536·67 + 1,275·2·67, 36,494,096 bytes,
11.27 B/diff, rounds_max 2. `full_kmp_kmers_giga_r1.05.json` gives 30/30 at 20.53 B/diff (4 partitions,
D = 202,557). The draft file gives 26/30. Both byte counts are one message's worth, so the numbers reproduce and the
accounting is like for like.

- **"Equal reliability": UPHELD.** One side has 30/30 on a subset, which per `power.md` bounds failure only below about
  10% (rule of three). The same setting failed 4/30 on the draft pair. The other side has a detect-and-resend scheme
  whose failure mode is a different kind of event. The two reliabilities were never measured on one scale.
- **Oracle decoder: UPHELD as a fact about the measured code.** `recon61.py:100-136` flags roots by membership in
  ownA and ownB in one process and checks the c sums on both flagged sets. S is not invertible, so a party without the
  other's set would get the other side's differences only as hash values (after factoring P, which is not
  implemented). The IBLT lists actual k-mers. Rule 8: **this point was raised by L5 alone** (no other attack mentions
  the two-set decoder), and it is correct.
- **Rounds: UPHELD, minor.** GigaDB took two rounds, and the bytes include round 2.
- **Why not REFUTED:** 11.3 and 20.5 both reproduce, and the sketch's byte count does not depend on who decodes.

**Repair.** Endorse dropping "of equal reliability" and "for the same reliability" in the abstract, Results and
Discussion (free). Endorse "over two rounds" (free) and "against 20.5 for an IBLT that succeeded in 30 of 30 trials on
4 of 64 partitions" (free). **OVERTURNED:** the attacker's sentence "Each side recovers its own differing k-mers and the
other side's hash values". It describes a two-party decoder that was never built or run, so it would install a new
unmeasured claim. Write instead: "Our decoder ran with both sets in one process. A two-party decoder, which would return
the other side's hashes rather than its k-mers, was not built." (free)

## F6 (C16): overflow detection is not deterministic. UPHELD on the defect; the repair is NARROWED

The defect holds. Whether an overflowing bucket passes the c check depends on the hash values. An A/B hash collision
cancels in every power sum and yields a wrong-but-consistent decode. Beyond that case, the c extra equations are
satisfied by an overflowing residue with small but nonzero probability over (x, λ). The authors' own
`recon61.py:157` says "with probability 1 − O(m/q)", and Minsky's check is probabilistic. So "deterministically" in the
Discussion (and in `kp10_vs_shibuya` §C and table D) is false.

The repair is wrong as written. "Bounded by the hash's collision bound" names a specific bound that has not been
derived. The check involves power sums up to degree m+c of S values, and the flagged sets themselves depend on (x, λ),
so the per-pair bound 2(k−1)/q is not the governing quantity. Endorse instead: "...is detected, except with negligible
probability over the random choice of x and λ, rather than resolved wrongly", and delete "deterministically" in the
Discussion (free). One correction on the attacker: the membership-based flagging in the measured runs makes detection
stronger than the c sums alone, not weaker, so F6's last bullet supports the claim's result while undercutting its
description.

## F7 (C13): bit balance and the zero. NARROWED, as filed

`stats_seeds.cpp:61` counts bits over all 736,547,799 occurrences. The iid SD is 1.8×10⁻⁵, and the expected maximum of
64 is about 4.5×10⁻⁵. All five hashes, very different in construction, land at 0.00088–0.00103, about 20× that. Under
`pseudo_replication.md` repeated k-mers are not independent draws. "0.00090 against 0.00088" therefore does not
discriminate hash quality below about 10⁻³.

For "none of 367 outside", 367·2/1,774 = 0.41 expected and P(zero) = 0.66. Under `power.md`, a modal zero is no
finding. 1,774 = 341+363+371+332+367 pools five methods, while the paper names four (the two packed-code gates are
folded into one). The attacker is also right that "so the long gaps reflect repeats" follows: 1,772/1,774 against
47.7% masked.

**Repair.** Deleting "none of the Eisenstein hash's 367 lay outside them" or giving P(zero) beside it is free. Saying
"five methods" is free. Rewording bit balance as "a sanity check, insensitive below about 10⁻³" is free. Recomputing it
on distinct k-mers costs a re-run.

## F8 (C17): orientation cause untested. NARROWED, as filed

There is no artifact that measures draft contig orientation: grep finds only the assertion in `kp10_vs_shibuya` §C, and
nothing in `bench/data/kmer_poly_work/`. Indirect evidence (7,931 canonical k-mer differences against 48 M strand
strings, cut to 33 M by canonicalizing) makes orientation plausible, but other string-level causes were not excluded.
137 KB is the 2⁸ × 67 × 8 floor (17.3 B/diff), and the source file says so.

**Repair.** Endorse "we attribute this to contig orientation; we did not test it directly" (free). Endorse naming the
method in the abstract ("km-peeler's extended-syncmer method") together with the Discussion's concession that a
canonical k-mer table is orientation-blind too (free).

## F9 (C9): Ts/Tv range mixes a filtered subset. NARROWED, as filed

`ideas_eisenstein_2026-11-10.md:21-26` gives three pairs at 1.584–1.669 (61.3–62.5%). 1.847 (64.9%) is KOMPa vs KOMPb
outside soft-masked repeats. The arithmetic of the range is right.

**Repair.** Endorse "61 to 63 percent (Ts/Tv 1.58 to 1.67), and 65 percent between the two haplotypes outside
soft-masked repeats" (free; "between the two haplotypes" is my addition, because the filtered row is one pair only).
The source also says the KOMPa/KOMPb comparison runs through cs10 positions, so "pairwise alignments" should note it is
projected onto cs10.

## F10 (C18): "measured, not proven" a false dichotomy. NARROWED (main leg OVERTURNED; GF(2) wording leg UPHELD)

**Premise test, main leg.** The attacker argues that ntHash with fixed seeds is "a fixed linear map from 2k = 62 bits to
64", so collisions are "decided exactly by one rank computation". That requires the base→seed map to be affine over
GF(2) in the 2-bit encoding, which needs seedA⊕seedC⊕seedG⊕seedT = 0. The library's seeds (`ntHash/src/internal.hpp:125-128`)
XOR to 0x047f4d133ed3484a, not zero. The per-strand hash is XOR-additive over positions but not linear in the 2-bit
input. Each position contributes one of 6 nonzero differences from a 3-dimensional space that misses one element, so
collision-freeness is a constrained kernel question, not a rank computation. No one has published such a proof, and
"their quality is established by measurement" accurately describes how ntHash's authors established it.

**What the attacker should have checked:** whether the seed table is affine, before calling the hash a linear map of
the 2-bit code. The paper's own phrase "linear over GF(2)" misled the attacker too (see the next leg). The suggested
abstract rewording ("an input that collides does so in every run; ours carries a bound over random seeds") is
acceptable but optional. It corrects nothing false.

**Wording leg: UPHELD.** The canonical value is min(f, r) for ntHash1 and f + r mod 2⁶⁴ for ntHash2, so it is not
GF(2)-linear. Endorse "hold for ntHash's per-strand hashes, which are XOR-additive over positions" (free). I amend the
attacker's "linear over GF(2)" given the seed result above. The constant-time substitution claim itself survives for
the canonical value, because both strands update in O(1).

## F11 (C19): reading of Marçais et al. NARROWED, as filed (the inference leg OVERTURNED)

The summary (`marcais_2024_knonical.md`) has the source saying "can negate" and "may not preserve". The paper's "show
that ... negates" upgrades that. The table shows fractional sampling's desert sizes moving in both directions (63.20 →
38.14 Mb, 24.68 → 57.52 Mb) with SD 34.7 Mb, and the caption excepts fractional from "orders of magnitude larger".
"No visible effect" overreads both: "no systematic effect" is what the source supports.

The attacker's third leg is OVERTURNED. "Fractional sampling has no window guarantee to negate, so the conclusion does
not follow" is the paper's own point. A sketch using fractional sampling loses nothing to canonicalization's
window-guarantee failure, so the conclusion does follow.

**Repair.** Endorse "can negate ... syncmer-style schemes; for fractional sampling, canonicalization had no systematic
effect on desert size, which varied widely in both directions" (free).

---

## Rule-8 record

- The two-set decoder in F5 was seen by L5 alone, and it is correct.
- The rank-computation premise in F10 is L5-only, and it is wrong.
- F1's weight-29 diagnosis appears in at most two other attack files by keyword, and it is UPHELD on re-run.
  Agreement there is not corroboration.

## Verdict table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | re-run: ours = (1−r)^29, ntHash2 = (1−r)^20 to 4 decimals; the test never compares h(s) with h(rc s); repair free with "matches a 20-position seed" |
| F2 | NARROWED | 1/20 is unsupported (break-even 1/20.2 ± 3–6% job variance), not refuted (0.13% at the endpoint); endorse "1/100 and below, break-even near 1/20" (free); reject the unmeasured "1/25" |
| F3 | NARROWED | 16,362 vs 16,512 orbits re-run: 150 homometric merges; complement leg redundant (implied by reverse within ds classes); caption F_{q²} fix free |
| F4 | NARROWED | 2.6e-17 and 1.3e-17 corrections UPHELD (free, also in the kp10 note); "agrees" narrowed, attacker's rewording OVERTURNED (F+R/F·R are a working negative control); 2.5-pair bound is context, not an error |
| F5 | NARROWED | numbers and bytes reproduce like for like; "equal reliability" and two-set oracle decoder UPHELD (oracle point L5-only); the "other side's hash values" repair OVERTURNED as an unmeasured protocol claim |
| F6 | UPHELD | detection depends on the hash (code comment, Minsky); delete "deterministically"; the "bounded by the collision bound" repair is NARROWED to "negligible probability over x, λ" (free) |
| F7 | NARROWED | bit balance pseudo-replicated (~20× iid SD for every hash); P(zero) = 0.66 for "none of 367"; five methods not four; the repeats conclusion follows |
| F8 | NARROWED | no artifact measures orientation; endorse "attributed, not tested" (free) and naming the method in the abstract |
| F9 | NARROWED | 1.847 is the haplotype pair outside repeats, not a fourth pair; endorse "61–63%, 65% between the haplotypes outside repeats" (free) |
| F10 | NARROWED | main leg OVERTURNED: ntHash seeds are not affine (XOR = 0x047f…), so no 62→64 linear map or rank proof; GF(2) wording for canonical value UPHELD, fix free |
| F11 | NARROWED | "can negate" and "no systematic effect" wording UPHELD (free); "conclusion does not follow" OVERTURNED, since no window guarantee is exactly the paper's point |
