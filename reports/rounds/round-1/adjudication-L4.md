# Round 1 adjudication of L4 (provenance)

**Integrity check.** `sha256(attack-L4-provenance.md)` = `4a8545eee90b58942ab2b23fdbd9240f999e0fe079c45af342d094e24ea00b8d`. It matches the sealed digest `4a8545eee90b` in `round.json`. I read the attack first and the claims second, as the brief fixes.

**Prior rulings.** `scripts/round history "strand-symmetric" --kind adjudication` returns a searched absence ("no prior mention ... in 1 round directories"). This is the first round on this paper, so no finding re-files a settled point.

**Method cards.** The cards live in `bench/repos/library/methods/`; `refs/methods/` does not exist in this repo. Two govern this attack:
- `substitution_scoring.md`, decision 3 (F9). It decides only part of F9; see below.
- `reproducibility.md` (F4). An uncommitted check is not evidence that the paper can point to.

No card governs F8, which is algebra. I decided it from the decoder code.

**What I re-ran.** Scratch scripts are in `/home/agent/.claude/jobs/0d0eac3e/tmp/ADJ-L4/`.
- I checked all of the attack's arithmetic.
- For F2, I recomputed the F·R and F+R counts exactly over Z[i] and Z[ω], with no field and no random point.
- For F3, I derived cycles/task-clock for all 51 runs from `launcher.log`.
- I read the following directly:
  - the identity-check code in `bench_e.cpp` and `bench_pack.cpp`;
  - the `bench_spaced_e.cpp` test loop;
  - the ntHash c26bd45 seed warning and the canonical sum;
  - the `recon61.py` decoder contract;
  - the km-peeler and KP-10 JSONs;
  - the Ts/Tv JSON and the MAF header;
  - the Cristea and Das conversions;
  - the Marçais, Minsky, Kazemi and Mohamadi summaries;
  - the dated negative searches in `literature/`.

---

## F1: the bound's number belongs to the other formula. UPHELD. The fix is free.

- 60/q = 2.602e-17 and 62/q = 2.689e-17. The body says "at most 2(k−1)/q ... about 2.7 × 10⁻¹⁷", and 2.60 does not round to 2.7. The source note (`kp7_eisenstein_2026-11-10.md:29-30`) has "≤ about 2k/q ≈ 2.7 × 10⁻¹⁷"; the body changed the formula and kept the number.
- CM31: 2k/p² = 1.344e-17 and 2(k−1)/p² = 1.301e-17. Only 64/p² = 1.388e-17 rounds to 1.4. The figure in `kp7_m31_simd_2026-11-10.md:13` and in Implementation line 65 is wrong under either formula.

**Whole-set check (rule 3).** The abstract's "about 2k/q, 2.7 × 10⁻¹⁷" is internally consistent, and it is a valid, slightly looser bound. It does not need to change. Only the body's number and the CM31 number are wrong.

**Endorsed repair.**
- Body: "about 2.6 × 10⁻¹⁷".
- Implementation, the CM31 source note, `kp10_vs_shibuya:43` and the ledger rows: "about 1.3 × 10⁻¹⁷".

The fix is free.

## F2: F·R merges more than reversal explains at k = 8. UPHELD. The fix is free apart from one small committed run.

**Burnside count.**
- The fixed points under reversal are 4^⌈k/2⌉ = 256 at both k = 7 and k = 8.
- Under reverse-complement there are none at odd k and 4^(k/2) = 256 at even k.
- Under complement alone there are none.
- So the counts are 4,160 at k = 7 and 16,512 at k = 8.

**My exact enumeration.** I computed over the ring itself, so no random collision is possible.
- Fourth-root encoding: F·R gives 4,160 and **16,362**.
- Eisenstein encoding: F·R gives 4,160 and **16,500**.
- F+R gives 729 and 6,561 under both encodings.

So the 150-value excess (fourth-root) and the 12-value excess (Eisenstein) are structural. They are homometric sequences, which share an autocorrelation without being reverses of each other.

**Scope.** The paper's sentence "cannot tell a k-mer from its plain reverse" is true. What is wrong is the implication that reversal is the whole mechanism, and the k = 8 table count does not match it.

**The Eisenstein 729.** It reproduces under Eisenstein. Only the committed artifact is missing.

**Endorsed repair.** Use the attacker's wording ("or from other k-mers with the same autocorrelation") and commit an Eisenstein run of the three functions. The cost is one sub-minute run.

## F3: 2.0 GHz is not the operating clock, and the spaced-seed timing breaks Methods line 73. UPHELD. The fix is free.

**Clock rate.** I derived cycles/task-clock for all 51 runs.
- Compute runs: about 2.5–3.00 GHz.
- The `__load_only__` runs, including run_1 at 2.47 GHz, which the attacker did not list: 2.22–2.47 GHz.
- `meta.tsv` gives `cpu_model ... @ 2.00GHz`.

So Methods' "Xeon Gold 6330 at 2.0 GHz" is a true description of the part. The abstract's "On a 2.0 GHz Xeon, the hash costs 24.7 cycles" invites a cycles-to-time conversion that is about 25–35% off. The note's "core clock cycles at the 2.0 GHz base frequency" (`kp7_eisenstein_2026-11-10.md:59`) is false. I agree that the cycle counts themselves are sound.

**Spaced-seed timing.** `paper_bench_2026-11-11.md:54` and `run.log:344` confirm that run_51 is a single invocation of `bench_spaced_e` using its own wall-clock rates. It contradicts "Every timing in the paper ... twice ... load subtracted".

**Endorsed repair.**
- Write "2.0 GHz base, turbo enabled; cycles are core cycles".
- Add the explicit exception sentence to Methods.

The fix is free.

**Correct minority (rule 8).** No other lens raised the turbo point; grep of the six other attacks for "GHz" finds no hits. This lone critic is correct.

## F4: the bit-for-bit check is not in the job and was never committed. UPHELD. Repair cost: one run.

- `bench_e.cpp` returns at `if (only) return 0;` (line 162) before its identity pairs. Every paper run used `only=`. The pair list has no NT_BATCH vs library pair in any case.
- The only check is in `bench_pack.cpp:276-295`. It compares an XOR of all hashes plus the count, and (position, hash) element by element on the first 2,000,000 k-mers of record 0.
- `results/` contains no output of it. The only "bit-identical" text is the label in `kp7_next_steps_2026-11-10.md:18`, which has no log behind it.

Under `reproducibility.md`, the claim currently rests on no committed evidence. Its wording ("bit for bit") is also stronger than "XOR checksum plus a 2 M prefix". "About twice" (2.07) holds.

**Endorsed repair.**
- Commit `bench_pack`'s identity output on full cs10 at ntHash c26bd45.
- Describe the check as the attacker words it.

Cost: one laila run. Better still, extend the element-by-element comparison to all k-mers, which costs nothing extra in code.

## F5: at exactly 1/20 the model predicts a loss. UPHELD on a minimal defect. The fix is free.

- 9.96 + 28.5/20 = 11.385 against 11.37. The break-even is 1/20.2 (1/20.3 in the source).
- The body (line 106) already says "about 1/20 and below". Only the abstract drops "about".
- The margin is far inside the repeat spread, so this is a wording defect, not a substantive one.

**Endorsed repair.** Put "about" in the abstract. The fix is free. The attacker's alternative ("1/100 and sparser, measured, 1.14×") is also correct and is the stronger statement, but it is optional.

## F6: the 0.747 vs 0.818 comparison. NARROWED, from the attacker's REFUTED.

**What holds.**
- 0.99²⁰ = 0.81791 against 0.81794 measured, and 0.99²⁹ = 0.74717 against 0.74722.
- The mask has weight 20, and its union with its mirror has 29 positions.
- The test (`bench_spaced_e.cpp:139-167`) compares s with a mutated s on the same strand only. It never tests hash(s) = hash(rc(s)).
- The ntHash2 warning is real (`src/seed.cpp`).
- The paper never tells the reader that 0.747 < 0.818 is a sensitivity cost, which is the figure of merit for a spaced seed.

**What does not hold.** The attack says the "exact double-stranded value" is "defined to match our construction". That is wrong in substance. A natural strand-consistent spaced hash must agree on s and rc(s). rc(s)'s care bases are s's mirrored positions, so such a hash depends on mask ∪ mirror(mask). The 29-position fraction is therefore the correct reference for strand consistency, not a tautology chosen to flatter. A strand-consistent hash on the intersection only is possible, but that is a different seed.

So the sentence is not a false claim. It is a true comparison stated without its cost and without a direct test of the consistency it implies.

**Repair.** I endorse the attacker's wording that states the cost (29 of 31 positions). Qualify "which is not strand-consistent" with "as its own warning states" until the direct test is run. The direct hash(s) = hash(rc(s)) test for both methods costs one short run, and it is what makes the juxtaposition with the warning earn its place.

## F7: "61 percent of it spurious" and "equal reliability". UPHELD on both defects. OVERTURNED on one part of the repair.

**The 61%.** The returned set is 3,238,094 true differences plus 1,986,807 spurious ones. The spurious share is 1,986,807 / 5,224,901 = **38.0%**; 61% is the ratio against the true difference. "Superset" is also false, because the method missed 7 true differences in all 5 trials (`full_kmp_ext_giga.json`).

**Equal reliability.**
- `full_kp10_eis_giga.json` records `rounds_max 2` and `splits 1275`. The bytes check: (65,536 + 2 × 1,275) × 67 × 8 = 36,494,096 exactly.
- The k-mer IBLT's 30/30 is one-shot on `parts [0,1,2,3]` (D = 202,557) and sized at 1.05 × the oracle D.
- The source itself says the whole-genome IBLT figure is extrapolated ("66 MB vs 36.5 MB extrapolated").

The regimes are not equal. One scope correction: KP-10's (b = 16, m = 64) was also chosen from a sweep on this pair (`kp10_eisenstein_b*_m*.json`). So knowledge of D is not one-sided, and the attacker's "sized from the oracle true D" should not be used as an asymmetry.

**Repair.**
- The 38% / 61% / missed-7 wording is endorsed and free.
- The attacker's replacement for "of equal reliability" is endorsed and free.
- I **OVERTURN** the claim that "the IBLT at n = D with one retry (19.6 B/diff, 15/30) is the like-for-like comparison". A 15/30 one-shot table with one retry succeeds about 75% of the time, not 100%. Its retry bytes are not counted, so its expected size is at least 19.6 × 1.5 ≈ 29 B/diff. That is neither equal reliability nor a lower cost.
- Do not install that comparison. If a like-for-like figure is wanted, measure the IBLT at whatever ratio gives 30/30 on all 64 partitions. That is a measurement cost, not free.

## F8: "detected deterministically". NARROWED on the diagnosis. OVERTURNED on the repair's number.

**The decoder contract** (`recon61.py:1-12`):
- each side flags its own set members that are roots of the Padé numerator or denominator;
- the flagged elements must reproduce the c extra sums;
- the flagged counts must equal δ.

**A Vandermonde argument shows when detection is deterministic.** A wrong acceptance needs the signed set (true differences minus flagged elements) to have vanishing power sums p₀ … p_{m+c}. Suppose the hashes are distinct and nonzero. Then a nonzero signed set of support t ≤ m+c+1 cannot have p₀ … p_{t−1} = 0.

- If every flagged element is a true difference, the support is at least n − m. Overflow is then detected **deterministically** for every bucket with n ≤ 2m + c + 1 differences: n ≤ 131 at m = 64, c = 2.
- With about 49 differences per bucket on average, essentially all overflows in the GigaDB run were in this regime.

**Where the attack is right.** Detection is not deterministic in general. It can fail in three cases:
- when a non-difference element is accidentally a root;
- when two elements share a hash;
- when n exceeds 2m + c + 1.

Each failure needs a coincidence over the random (x, λ). So "deterministically" (Discussion) and the unqualified "is detected rather than resolved wrongly" (Methods) overstate.

**Where the attack is wrong.**
- **The citation.** It cites Minsky's probabilistic check as "the same". That check tests a rational function at random evaluation points (unknown-bound protocol). Here the extra sums are fixed functions of the hash values, which is a different mechanism.
- **The bound.** Its proposed bound "(m + c)·2k/q per bucket" is unproven. The decoder's output depends on the random values, so the failure event is not one fixed polynomial identity, and a union over the possible flagged subsets is far larger.

**Endorsed repair.** "an overflowing bucket is detected deterministically when it holds at most 2m + c + 1 differences and the hashes in it are distinct, and otherwise except with negligible probability over the random hash". Alternatively, drop "deterministically" and add "except with negligible probability over the random hash". The fix is free.

Do not print the attacker's bound. Several lenses filed a "deterministic" point; on the bound, this ruling binds only L4.

## F9: the Ts/Tv range and the alignment description. NARROWED. The attacker's repair is OVERTURNED in part.

**What holds.**
- The six values in the JSON are:
  - all columns: 1.585, 1.584 and 1.669, or 61.3%, 61.3% and 62.5%;
  - uppercase only: 1.780, 1.769 and 1.847, or 64.0%, 63.9% and 64.9%.
- "61 to 65" spans all six, so the arithmetic is right.
- KOMPa vs KOMPb is projected through cs10 coordinates and is not a pairwise alignment (`kp1_tstv.py` docstring).
- The MAFs are `scoring=blastz`.

**The card.** `substitution_scoring.md` decision 3 is about estimating a matrix M from alignments made under M. It applies by analogy only, because the paper reports a descriptive spectrum and does not estimate a scoring matrix. A one-clause caveat is reasonable. The card does not make the figure circular.

**Whole-set check (rule 3) against the attacker's repair.** The proposed sentence says 61% (1.58) "rising to 65 percent (1.85) between two haplotypes outside soft-masked repeats". That attributes the rise to the haplotype comparison, but cs10 vs KOMPa outside repeats is already at **64.0% (1.78)**. Almost all of the rise comes from excluding soft-masked columns, not from the pair. The proposed text would publish a new misattribution.

**Endorsed repair.** "Between cs10 and two other assemblies, and between two haplotypes of one plant compared at shared cs10 positions, 61 to 63 percent of aligned substitutions are transitions over all columns (Ts/Tv 1.58–1.67) and 64 to 65 percent outside soft-masked repeats (1.77–1.85); the alignments were made with lastz's transition-favouring scoring." The fix is free.

## F10: the Marçais reading. UPHELD. The fix is free.

- The summary quotes "can negate" (lines 89-90) and "may not preserve" (line 84).
- Its rebuilt Table 2 shows fractional sampling shrinking on chromosome 1 (63.20 → 38.14 Mb) and growing on chromosome 10 (24.68 → 57.52 Mb), with the paper's own SD of 34.7 Mb.
- The copy read was the preprint.

"Negates" and "no visible effect" both harden the source. The attacker's wording is endorsed, together with citing the preprint or checking the journal version. The wording fix is free; checking the journal version is a single read.

## F11: GF(2)-linearity and the attribution. UPHELD. The fix is free.

- ntHash c26bd45 `internal.hpp`: `canonical(fwd, rev) { return fwd + rev; }`. This is integer addition mod 2⁶⁴, which is not GF(2)-linear. Only the per-strand XOR-of-rotations values are.
- The Kazemi summary quotes ntHash2 replacing "the minimum of the hashes, which was used in the previous version" (lines 95-97). So combining strands is ntHash's idea, and the sum is ntHash2's.

The attacker's two replacement sentences are endorsed. The fix is free.

## F12: two Supplement descriptions. NARROWED.

**Upheld.**
- The de Bruijn test ran on "the first 20 Mb of cs10" (`ideas_eisenstein_2026-11-10.md:77`), not on "a cs10 chromosome's graph".
- 16,316,892 / 38.54 M = 0.42 of |A|, while the share of the union is 16.32 / 46.74 = **34.9%**. The wording "differ in 42 percent of their 31-mers" invites the union reading.

**Overturned.** The machine-dependence leg. "Four to six times" traces exactly to the laila table (5.9×, 4.3× and 4.2× at k = 31, 47 and 63). A local 1.7× from an unspecified machine does not make a laila-sourced number a provenance defect.

**Endorsed repair.** The attacker's two wordings. The fix is free.

## F13: the C7 literature readings. NARROWED.

**Upheld: (a), (b) and (e).**
- **(a)** Das's sentence (conversion lines 131-132) says "Purines -G and pyrimidines-C have identical imaginary but opposite real parts". At θ = π/3, G = ½ + j√3/2 and C = −½ − j√3/2, so the imaginary parts are opposite. The paper's paraphrase silently repairs it.
- **(b)** Cristea's purine/pyrimidine remark sits under equation (3), the conjugation map, at conversion lines 777-783.
- **(e)** The paper cites `das_2019_kaiser-exon` as a source while saying it is unread. The "later comparison" is `gulocak_2022_mapping-effect` (its ref. 22, per `eisenstein_wantlist.md:93-96`), which is held in the library and not cited.

**Narrowed to free wording: (c) and (d).**
- **(c)** Equation (4) is (−1 + j)/√2 · √2. It is the θ = π/4 point up to a positive scale, which is immaterial to the reading. Adding "up to scale" is free; the reading is not wrong.
- **(d)** "reviewed by Das and colleagues" can grammatically attach to "the trigonometric family", which das_2021 does review. The ambiguity is real, and reattaching the citation is free.

The whole repair is free, apart from any decision to fetch Das 2019.

## F14: "tens of gaps per seed". UPHELD. The fix is free.

The paper's own source, `kp7_stats_seeds_eisenstein_2026-11-11.md:19`, gives ntHash "341 (8–96)", so one seed left 8. The definition of "in a repeat" (at least half the spanned bases soft-masked) is in `kp7_next_steps:60` and missing from the paper, which matters when 47.7% of cs10 is soft-masked.

The attacker's wording "up to tens (8 to 96)" and the added definition are endorsed. The fix is free. This point was also raised by L2.

## F15: the absence claims omit their dated searches. UPHELD. The fix is free.

- The dated OpenAlex and web searches exist: `literature/wantlist.md:69-83` and `eisenstein_wantlist.md:57-83`, both dated 2026-11-10.
- `benbassat_2014_incremental-hashing` and `chen_2026_hash-functions` are held in the library, and the paper cites neither (grep count 0).

This partly restates what the claim list already notes for C5 and C6. Its added value is pointing at the searches that exist. The dated qualifier is endorsed. "Cite or dismiss chen_2026" is endorsed and costs one read. I agree with the attacker that pre-emption is L6's question.

## The attacker's notes on C3, C4, C8, C17 and C20

- **C3 (survived).** I concur; the recurrence derivation is correct.
- **C4.** UPHELD. "About twice as many multiplications" disagrees with the project's own count of 3 against 7–8 (2.3–2.7×). The cycle ratio is 1.92. "About twice the cycles" is endorsed. The fix is free.
- **C8 (survived).** I concur: |1−ω|² = 3, |1+ω|² = 1, |2|² = 4.
- **C17.**
  - The 137 KB = 2⁸ × 67 × 8 = 137,216 traces.
  - On the cause leg, the attack is **OVERTURNED**. The orientation explanation is very nearly deductive, not merely "inferred". The canonical k-mer difference is only 7,931, so the two assemblies have the same sequence content up to strand. Extended syncmers are read off each assembly as written (and upper-cased), so 48 M differences of about 55 M per side can only come from strand, with contig junctions contributing a negligible share.
  - The attacker should have reasoned from the 7,931 before calling it an inference.
  - Counting reversed contigs would make it a measurement and is cheap, but the text is not wrong.
- **C20.** UPHELD on wording. "Block sizes from 1 to 4,096" describes the discrete sets {1, 3, 5, 64, 4096} and {1, 3, 64, 4096} (`test_kp7e.c:149,197`). The parenthetical listing is endorsed. The fix is free.

---

## Verdict table

| finding | verdict | why |
|---|---|---|
| F1 | UPHELD | 60/q = 2.60e-17 ≠ 2.7; CM31 is 1.30–1.34e-17, not 1.4; the abstract's 2k/q is fine; the fix is free |
| F2 | UPHELD | exact enumeration gives 16,362 (fourth-root) and 16,500 (Eisenstein) against 16,512 from reversal alone; homometric merges; the fix is free plus one tiny run |
| F3 | UPHELD | runs at 2.22–3.00 GHz (turbo); the spaced-seed timing is one wall-clock run; cycle counts are sound; the fix is free; lone correct critic |
| F4 | UPHELD | identity check not run in the paper's job (`only=` returns first), checksum plus 2 M prefix only, no output committed; repair costs one run |
| F5 | UPHELD | the model gives break-even 1/20.2; the abstract lacks the body's "about"; minimal; the fix is free |
| F6 | NARROWED | numbers and the missing direct rc test hold; the "exact" reference is the correct strand-consistency target, not a construction artifact; state the cost |
| F7 | UPHELD | 38% not 61% and not a superset; regimes unequal; the attacker's "IBLT with retry is like-for-like" is OVERTURNED |
| F8 | NARROWED | detection is deterministic for overflows up to 2m+c+1 with distinct hashes and probabilistic only via coincidences; the attacker's bound and Minsky analogy are OVERTURNED |
| F9 | NARROWED | projection and range composition hold; the card applies by analogy only; the attacker's fix misattributes the rise to haplotypes (cs10–KOMPa uppercase is already 64%) |
| F10 | UPHELD | the source says "can negate", and fractional sampling moves both ways with SD 34.7 Mb; the fix is free |
| F11 | UPHELD | the canonical value is the integer sum fwd+rev; ntHash 1 already combined strands by minimum; the fix is free |
| F12 | NARROWED | 20 Mb and 42%-of-A (35% of the union) hold; the machine-dependence leg is overturned (4–6× traces to laila) |
| F13 | NARROWED | (a), (b) and (e) hold against the PDFs; (c) and (d) are scale and attachment wording, not misreadings |
| F14 | UPHELD | the paper's own source gives ntHash 8–96 per seed; define "in a repeat"; the fix is free |
| F15 | UPHELD | the dated searches exist but are uncited; benbassat and chen_2026 are held and uncited; the fix is free |
| C3 | UPHELD | concur SURVIVED; the recurrence re-derived |
| C4 | UPHELD | "about twice the multiplications" against a counted 2.3–2.7×; say "cycles"; the fix is free |
| C8 | UPHELD | concur SURVIVED; the distances are 1, 3 and 4 |
| C17 | NARROWED | 137 KB traces; the orientation cause follows from the 7,931 canonical difference, so the attacker's "inferred" leg is overturned |
| C20 | UPHELD | "1 to 4,096" is a five-point set; list it; the fix is free |
