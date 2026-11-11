# Round 1 adjudication of L1 (correctness)

**Seal check.** `sha256(attack-L1-correctness.md)` = `1632e431b0260f60…`, which matches the sealed `1632e431b026`. The attack read is the one sealed.

**Order.** I read the attack first and then the claims and their evidence, as the brief fixes.

**History.** I ran `scripts/round history` for "strand-symmetric", "overflow", "spaced", "Ts/Tv" and "bit for bit". It found no prior adjudication. This is round 1 for this repository, so no finding re-files a settled point.

**Method cards.** `refs/methods/` does not exist in this repository. No card governs Schwartz–Zippel, power-sum decoding or benchmark protocol, so every dispute below is ruled on the code and the data.

**Scratch.** `/home/agent/.claude/jobs/0d0eac3e/tmp/ADJ-L1/` holds copies of the attacker's scripts, which I re-ran, and `bp.out`, a fresh `bench_pack` run.

---

## F1: λ derived from x through one splitmix64 seed (C1, C18)

**Premise.** Re-run, and it holds.
- `hash/kp7e.c:59-61` draws `x = rand_q(&s); lam = rand_q(&s)` from one 64-bit state.
- `hash/kp7m31.c:83-84` does the same.
- `sm.py` inverts the mixer for seed `0x5eedE`. Given x alone it finds 8 candidate λ, and the true λ is among them.
- A 64-bit seed cannot supply a uniform pair over F_q², which has about 2^122 points.

So the Schwartz–Zippel premise, x and λ independent and uniform, is not what the code implements. For a fixed x with a(x) ≠ 0, exactly one λ collides, and whether the derived λ lands on it is a question about splitmix, not about polynomials. The bound is therefore heuristic for the code as written. The paper's "drawn independently of it" refers to the input, which does hold; independence of x and λ from each other is the unstated gap.

**Scope.** "The polynomial argument does not cover … every reported number" is OVERTURNED.
- No reported number is a measurement of the bound.
- The timings do not depend on it.
- The exhaustive k = 7/8 checks support their real inference, that no two classes have identically equal S, at any single draw, whatever the sampler. An identical pair would collide at every (x, λ).

The attacker should have separated the *theorem* (unaffected) and the *experiments* (unaffected) from the *implementation's claim to inherit the theorem* (affected).

**C18 leg.** NARROWED. The experiments do use fixed seeds, so "a collision in one run occurs in every run" is true of these runs too. The paper's contrast is between designs, where the seed is a parameter here and baked in for ntHash. The attacker's second sentence ("The experiments use fixed seeds; the bound applies to a deployment that draws them at random") fixes this, and it is free.

**Repair.**
- The text fix is endorsed and free: state that the bound assumes x and λ independent and uniform over F_q, and that the experiments use fixed seeds.
- Drawing x and λ from ≥ 122 bits of fresh entropy is correct but **not free**. It changes the code in `kp7e.c` and `kp7m31.c`, every hash value, and any test that pins values, so it needs a test re-run. Timings should not move. Until it is made, the paper should say the implementation seeds both from one 64-bit value.

**Correct minority.** Only L1 filed this. No other attack mentions the shared splitmix derivation.

**Verdict: NARROWED.** The implementation does not meet the bound's premise. The theorem, the exhaustive checks and the timings are untouched.

## F2: 2.7×10⁻¹⁷ is 2k/q, not 2(k−1)/q; CM31 1.4×10⁻¹⁷ is neither (C1)

**Recomputed.**

| bound | over q = 2^61 − 1 | over p² = (2^31 − 1)² |
|---|---|---|
| 60/… = 2(k−1) | 2.602e-17 | 1.301e-17 |
| 62/… = 2k | 2.689e-17 | 1.344e-17 |
| 64/… | — | 1.388e-17 |

- `paper.md:55` says "at most 2(k−1)/q … about 2.7×10⁻¹⁷", and that pairing is wrong in the printed digit.
- `paper.md:65` says "about 1.4×10⁻¹⁷", which is neither 2(k−1)/p² nor 2k/p². It traces to `kp7_m31_simd_2026-11-10.md:12`, which writes "2k/p² ≈ 1.4e-17", itself a slip.
- The degree argument holds. The λ² terms cancel, so the total degree is max(k, 2k−2) = 2k−2.

**Where the attack is weaker than it reads.** The abstract's "about 2k/q, 2.7×10⁻¹⁷" is internally consistent, and 2k/q is a valid, looser upper bound. Changing the abstract is a consistency choice, not a correction.

**Whole set (rule 3).** The 1.4 also appears in:
- `results/kp10_vs_shibuya_2026-11-11.md:43` and `:129`
- `results/kp7_m31_simd_2026-11-10.md:12`
- `manuscript/ledger.md:17`

All must change together.

**Repair is free.** Methods: "2(k−1)/q, about 2.6×10⁻¹⁷". Implementation: "about 1.3×10⁻¹⁷". Alternatively, keep 2.7 and 1.34 by stating the 2k bound everywhere.

**Verdict: UPHELD.**

## F3: "determines the unordered pair … the double-stranded k-mer" (C1)

The paper's sentence opens "As a polynomial in λ" and the next sentence moves explicitly to "a nonzero polynomial in (x, λ)". Read as polynomials, which the text signals, it is correct. The attacker is right that it can also be read as values, where it is false, and the rewrite removes that reading.

"As the factorization argument requires" is loose. The argument predicts distinct polynomials; it requires nothing of one draw. The attacker's replacement ("so no pair of classes has identically equal S") says what the check shows.

I reproduced 8,192 / 32,896 and found F injective on all 4^k k-mers (`exh.py`). Nothing in the paper is false here.

**Both rewrites are endorsed, and the fix is free.**

**Verdict: NARROWED.** This is a clarity defect, not an error.

## F4: overflow detection is not deterministic (C16)

**Re-run.** `ovf.py` calls the paper's own `recon/recon61.py:decode` with m = 64, c = 2, on a 140-difference bucket built from 70th roots of unity. All 66 power-sum differences are zero, and `decode` returns `(0, 0)`, not `None`. "Deterministically" is false as mathematics.

**The attack understates its case.** In `reconcile()` (`recon61.py:157`) a bucket with equal counts and identical sketches is skipped before `decode` is ever called. The code's own comment says this holds "with probability 1 − O(m/q)". The constructed bucket would therefore pass silently in the real protocol, and the code already concedes probabilistic behaviour.

**Repair.**
- Deleting "deterministically" (Discussion) is endorsed, and the fix is free.
- "Except with negligible probability over the random seed" (Methods) is NARROWED. No bound is derived anywhere, and with fixed seeds (F1) it is not a probability over anything the user controls. Use "detected with high probability", or give a derived bound.
- Whole set: `results/kp10_vs_shibuya_2026-11-11.md:97`, `:124` and `:133` all say "deterministic" and need the same edit. That is free.

**Verdict: UPHELD.**

## F5: the symmetric spaced seed depends on 29 positions, not 20 (C14)

**Verified.**
- `spaced_e.h:4` defines R_m = F_m(rc(window)).
- The mask `1101100111010111011100101110101` has weight 20, and its union with its mirror has 29 positions.
- Measured against predicted: 0.74722 vs 0.99²⁹ = 0.74717, and 0.81794 vs 0.99²⁰ = 0.81791.
- The "exact" reference in `bench_spaced_e.cpp:7-8,149` is defined over the care positions of *both strands*. So "the exact double-stranded value" is accurate by its own definition, but the paper never says that definition costs 9 don't-care positions.

The loss is inherent to any strand-symmetric hash of a masked window. It is a cost the paper should state, not a bug.

**Repair.**
- The first fix sentence ("depends on the mask and its mirror, 29 positions … 0.747 = 0.99²⁹") is endorsed, and it is free.
- The second sentence is OVERTURNED as worded. "SeedNtHash … reads only the forward mask's 20 positions, which is why it is not strand-consistent" installs a causal claim the attacker admits was never tested, since h(w) and h(rc w) were not compared. Say instead "0.818 = 0.99²⁰, consistent with reading only the forward mask's 20 positions". Better still, test h(w) = h(rc w) directly; that costs one run.

**Verdict: NARROWED.** The defect is upheld; the prescribed causal repair is not.

## F6: 1/20 is just above the break-even of 1/20.3 (C11)

**Recomputed.** 9.96 + 28.48 × 0.05 = 11.384 against 11.37 for batched ntHash.
- At exactly 1/20 the sampler ties ntHash, or is 0.1% slower, inside the 1.3% repeat spread.
- I checked that no other mode rescues 1/20: branchy, interpolated between 1/100 and 1/16, gives about 12.5.
- The extrapolated line predicts 10.245 at 1/100 against a measured 10.22, so the line itself is sound.
- Branchy (inline) slope: (10.01 − 9.40)/0.009 ≈ 68 cycles per kept k-mer. So "a kept k-mer costs about 28.5 cycles" is the collected-mode figure only.

**Both fixes are endorsed, and the fix is free.** Abstract: "below about 1/20". Results: "with collected evaluation".

**Verdict: UPHELD.**

## F7: "61 percent of it spurious" (C15)

The source says "61% of D". Of the superset it is 1,986,807 / (3,238,094 + 1,986,807) = **38.0%**. The paper's "of it" points at the superset, so the sentence is wrong.

**Repair.** The attacker's fix is correct and free. One alternative preserves comparability with Shibuya's "about 3%", which the results file sets against the same per-D basis: "spurious k-mers numbering 61 percent of the true difference". Either is free.

**The rounds point.** The KP-10 row records 2 rounds. Disclosing "two rounds against one" next to "equal reliability" is endorsed, and it is free.

**Verdict: UPHELD.**

## F8: the fourth-root hash takes 8 multiplications against 3 (C4)

- `kp7.c:49` (`fmul_real`) costs 2 multiplications per strand.
- `kp7_sym` uses `fmul_p`, and `karatsuba = 0` by default (`kp7.c:117`). `bench_e.cpp` never sets it.
- The measured configuration is therefore 4 + 4 = **8** multiplications against 3, a ratio of 2.7, not "about twice".
- The cycle ratio, 47.37/24.69 = 1.92, is about twice. The congruences hold.

**Repair is endorsed with one sharpening, and it is free.** "8 multiplications per k-mer (7 with Karatsuba) against 3, at about twice the cycles."

**Verdict: UPHELD.**

## F9: F·R collapses beyond reversal; the caption names the wrong field (C2)

**Re-run** (`exh.py`, `fr.py`):
- At k = 8 there are 16,512 orbits under {reverse, complement, reverse complement}, but F·R takes 16,362 values.
- There are 146 extra collision groups, identical at three draws of x, so they are polynomial identities. One example is AAACCCCT against CAAATCCC.
- Under the Eisenstein encoding, F·R takes 16,500.
- `prototype/kp7_symmetric_check.py:1` computes over Z[i]/(q), q = 2³¹ − 1, so the caption's "q = 2³¹ − 1" implies a field the computation did not use.
- 729 and 6,561 reproduce.

The text is true but incomplete; the caption is wrong.

**Both fixes are endorsed, and the fix is free.**

**Verdict: UPHELD.**

## F10: ntHash's canonical output is not GF(2)-linear (C18)

`hash/ntbatch.h` emits `f + r` in integer arithmetic mod 2⁶⁴, which is not GF(2)-linear. The strand hashes, built from rotations and XOR, are linear. The paper's properties still hold, as the attacker says himself, through the (f, r) pair.

**Rewording endorsed, and the fix is free.**

**Verdict: UPHELD**, scoped as the attacker scoped it.

## F11: the bit-for-bit check: where it is and what it covers (C12)

**Verified.**
- The paper job runs only `bench_e` with `only=`, which returns before any identity check (`bench_e.cpp`), and `bench_e`'s identity pairs never include ntbatch against the library.
- The check lives in `bench_pack.cpp:274-296`: an XOR-plus-count over all k-mers, and element-wise (position, hash) pairs only for the first 2,000,000 k-mers of record 1.

**Re-run.** The attacker did not test whether the check passes, so I ran the repository's built `hash/bench_pack` on `cs10_head6m.fa` (`bp.out`).
- "ntHash batched identical to ntHash canonical: yes"
- "positions and hashes identical on the first 2000000 k-mers of NC_044371.1: yes"

No counterexample to "bit for bit" exists. A gap the attacker missed: **no output of the genome-wide check is recorded anywhere in `results/`**, and the ledger has no row for it. "About twice" holds (2.07).

**Repair.** The scoped wording is endorsed, and it is free. Recording one full-cs10 `bench_pack` identity output in `results/` costs one run and is owed.

**Verdict: NARROWED.** The claim's scope is overstated; the claim is not refuted.

## F12: the upper Ts/Tv end comes from a repeat-masked subset (C9)

`results/ideas_eisenstein_2026-11-10.md`:
- The three full alignments give Ts/Tv 1.585, 1.584 and 1.669, which is 61.3–62.5% transitions.
- 1.847 (64.9%) is "uppercase only" for the haplotype pair alone.
- The ledger's "1.58–1.85 / 61.3–64.9%" mixes strata without saying so.

**Repair endorsed, with one adjustment, and it is free.** Write "61 to 62.5 percent (Ts/Tv 1.58 to 1.67) … 65 percent outside soft-masked repeats, in one pair". Rounding 62.5 to 63 is acceptable but not needed.

**Verdict: UPHELD.**

---

**Minority record.** F1 was filed by L1 alone and is narrowed, not discarded: its implementation premise is the most consequential lone finding in this attack. The overflow point (F4) and the 1/20.3 point (F6) appear in several other attacks' texts. Same-family agreement is not corroboration; both were re-run here.

**Free repairs owed** (wording only, no re-measurement): F2, F3, F4 (wording), F5 (first sentence only), F6, F7, F8, F9, F10, F11 (wording), F12, and the F1 text sentence.

**Not free:**
- F1 code reseeding: code change, new hash values, test re-run.
- F5 direct strand test of SeedNtHash: one run.
- F11 recorded full-genome identity output: one run.

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | x and λ do come from one 64-bit seed, so the code does not inherit the bound; "every reported number" is wrong, since the exhaustive checks and timings do not depend on it |
| F2 | UPHELD | 2(k−1)/q = 2.6e-17 and the CM31 1.4e-17 matches neither bound; the abstract's 2k/q is valid; free fix across paper, two results files and the ledger |
| F3 | NARROWED | the text says "as a polynomial" and is correct read that way; the rewrite is a free clarity fix |
| F4 | UPHELD | the counterexample decodes 140 differences as zero, and reconcile() skips equal sketches anyway; delete "deterministically" (free); "negligible" is unquantified |
| F5 | NARROWED | the 29-position dependence is real and undisclosed; the "which is why" causal repair is untested |
| F6 | UPHELD | at 1/20 the sampler ties ntHash (11.38 vs 11.37); 28.5 cycles is collected mode only; free fix |
| F7 | UPHELD | spurious share of the superset is 38.0%, not 61%; disclose two rounds against one; free |
| F8 | UPHELD | measured configuration is non-Karatsuba, 8 multiplications against 3; cycles about 2×; free |
| F9 | UPHELD | 146 homometric collision groups beyond reversal at k = 8; the caption's field is wrong; free |
| F10 | UPHELD | the canonical f + r mod 2⁶⁴ is not GF(2)-linear; the strand hashes are; free rewording |
| F11 | NARROWED | the check is a checksum plus a 2M prefix in another binary, but it passes on re-run; scope wording is free; a recorded run is owed |
| F12 | UPHELD | 1.847 / 64.9% is the uppercase-only subset of one pair; free rewording |
