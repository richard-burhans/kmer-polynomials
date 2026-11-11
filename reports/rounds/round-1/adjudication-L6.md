# Round 1 adjudication of L6 (prior-art)

Adjudicator: independent skeptic, 2026-11-11. The attack file `attack-L6-prior-art.md` hashes to `5451ee912f9f`, which matches the sealed digest. `ATTACKED-paper.md` is identical to `manuscript/paper.md` (`diff -q` is silent).

**Prior rulings.** `scripts/round history "overflow" --kind adjudication` reports a searched absence, and round 1 is still open. This is round 1, so no earlier ruling bears on any finding.

**Method cards.** The attacker reports that the repo has no `refs/methods/`, and none of these findings is a statistical-method dispute. No card governs any of them.

**My own checks (2026-11-11).**
- The full-text conversions in `library/papers_md/`: minsky_2003, naumenko_2019, mohamadi_2016, kazemi_2022, marcais_2024.
- The summaries of dodis_2008, minsky_2002, katoh_2002 and das_2017.
- The ntHash checkout at c26bd45: `src/internal.hpp`, `src/seed.cpp`.
- The decoder the paper's sketch numbers came from: `recon/recon61.py`.
- WebFetch of `raw.githubusercontent.com/sipa/minisketch/master/include/minisketch.h`.
- One WebSearch for a published strand-symmetric product hash. It found nothing; the only hits were canonical-k-mer documentation from OBIKmers.

## F1: overflow detection is not deterministic (C16)

**Diagnosis: UPHELD.** I tested the premise against the code, not only against the citations.

The decoder (`recon61.py:100-130`) accepts a bucket only when three things hold:
1. Each side's flagged elements are roots of P or Q and are drawn from **its own** set.
2. The flagged elements reproduce the c extra power sums S_{m+1..m+c}.
3. The flagged counts differ by delta.

Condition 1 makes this stronger than Minisketch's check, because a false decode must be built from real elements of each side's bucket. It does not make the check deterministic. A wrong acceptance needs the signed multiset (true difference − accepted difference) to have zero power sums at orders m+1 … m+c. With nonzero distinct values, a Vandermonde argument excludes that only when its support has at most c elements. Any overflow larger than that is caught only because hash values are unlikely to satisfy c polynomial equations by chance. That is a probabilistic event over x and λ.

There is a deterministic case: the `abs(delta) > mbar` early return. It covers only overflows visible in the size difference.

The sources say the same. Minsky 2003 (conversion lines 525-531) says such a test "must be probabilistic in nature in order to be efficient". Erlay (lines 532-537) says decoding is "very likely to fail" beyond capacity. minisketch.h says an oversized sketch may "still decode, but the result will be nonsense", and it sizes capacity to a 2^-fpbits false-positive rate.

**Repair: NARROWED.**
- **Endorsed, free:** delete "deterministically" from the Discussion. Replace "detected rather than resolved wrongly" with "detected with high probability: a wrong decode must reproduce the c extra sums by chance".
- **Not endorsed: the attacker's "about q^{-c}".** That figure is a uniform-random-hash heuristic, of the kind behind Minisketch's fpbits. Minsky's bound ((D−1)/|E|)^k is for independently drawn random evaluation points, which this scheme does not use, so it does not transfer. A bound that is proven for this hash would need Schwartz–Zippel on polynomials in (x, λ) of degree about 2(k−1)(m+c). It would also need a union over candidate subsets of each bucket, and that union is not obviously small. Printing "q^{-c}" would put an unproven number into a paper whose selling point is proof.
- If the authors want a number, they must derive it or label it heuristic. That costs a derivation, not a re-measurement.

## F2: power-sum sketch attribution (C15, C16)

**UPHELD on the defect; NARROWED on scope.**

The defect is real:
- Minsky 2003 transmits evaluations of the characteristic polynomial and recovers the difference by rational interpolation. A grep for "power sum" or "Newton" in the conversion returns nothing.
- The power-sum sketch is PinSketch (dodis_2008, Construction 6), deployed by Erlay.
- Split-and-resend is Minsky 2002's PARTITION-RECON and Erlay's bisection.
- All three are entries in `references.bib`, and each is cited 0 times in `paper.md`.

The scope is narrower than the attack says:
- The paper's decoder is not PinSketch as published. PinSketch works in characteristic 2, sends odd sums only, and cannot tell which side an element came from. The signed P/Q split by Padé, "say which set each came from", is Minsky's rational-function idea transported to power sums. Keeping minsky_2003 for that step, as the attacker proposes, is right.
- The attacker's remark that the "Fourth" contribution should disclose prior art is weak. That contribution claims a *comparison*, not a new sketch.
- The defect does not touch C15's numbers, as the attacker concedes.

**Repair endorsed, free:** cite dodis_2008 and naumenko_2019 at the power-sum sentences (Introduction, line 23; Methods, line 59), and cite minsky_2002 and naumenko_2019 at "split and resent". All of these are bib entries already.

## F3: (F+λ)(R+λ) is a characteristic-polynomial fingerprint (C6)

**NARROWED, as the attacker filed it. I uphold that narrowing.**

The algebra is correct. (λ+F)(λ+R) is the characteristic polynomial of the multiset {F, R} evaluated at −λ. Minsky 2003 (lines 124-129) says this representation "is not novel and was first proposed in [18]" (Lipton, STACS 1990, pp. 207-215, line 897).

The per-strand (k−1)/q bound (Introduction, line 21; Methods, line 45) is stated without a citation. benbassat_2014 is held in the library but is not in the bib. chen_2026 is in the bib and uncited.

What survives:
- The *composition* applied to DNA strand symmetry, with a two-variable Schwartz–Zippel bound over (x, λ), has no located precedent. That covers the attacker's searches 1 and 5, the project's 2026-11-10 OpenAlex queries, and my own search today.
- The paper still lacks the dated search marker the claim itself flags.

**Repair endorsed:**
- One attribution sentence citing minsky_2003, with Lipton named through it. Free.
- A citation for (k−1)/q: chen_2026 is free; benbassat_2014 needs one bib entry.
- A dated search statement. Free.

Citing Lipton 1990 directly would mean acquiring a paper the library does not hold. That has a small cost and is optional.

## F4: ntHash history and "linear over GF(2)" (C18)

**NARROWED. The attacker filed REFUTED; I overturn that scope.**

The two attribution defects are real:
- kazemi_2022 (lines 95-97) says the sum replaced "the minimum of the hashes, which was used in the previous version".
- mohamadi_2016 (lines 85-97) already computes the reverse-complement hash and the canonical form.
- So the Discussion sentence "Strand symmetry by combining the two strand hashes is ntHash2's idea" is wrong: the combining is ntHash 2016's idea, and the sum is ntHash2's change.
- The Introduction's "ntHash … by adding" is true of the ntHash *library* the paper benchmarks (`internal.hpp:27-29`, `return fwd + rev;`). Its citation of mohamadi_2016 for "adding" is misleading, and it contradicts the Discussion.

The GF(2) phrase is imprecise, not false in substance:
- The per-strand roll is XOR of rotated table entries, so it is GF(2)-linear. The listed properties hold because each strand is maintained separately.
- The same is true of this paper's own S, which is not linear either and is maintained through F and R.
- The canonical sum (and the old min) are indeed not GF(2)-linear.

C18's core holds: the seeds are compile-time constants (the attacker concedes this), ntHash has no proven collision bound, and the 8-mer test was not attacked. "REFUTED as written" therefore overstates the finding.

**All three of the attacker's replacement sentences are endorsed and free**, with one change: say "ntHash 2016 took their minimum" rather than "ntHash … takes", since the current ntHash library adds.

## F5: Marçais et al. on fractional sampling (C19)

**NARROWED.**

Upheld on the wording:
- The source hedges: "can negate" and "may not preserve" (lines 84-90).
- Its theory covers context-free methods in general, not only syncmer-style schemes.
- "No visible effect" overstates Table 2. The library's rebuild shows fractional sampling moving both ways (63.20 → 38.14 Mb on chr1, 24.68 → 57.52 Mb on chr10), within an SD of about 30 Mb (line 747-749, verified).

Overturned on the headline "cites them for the opposite":
- Marçais say fractional sampling has the weakest guarantee, with deserts before and after canonicalization. That supports the paper's actual point that long gaps are "not a weakness of either" hash: the deserts belong to the sampling scheme and the repeats, not to the hash or the canonicalization.

**Repair NARROWED.**
- Endorsed, free: "canonicalization can negate the window guarantee of context-free schemes", plus "fractional sampling has no window guarantee and large sketching deserts with or without canonicalization".
- Not endorsed: "so the long gaps we see are expected of any fractional sampler". Marçais measured deserts of at least 75 k-mers at density 1/(k−s+1) on human chromosomes, which is a different regime from our gaps of more than 200 k-mers at our rates. The paper's own causal evidence is the soft-mask overlap (C13). Say "consistent with" rather than "expected of any".

## F6: the ntHash2 warning comes from the code, not the paper (C14)

**NARROWED, as filed. I uphold it.**

The warning is real: `seed.cpp:96-102`, "is not symmetric, reverse-complement hashing will be inconsistent". kazemi_2022 (lines 113-114) states reverse-complement spaced-seed hashing without that restriction. I did not check the ntHash2 supplement, and neither did the attacker; the bib cites only the paper. The sentence is true of the software and mis-sourced.

**Repair endorsed:** cite the library ("ntHash v2.4.0, SeedNtHash"). This costs one software bib entry and no re-measurement.

## F7: the FFT diagonal finder is MAFFT's FFT step (C21)

**NARROWED, as filed. I uphold it.**

The katoh_2002 summary confirms the method: correlation by FFT, the top 20 peaks, and a 30-site sliding window. The project's own `literature/wantlist.md` tags katoh_2002 for "KP-13 … the FFT step", and it is absent from the bib. Sibling bullets in Supplement S1 cite their prior art (Voss, Faloutsos), so this omission is inconsistent with the section's own practice. The supplement claims no novelty explicitly, which is why this is NARROWED rather than more. The numbers were not attacked.

**Repair endorsed:** cite Katoh 2002. One bib entry for a paper already held; nearly free.

## F8: C5 and C7 survived the attack

**C5: SURVIVED, and I agree.**
- The identity and "no fourth-root map" hold: distinct non-opposite fourth roots are orthogonal, so Re(A C̄) = 0.
- The attacker's search is thin, a single WebSearch, but the project's OpenAlex queries back it.
- The owed fix is the dated search marker. Free.

**C7: SURVIVED, and I agree.**
- I checked the das_2017 summary: θ = π/3, the rationale, no distances, and 1218/1251 = 97.36% as the Eq. 8 value.
- The attacker's wording nit is weak. The paper already says "what was later called the trigonometric family", and "fixed the family's angle" is accurate for a paper that introduced θ and set it. The rewording is optional and free.

## Lone or shared points

Every finding here was filed by L6 alone, and I cannot see the other lenses' rulings. The ones I uphold on a single attacker's evidence are F1 (C16 not deterministic), F4's attribution, and F6's code-vs-paper source. They are recorded as lone correct critiques.

## Verdict table

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | not deterministic: UPHELD from the decoder code and three sources (delete "deterministically", free). The prescribed "about q^{-c}" is OVERTURNED: a heuristic, not a bound for this scheme |
| F2 | NARROWED | power sums and split/resend belong to PinSketch, Minsky 2002 and Erlay, and must be cited (free). The signed Padé recovery does rest on Minsky 2003; the "Fourth contribution" complaint is overreach |
| F3 | NARROWED | (F+λ)(R+λ) is the Lipton/Minsky characteristic-polynomial fingerprint and (k−1)/q is uncited. The composition for strand symmetry survives, with no precedent found |
| F4 | NARROWED | "combining strands is ntHash2's idea" is wrong (it is ntHash 2016's) and Introduction/Discussion conflict (free fix). The seed and bound claims hold, and the GF(2) phrase is imprecise rather than false, so not REFUTED |
| F5 | NARROWED | the hedge, scope and "no visible effect" are overstated (free fix). Marçais support rather than contradict "not a weakness of either"; "expected of any fractional sampler" is not endorsed |
| F6 | NARROWED | the warning is in the ntHash v2.4.0 code, not the cited paper; cite the software (one bib entry) |
| F7 | NARROWED | the FFT diagonal finder is MAFFT's step, held and want-listed but uncited; add Katoh 2002 (nearly free) |
| F8 | UPHELD | the attacker's SURVIVED verdicts on C5 and C7 are correct (the claims stand). Owed: C5's dated search marker (free); the C7 wording nit is optional |
