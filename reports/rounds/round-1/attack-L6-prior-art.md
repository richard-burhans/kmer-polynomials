# Round 1 attack, slot L6, lens: prior-art

Attacker: L6 (prior-art). Date: 2026-11-11. Target: `reports/rounds/round-1/ATTACKED-paper.md`.

**Where the literature is.** The brief has no appendix, so I located the sources myself:
- summaries in `../library/papers_md/summaries/<slug>.md`
- full-text conversions in `../library/papers_md/<slug>.md`
- the ntHash source the benchmarks used, in `bench/data/kmer_poly_work/ntHash` (v2.4.0, commit c26bd45)
- the project's own novelty searches, in `literature/wantlist.md` and `literature/eisenstein_wantlist.md`

The project has no `refs/catalog_index.json`. I used `manuscript/references.bib` and the summaries directory instead.

**Method cards.** I found none that govern this repository: it has no `methods/` directory, and none of this lens's claims is a statistical method. **Prior rounds.** This is round 1. `scripts/round history "strand-symmetric"` returns "no prior mention … a searched absence". No prior adjudication applies to any finding below.

**Live searches (all 2026-11-11):**
1. WebSearch: "strand-symmetric canonical k-mer hash polynomial finite field reverse complement provable collision bound". Found generic polynomial-hash bounds only. No strand-symmetric polynomial hash.
2. WebSearch: "multiset fingerprint product (x - a_i) random point Lipton 1990 …". Found the multiset-polynomial fingerprint technique. The year was not confirmed online; Minsky 2003's ref. [18] gives it as STACS 1990.
3. WebSearch: "Minisketch power sums BCH decode more than capacity …". Nothing on Minisketch itself.
4. WebFetch of github.com/sipa/minisketch (README) and raw `include/minisketch.h`. The quotes are in F1.
5. WebSearch: "\"reverse complement\" invariant hash k-mer symmetric function forward reverse hash product sum canonical without minimum". Found k-nonical, the Biostars 8219 thread (combine the two hashes with "any other symmetric function", no proof) and the Wittler canonical MPHF. None pre-empts (F+λ)(R+λ).
6. WebSearch: "complex numerical mapping DNA nucleotides transition transversion distance complex plane angle …". No complex-map transition/transversion criterion found.

## F1 -- the cited source says overflow detection cannot be deterministic (C16)

The paper says "The $c$ extra sums check the result, so a bucket holding more than $m$ differences is detected rather than resolved wrongly". The Discussion adds that the sketch "detected its own overflow deterministically".

Three held or primary sources say otherwise. The first is the paper's own citation.
- **Minsky, Trachtenberg & Zippel 2003** (`library/papers_md/minsky_2003_set-reconciliation.md:529-530`): "Any such test must be probabilistic in nature in order to be efficient". The c extra values are exactly their "extra evaluation points" (lines 556-573). The probability of accepting a wrong result is bounded by a power of (D−1)/|E| (lines 537-541), and that probability is not zero.
- **Naumenko et al. 2019, Erlay** (`library/papers_md/naumenko_2019_erlay.md:532-537`). This is the deployed PinSketch power-sum protocol. It says decoding "will always succeed when the size of the difference … does not exceed the capacity … Otherwise, the procedure is very likely to fail." "Very likely" is not "always".
- **Minisketch** (`include/minisketch.h`, fetched 2026-11-11). It says an oversized sketch "it is possible that it will still decode, but the result will be nonsense". Capacity is chosen so that such a sketch has "a chance no higher than 2^-fpbits of being decoded incorrectly."

The mechanism is easy to see. If a bucket's signed difference D has |D| > m, a false decode needs a signed set D′ with |D′| ≤ m whose first m+c power sums match those of D. In other words, the signed multiset D − D′ must have zero power sums up to order m+c. Vandermonde rules this out only when |D|+|D′| ≤ m+c. Beyond that, it is an event over the random hash values, so its probability is small but not zero.

**Correction.** Replace "detected rather than resolved wrongly" with "detected except with probability about q^{-c} per overflowing bucket (Minsky et al. 2003; cf. Minisketch)". Delete "deterministically" from the Discussion. **C16: REFUTED.**

## F2 -- the power-sum sketch is PinSketch, not Minsky's; the split-and-resend is Minsky–Trachtenberg/Erlay (C15, C16)

The Introduction and Methods credit "Power sums of the hashes … recovered exactly with Newton's identities" to `@minsky_2003_set-reconciliation`. Minsky 2003 sends **evaluations of the characteristic polynomial** and recovers by rational interpolation (summary, Method bullets 1-3). It does not send power sums.

The power-sum sketch is PinSketch, from Dodis, Ostrovsky, Reyzin & Smith 2008, Construction 6 (`summaries/dodis_2008_fuzzy-extractors.md`: "the sketch is the vector of odd power sums"). Erlay and Minisketch deploy it. Both papers are held in the library and are already entries in `manuscript/references.bib`, yet `grep -c "dodis\|naumenko" manuscript/paper.md` returns 0.

The same holds for the paper's overflow handling: "It is then split and resent, at the cost of one more round". This is Minsky & Trachtenberg 2002's partitioned protocol, also in the bib and uncited, and Erlay's bisection fallback (`naumenko_2019_erlay` summary, "Failure handling").

The k-mer application may be new: the project's OpenAlex search of 2026-11-10 is in `literature/wantlist.md`, and search 1 here found nothing either. The sketch construction itself is prior art, and the "Fourth" contribution should say so.

**Correction.** Cite PinSketch (dodis_2008) for the power-sum sketch, and minsky_2002 plus naumenko_2019 for partition/bisection on overflow. Keep minsky_2003 for the characteristic-polynomial idea and its verification bound. **C15: NARROWED** (on attribution only; I did not attack the byte counts).

## F3 -- (F+λ)(R+λ) is the classical characteristic-polynomial fingerprint of the multiset {F, R} (C6, and the claim C1 rests on)

S = (λ+F)(λ+R) is the characteristic polynomial of the two-element multiset {F, R}, evaluated at −λ. The step "factors uniquely … therefore determines the unordered pair" is that polynomial's defining property: its zeros are the multiset.

The paper's own citation says this representation is old. Minsky 2003, §I-C (`papers_md/minsky_2003_set-reconciliation.md:126-130`): "Central to our approach … is the representation of sets by their characteristic polynomials. This representation is not novel and was first proposed in [18] and later by Blum and Kannan [19] as part of a probabilistic set equality test". [18] is R. J. Lipton, "Efficient checking of computations", STACS 1990 (line 897). Using a random evaluation point to test multiset equality is Lipton's fingerprint (search 2, 2026-11-11).

The other half is the per-strand hash, a polynomial over F_q evaluated at a random x. That is Karp–Rabin hashing. Ben-Bassat & Chor 2014 apply it to DNA reads (held: "A number of IH functions are known, most of which use univariate polynomials over finite fields", lines 153-154), and Chen 2026 reviews it. The paper cites neither for the (k−1)/q bound in the Introduction and in "A polynomial hash over a prime field". chen_2026 is in `references.bib` and uncited in the text.

The "Second" contribution is therefore a composition of two classical fingerprints: Karp–Rabin on each strand, then Lipton's multiset fingerprint on the pair. What may be new is the composition, applied to DNA strand symmetry, with base values in F_q via ω. In six searches I found no paper that publishes that composition (searches 1 and 5 here, plus the project's 2026-11-10 OpenAlex queries in `literature/eisenstein_wantlist.md`). Biostars 8219 suggests "any symmetric function" of the two strand hashes but gives no construction or bound. The paper itself carries no dated search marker for this absence.

**Correction.** In the Methods section, add: "(F+λ)(R+λ) is the characteristic-polynomial fingerprint of the multiset {F, R} [Lipton 1990; Minsky 2003], composed with Karp–Rabin hashing of each strand". Cite Karp–Rabin, or Ben-Bassat & Chor 2014 and Chen 2026, for (k−1)/q. Recast the Second contribution as the composition and its application, with a dated search statement. **C6: NARROWED.**

## F4 -- ntHash's history is misattributed, and "linear over GF(2)" is false for its canonical hash (C18)

Three sentences are wrong when checked against the held papers and the code the benchmarks ran.

1. The Introduction says "ntHash rolls such a hash … by adding hashes of the two strands … [@mohamadi_2016_nthash; @kazemi_2022_nthash2]". ntHash (2016) did not add the two strand hashes. ntHash2 says it introduced the sum "as a replacement for the minimum of the hashes, which was used in the previous version" (`papers_md/kazemi_2022_nthash2.md:95-97`).
2. The Discussion says "Strand symmetry by combining the two strand hashes is ntHash2's idea." ntHash 2016 already computes the reverse-complement hash alongside the forward one, "as well as the canonical form … without actually reverse-complementing" (`papers_md/mohamadi_2016_nthash.md:85-97`, eq. 4). Combining the two strand hashes is ntHash 2016's idea. ntHash2's change was to combine them by sum rather than by minimum. The Introduction and the Discussion also contradict each other about who did what.
3. The Discussion says "… hold for ntHash too, which is linear over GF(2)". Only the single-strand hash is GF(2)-linear (rotations and XOR). The canonical hash the paper compares against is the integer sum mod 2^64, `return fwd + rev;` (`ntHash/src/internal.hpp:27-29`), and carries make that non-linear over GF(2). ntHash 2016 used min, which is also non-linear. The extra hashes are `hash_array[0] * (i ^ k * MULTISEED)` (`internal.hpp:115`; mohamadi_2016 line 102, "multiplication, shifting and XOR"). The properties listed may still hold if each strand is updated separately and the two are recombined. The stated reason does not hold.

**Correction.**
- Introduction: "ntHash rolls hashes of both strands in constant time and takes their minimum; ntHash2 takes their sum and extends this to spaced seeds."
- Discussion: "Strand symmetry by combining the two strand hashes is ntHash's idea [2016]."
- Discussion: "… hold for ntHash too, whose per-strand hash is linear over GF(2)."

**C18: REFUTED** (as written). The seed-fixing statement itself holds: the seed tables and MULTISEED are compile-time constants.

## F5 -- Marçais et al. say fractional sampling has the WORST deserts; the paper cites them for the opposite (C19)

The paper says: "Marçais and colleagues show that folding a k-mer with its reverse complement negates the window guarantee of syncmer-style schemes, but it has no visible effect on fractional sampling [@marcais_2024_knonical], which is how a sketch would use this hash." It uses this to support "Long sampling gaps are not a weakness of either."

The source says otherwise on three points.
- **Hedged verb.** The source says "can negate" (preprint line 90) and "may not preserve" (line 84), not "negates". Its theory covers context-free methods in general, Mykkeltveit and Champarnaud sets included, not only syncmer-style schemes (lines 145-147, 279-281).
- **Fractional sampling is the worst case.** Lines 735-737: "The fractional set method has the weakest guarantee, with significant sketching deserts before and after canonicalization. The variance due to the choice of the random permutation … is also high (standard deviation of 34.7 Mb for φ and 28.3 Mb for φc on Chr 1 with k = 31)." Lines 506-509: fractional sampling "has no window guarantee, and leaves one very large strongly connected component … a random set of k-mers does not approximate a set with a window guarantee."
- **"No visible effect" overstates the table.** Table 2's caption excepts fractional from "orders of magnitude larger", but the library's rebuild from PDF p. 8 shows it moving both ways: 63.20 → 38.14 Mb on chr1 and 24.68 → 57.52 Mb on chr10 at k = 31 (`summaries/marcais_2024_knonical.md`). The supportable statement is "no systematic effect, within a run-to-run SD of ~30 Mb".

The paper's sampler is exactly a fractional, hash-threshold method. Marçais is therefore evidence that long gaps ARE a property of this sampling scheme. They come from the de Bruijn SCC meeting repeats, which fits the paper's own soft-mask finding in C13. They do not show that gaps are harmless.

**Correction.** "Marçais and colleagues show that canonicalization can negate the window guarantee of context-free schemes. Fractional sampling has no window guarantee to lose and large sketching deserts with or without canonicalization, so the long gaps we see are expected of any fractional sampler." **C19: NARROWED.**

## F6 -- the ntHash2 "warning" is a library runtime message the cited paper does not contain (C14)

The paper says "ntHash2's spaced-seed hash warns that reverse-complement hashing becomes inconsistent for a non-palindromic mask". The warning exists, but in the code: `ntHash/src/seed.cpp:96-102`, v2.4.0, reads "Seed … is not symmetric, reverse-complement hashing will be inconsistent".

The only ntHash2 source in the bibliography is the paper, and the paper states RC hashing for spaced seeds with no restriction: "ntHash2 can perform reverse-complement hashing without requiring extra iterations by swapping the corresponding indices in the blocks" (`papers_md/kazemi_2022_nthash2.md:113-114`). A reader who follows the citation finds the opposite of the paper's sentence.

**Correction.** Cite the software: "the ntHash library (v2.4.0, SeedNtHash) warns …". Optionally add that the ntHash2 paper does not state this restriction. **C14: NARROWED** (on the citation; I did not attack the 0.747/0.818 numbers).

## F7 -- the FFT diagonal finder is MAFFT's FFT step; the paper holds MAFFT and does not cite it (C21)

Supplement S1, "FFT correlation as a diagonal finder", presents correlation peaks as an idea "the same algebra suggested". This is the FFT step of Katoh et al. 2002 (`summaries/katoh_2002_mafft.md`): correlation at every lag by FFT, then segments located from the top 20 peaks with a 30-site sliding window. MAFFT is held in the library. It is on this project's own want-list tagged for this very idea ("katoh_2002_mafft … KP-13 … the FFT step", `literature/wantlist.md`). It is absent from `manuscript/references.bib`.

The negative outcome (indels fragment diagonals) is consistent with MAFFT's design, which merges segments and then runs DP between them.

**Correction.** Cite Katoh 2002 in that bullet: "MAFFT's FFT step [Katoh 2002] used as a stand-alone diagonal finder…". **C21: NARROWED** (attribution; I did not attack the numbers).

## F8 -- C5 and C7: attacked, held

**C5.** The identity is the polarization identity. Within Das's family, Re(A C̄) = cos²θ − sin²θ = cos 2θ, so the criterion reduces to "θ ≠ π/4". This is simple, but I found no publication of it. My live search 6 (2026-11-11) found none. The project's OpenAlex queries of 2026-11-10 (`literature/eisenstein_wantlist.md`, "Queries that found nothing relevant") found none. Yu 2018 (`papers_md/yu_2018_encoding-survey.md:591-600`) treats only binary codes, and the paper already credits it. Cristea 2002 discusses "interchangeable bases … on the edges parallel to the pyrimidines and purines directions" (`papers_md/cristea_2002_genomic-signals.md:705-713`), but for the 3-D tetrahedron, where all base distances are equal, so it is not the criterion. One defect remains: the paper states this priority claim with no dated search marker. Add one. **C5: SURVIVED.**

**C7.** Each reading checks against the full-text conversions:
- Cristea eq. (4) is the θ = π/4 point and eq. (5) the negated fourth-root map (lines 810-835).
- Das 2017 gives "θ value is taken as pi/3", with the purine/pyrimidine and sum-to-zero rationale and no distances (summary, confirmed there on PDF p. 2).
- The printed 97.36% equals 1218/1251, the Eq. 8 "specificity" formula. 1218/1262 = 96.5%.
- Anastassiou makes complement a conjugation (t = a*, g = c*).
- Mendizabal's Table 6 shows spectral identity for the complementary and reverse-complementary sequences (`papers_md/mendizabal_2017_dna-numerical.md:1327-1331`).

One wording nit: "Das, Nanda and Das fixed the family's angle" implies the family existed before them. Das 2017 introduced the parametrised family and fixed θ in the same paper. "the θ = π/4 point of what was later called the trigonometric family" is accurate. I did not open the Das PDF myself. **C7: SURVIVED.**

## Steering note

The claim list does not list the Introduction's ntHash sentence ("by adding hashes of the two strands", F4.1) as a claim. It also omits the uncited (k−1)/q bound (F3) and the Supplement's FFT idea attribution (F7). All three are prior-art defects that only surface by reading the text rather than the claim list.

## Verdicts

| C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 | C9 | C10 | C11 | C12 | C13 | C14 | C15 | C16 | C17 | C18 | C19 | C20 | C21 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| UNCHALLENGED | UNCHALLENGED | UNCHALLENGED | UNCHALLENGED | SURVIVED | NARROWED | SURVIVED | UNCHALLENGED | UNCHALLENGED | UNCHALLENGED | UNCHALLENGED | UNCHALLENGED | UNCHALLENGED | NARROWED | NARROWED | REFUTED | UNCHALLENGED | REFUTED | NARROWED | UNCHALLENGED | NARROWED |
