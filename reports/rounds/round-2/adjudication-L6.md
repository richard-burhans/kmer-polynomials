# Round 2 adjudication of L6 (prior-art)

Adjudicator: skeptic for L6. I read the attack first and the claims second, as the brief orders. I read no other
round-2 attack or adjudication.

**Seal.** `sha256sum attack-L6-prior-art.md` begins `78b04ebfa257` and matches `round.json:176`. The file is the sealed
one.

**What I checked.**
- Every quotation the attack takes from the library: Mitzenmacher & Pagh 2017 :310-322, Eppstein 2011 :1470-1475 and
  :1602-1608, Boral 2014 :73-84, Shibuya 2022 :377-383, :672-678 and :740-783, Minsky 2002 :282-330, and the ntHash2
  summary.
- The km-peeler sizing code in `bench/data/kmer_poly_work/km-peeler/ibflib.c:16,121` at b8ec39b.
- `recon/kp10c.c` (decode, `flag_roots`).
- `results/kp10_vs_shibuya_2026-11-11.md` sections 1 and 3-5, and `results/kp10_vs_shibuya_2026-11-10_prior.md:98-110`.
- Round-1 history on "spurious", "resent whole" and "canonicalization".
- One web search on 2026-11-12, for Yang, Gilad & Alizadeh. It confirms "Practical Rateless Set Reconciliation", SIGCOMM
  2024, pp. 595-612, DOI 10.1145/3651890.3672219, arXiv 2402.02668. It reports about 1.35d coded symbols
  asymptotically. The library does not hold the paper. Its two `yang_2024_*` bib keys are an unrelated paper.
- I ran arithmetic only, locally: 66,477,002/16/3,238,101 = 1.2831, 1.05 × 1.222 = 1.2831,
  67 × 8 / (3,238,101/131,072) = 21.70, and 20.53 × 1.875 = 38.49. I ran no job.

**Method cards.** None governs IBLT protocols, reconciliation sizing or citation. I agree with the attacker on that.

---

## F1: a held, uncited paper extends a failed IBLT for no extra bytes (C2, C13)

**Diagnosis on the Introduction sentence: UPHELD.**
- The paper says "an invertible Bloom lookup table that fails must be resent whole" (paper :23). That is a general
  statement about IBLTs, and a held paper contradicts it in so many words. Mitzenmacher & Pagh 2017 :318-322 say:
  "double the IBLT size by adding one additional lower-order bit ... the total number of cells transmitted is the same
  as if the final IBLT had been transmitted initially".
- The paper cites neither that paper nor Yang 2024 (0 hits in `ATTACKED-paper.md`).
- The results file already says the IBLT protocol is the weak point. At :166 it lists a smarter protocol as "not
  tested".
- Round history shows "resent whole" was never raised, so this is new.

**Diagnosis on C2: NARROWED, not REFUTED.**
- The numbers reproduce: 38.49 = 20.53 × (1/8 + 1/4 + 1/2 + 1), from `iblt_protocol.tsv:6`.
- The Results sentence "the table had to be resent at double size" is true of km-peeler. km-peeler cannot reuse a
  partial listing (results :163).
- What fails is the abstract's bare "against 38.5 for the table". It reads as a property of IBLTs, and at D/8 it is
  not one.
- Under the held doubling protocol the table's total is one table at n = 1.05 D. That is about 20.5 B/diff, against the
  sketch's 20.76.
- That figure is ARITHMETIC, not a measurement. It assumes a refinement-consistent hash with power-of-two chunks, which
  km-peeler lacks (`ibflib.c:121`, `chunk_size = ceil((1.222+ε)n/r + 1)`). It also transfers the 30-of-30 success at
  1.05 D from km-peeler's own hashing.
- The sketch's lead survives at D/4 (19.40) and D/2 (16.69), against about 20.5.

**Leg not sustained.** The attack's correction says "Because an overflowing bucket is resent alone" cannot stand. That
is right only for the D/8 comparison with a table.
- As a statement of why the SKETCH's own cost stayed at 20.8 rather than about 2× its true-D cost, it is what the design
  isolates (results :156-158).
- The Discussion's "Because an overfull bucket is resent on its own, the sketch spent 11.3 ... where a one-round table
  spent 20.5" is NOT touched by F1. Doubling never sends fewer cells than the final table, so no doubling protocol
  reaches 11.3. That gain is per-bucket resending, which is exactly what the paper says. (F2 bears on that sentence.
  F1 does not.)

**The prescribed fixes.**
- **Introduction rewrite: ENDORSED, free.** "...whereas km-peeler's table, if it fails, is resent whole (an IBLT can
  instead be doubled by sending only the new cells [mitzenmacher_2017_multiparty-reconciliation], or made rateless
  [Yang et al. 2024])."
  - The mitzenmacher key is already in `references.bib:7390`.
  - Yang 2024 needs a new bib entry. That is a small acquisition, not a run.
- **Abstract option (b): ENDORSED in a stronger-scoped form, free.** "...against 38.5 for km-peeler's table resent at
  double size; a table doubled by sending only its new cells would total about 20.5 [mitzenmacher_2017], by arithmetic,
  not run".
  - The attacker's "a protocol the IBLT literature improves on" is vaguer than the number we already have.
  - The Results sentence and the Figure 4B caption need the same scope. Fixing the abstract alone leaves the Results
    and Discussion "against 38.5" (paper :146, :156) carrying the old reading. Rule 3 of the brief: fix the whole set.
- **Option (a), a measured doubling run: not required.** It is owed only if the paper wants to keep an underestimation
  advantage at D/8 as a measured claim. It costs a km-peeler code change plus one laila job.
- **The "Because" deletion: OVERTURNED except at D/8** (see above).

**Ruling.** NARROWED. "Resent whole" is false of IBLTs and must name km-peeler. The 38.5 stands as km-peeler's figure
but not as "the table's". The sketch's D/8 advantage over a held, uncited protocol is about zero by arithmetic.

## F2: the power-sum vs IBLT gap is a known 2011 result; C4 holds only for m = 64

**Citation leg: UPHELD, with a correction to the attacker's numeric parallel.**
- The quotations are verbatim. Eppstein 2011 :1470-1475: about 24 B per key for the IBF against about 10 B per element
  for CPISync. Boral 2014 :73-76: "a constant factor less communication than using IBLTs".
- The paper cites Eppstein only as a source for km-peeler (:135) and never cites Boral.
- The attacker's statement that 1.82× "is in line with Eppstein's 24/10" compares unlike units, and the agreement is
  coincidence:
  - Eppstein's IBF uses three 4-byte fields at 2× inflation.
  - Its CPISync reconciles 32-bit keys.
  - Here, km-peeler has 16-byte cells at 1.28 cells per difference, and the sketch has 8-byte field elements.
- Eppstein's comparison is also one round, with one large characteristic polynomial. The paper's own one-round sketch
  (21.70) is NOT smaller than the table (20.53). So the literature's "constant factor less" does not transfer to this
  implementation at m = 64, and the paper's measured attribution to the second round is honest for what it ran.
- The attacker's "presents it as a finding ... attributes it to its two-round design" therefore overreaches. The
  attribution is measured.

**C4 leg: NARROWED.**
- The arithmetic is verified: 3,238,101/2^17 = 24.7 per bucket and 67 × 8/24.7 = 21.70. The one-round headroom is the
  Poisson tail at small m.
- Minsky 2002 :312-315 frames the m-versus-computation trade, as the attacker says.
- "Comes from the second round" is true at m = 64 and is unmeasured at other m.

**The prescribed fixes.**
- **C4 with "with m = 64" / "at this bucket capacity": ENDORSED, free.**
- **The added clause "with larger buckets a single round costs fewer bytes and more decoding": ENDORSED only if labelled
  as an expectation with the citation.** For example: "larger buckets should lower the one-round cost at a higher
  decoding cost [minsky_2002_practical-reconciliation]".
  - The attacker's 9.7 B/diff at m = 1,024 is back-of-envelope and must not be printed as a result. Under the paper's
    own 0.8m rule the load would be about 819, not 850.
  - A measured figure is owed if wanted. It would be a whole-genome laila run, and at m = 1,024 the root flagging costs
    about 16× more per element (see F6).
- **Citing Eppstein and Boral in Results/Discussion: ENDORSED, free,** worded as: "characteristic-polynomial
  reconciliation is known to need a constant factor less communication than IBLTs [boral_2014; eppstein_2011, on
  32-bit keys]". The bare "about 10 against 24 B/diff" beside our 11.3/20.5 would invite a like-for-like reading it
  does not support. **REJECTED in that form.**
- **"Drop the Because causal clause" in the Discussion: OVERTURNED.** That clause states the measured
  two-round-versus-one-round source of the 11.3, which the 21.70 row supports at m = 64.

**Ruling.** NARROWED. The uncited prior result is real and free to add. The "comes from the second round" attribution
stands, scoped to m = 64.

## F3: "k-mer table at 1.05 times the difference" is 1.283 cells per difference (C7, C1/C6 label)

**Diagnosis: UPHELD.** All three legs check:
- The code: `ck_table[3] = 1.222` and `chunk_size = ceil((ck + ε) n / r + 1)`.
- The bytes: 66,477,002/16/3,238,101 = 1.2831 = 1.05 × 1.222, so ε = 0.
- The literature: Goodrich & Mitzenmacher's c_3 ≈ 1.222, as carried by Shibuya eq. (1).

A reader of the IBLT literature would take "1.05 times the difference" as cells per key. Nothing in C7 beyond the label
is challenged, and C6's numbers stand.

**The prescribed fix: ENDORSED, free.** In the Table 4 row and the C7 Methods sentence: "k-mer table for n = 1.05 × the
difference (1.28 cells per difference: km-peeler's default c_3 = 1.222, ε = 0)".

## F4: "61 percent against about 3 percent" compares unlike units and regimes (C8)

**Unit leg: UPHELD.**
- The source counts 1,986,807 spurious k-mers (`_prior.md:104`), and Shibuya's 3% is spurious k-mers too (:745, :783).
- The paper prints "spurious syncmers".
- Round 1 adjudication-L1 F7 endorsed "spurious k-mers numbering 61 percent". The paper carries a mis-applied round-1
  repair, which makes it a provenance defect as well as a wording one.

**Regime leg: UPHELD on the reading, but NOT new.**
- Shibuya's 3% comes from SARS-CoV-2 and simulated 30-kb sequences (:672-678). Their §3.5 conditions exactness on
  "each k-mer occurs once" (:377-383).
- The attacker says "what follows is new". It is not. Round-1 adjudication-completeness Gap 5 already found the
  61%-against-3% discrepancy real, noted that the results file explains it by regime, and endorsed reporting both
  figures (free). The attacker's history search found L1/L4 F7 but missed Gap 5.
- The attacker also says "the 7 missed differences remain untraced and are the harness question still open". Gap 5
  already attributed them to the 2 unlisted extended syncmers (507,219 of 507,221 listed every trial). Listing them
  remains optional.

**The prescribed fix: NARROWED.**
- Endorsed, free: the unit, plus "its authors report about 3 percent on viral and simulated 30-kb genomes, and note
  that the output is exact only when each k-mer occurs once".
- REJECTED: "as the unplaced scaffolds' do here".
  - It installs a causal attribution that nobody traced. The results file asserts it ("the unplaced scaffolds share
    much sequence with the chromosomes") but does not measure what fraction of the 1,986,807 occur in repeats.
  - Use instead "which may apply here; we did not trace it".
- The attacker's "the round-1 open question is largely answered by the cited paper" is the same overreach. The paper
  supplies a mechanism, not a measurement on this pair.

**Ruling.** NARROWED. The unit error is real and free to fix. The regime context is owed but was already endorsed in
round 1. The causal "as here" is not supported.

## F5: "needs no canonicalization" is a property of any strand-symmetric hash (C3)

**OVERTURNED.**
- The reading of the literature is right. ntHash2 sums the strand hashes (`kazemi_2022_nthash2` summary :13, quoting
  lines 95-97), so a power-sum sketch over ntHash2 values reduced into F_q would also need no canonical k-mer code.
- But C3's refutation condition is whether "needs no canonicalization" OVERSTATES WHAT WAS SHOWN. The attack shows
  only that the property is not exclusive, and the abstract sentence claims no exclusivity.
- It sits one sentence after the abstract's own statement that ntHash and ntHash2 "carry no proven collision bound".
  The Discussion (:156) already concedes "the hash supplies the symmetry, not that it beats canonicalization".
- The proposed "...and here with a collision bound" repeats a bound the abstract stated two sentences earlier.

**What the attacker should have checked.** Whether any surface of the paper credits the absence of canonicalization to
S rather than to strand symmetry in general. None does: round-1 completeness Gap 1 already removed the one sentence that
did. A non-exclusivity point needs a sentence claiming exclusivity.

**The fix: not owed.** It is harmless if the authors want it.

## F6: the literature puts the cost in decoding, but the paper's cost sentence counts encoding (C12)

**Wording leg: UPHELD.**
- "Recovery costs more computation ... The sketch needs about m + c = 66 field multiplications per element" gives the
  power-sum (build) count as the measure of recovery.
- The measured decode gap is real and unreported: 110.3 against 3.5 CPU-s on giga (results :174), about 31×.

**Mechanism leg: OVERTURNED.**
- The attacker imports Eppstein's and Boral's cubic decoding (Gaussian elimination in CPISync) and proposes "bounding
  each bucket at m = 64 keeps the cubic decoding ... per bucket".
- KP-10 does not decode that way. It uses Newton's identities and a Padé step (`kp10c.c:16-24`), and that per-bucket
  algebra is small: 65,536 buckets at m = 64.
- The decode time is dominated by `flag_roots` (`kp10c.c:228-237`). That step evaluates a degree-up-to-m polynomial at
  every element of every differing bucket on both sides, as results :184 states ("evaluates degree-~50 polynomials at
  every element"). So decoding costs about m multiplications per element as well, and the per-element framing is
  closer to the truth than the attack's cubic one.
- The attack also omits the end-to-end figure, which runs the other way in these implementations: 326 CPU-s for KP-10
  against about 460 for km-peeler (results :180-182). Printing "decoding was about 30 times slower" without it would
  mislead in the opposite direction.

**The prescribed fix: NARROWED.**
- Endorsed, free: "Building the sketch costs about m + c = 66 field multiplications per element, against 3 cell updates
  for the table. Decoding took about 30 times as long as the table's listing (110 against 3.5 CPU-s), because each
  side evaluates a polynomial of degree up to m at every element of a differing bucket. End to end, these
  implementations took about 326 and 460 CPU-s."
- REJECTED: the "cubic decoding ... per bucket" clause. It describes the wrong algorithm.

**Ruling.** NARROWED. The sentence counts the build phase and should say so, with the measured decode ratio. The
attack's cubic mechanism does not describe this decoder.

## Not reached

The attacker filed no finding on C5, C9, C10, C11, C14 and C15. I have nothing to rule on there.

## Correct minority

I did not read the other round-2 attacks, so I cannot say which of these points one attacker alone saw. Two are first
raised here:
- F1's held-but-uncited doubling protocol: round history shows "resent whole" and Mitzenmacher never raised.
- F3's cells-per-difference label.

Both are UPHELD on diagnosis. F1 is the most consequential prior-art finding of this lens and should be recorded even
if no other lens repeats it.

## Owed repairs (summary)

| repair | status | cost |
|---|---|---|
| Intro: "km-peeler's table ... is resent whole" + cite doubling [mitzenmacher_2017] and rateless [Yang 2024] | endorsed | free (Yang 2024 needs a bib entry) |
| Abstract, Results :146, Discussion :156, Fig 4B caption: 38.5 scoped to km-peeler; doubling ≈ 20.5 by arithmetic | endorsed | free |
| Measured doubling-IBLT run | optional, owed only to keep a D/8 advantage | km-peeler code change + one laila job |
| C4 "at m = 64"; larger-bucket expectation labelled as such [minsky_2002] | endorsed | free; a measured m = 1,024 point costs a laila run |
| Cite Boral 2014 / Eppstein 2011 for the char-poly vs IBLT constant factor (not the 24/10 numeric parallel) | endorsed | free |
| Table 4 / Methods: 1.28 cells per difference, c_3 = 1.222, ε = 0 | endorsed | free |
| C8: "spurious k-mers"; authors' regime and "exact only when each k-mer occurs once"; "may apply here; not traced" | endorsed | free |
| C12: build 66/element; decode ~30× (110 vs 3.5 CPU-s), root flagging; end to end 326 vs ~460 CPU-s | endorsed | free |

## Verdict table

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | "resent whole" is false of IBLTs per a held uncited paper (fix free); 38.5 stands for km-peeler only; the D/8 advantage vanishes by arithmetic, not measurement; the 11.3 "Because" is untouched |
| F2 | NARROWED | uncited prior result is real (free to cite); the 24/10 parallel compares unlike units; "comes from the second round" holds, scoped to m = 64 |
| F3 | UPHELD | 66,477,002/16/D = 1.2831 = 1.05 × 1.222; label fix is free |
| F4 | NARROWED | "syncmers" is a mis-applied round-1 repair (free fix); the regime point was already settled in round-1 Gap 5; the "as here" causal clause is untraced and rejected |
| F5 | OVERTURNED | the sentence is true and claims no exclusivity; the abstract and Discussion already scope it; no fix owed |
| F6 | NARROWED | 66 is a build count and the decode gap is unreported (free fix); the cubic-decoding mechanism is wrong for kp10c, whose decode is per-element root flagging, and end to end KP-10 was faster |
