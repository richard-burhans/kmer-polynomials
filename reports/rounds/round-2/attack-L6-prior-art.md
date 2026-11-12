# Round 2 attack, lens L6: prior-art

Attacker: prior-art lens. Frozen text: `reports/rounds/round-2/ATTACKED-paper.md`, `ATTACKED-ledger.md`.

**Where the literature is.** The brief declared no appendix. `refs/catalog_index.json`, which the lens paragraph names,
does not exist in this repo or in `../library`. I used `../library/references.bib`, the full texts in
`../library/papers_md/<slug>.md` and the summaries in `../library/papers_md/summaries/`. The method cards are in
`../library/methods/`. None of them governs IBLT sizing, set-reconciliation protocols or citation of prior work, so no
finding below rests on a card. The statistical claim (C5) belongs to the statistics lens.

**Prior rounds.** I ran `scripts/round history` on "resent whole", "Mitzenmacher", "CPISync", "Eppstein", "rateless",
"1.05", "3 percent", "syncmer", "km-peeler", "Shibuya", "PinSketch", "minisketch" and "canonicalization".
- **Never raised:** "resent whole", Mitzenmacher, CPISync, Eppstein and rateless have no prior mention in any round.
- **Raised before:** round 1's adjudication-L1 F7 and adjudication-L4 F7 settled the 38%/61% arithmetic, and round 1's
  adjudication-L6 F2 settled the PinSketch, Minsky 2002 and Erlay attribution.
- **How this attack relates:** F4 below builds on round 1's F7. What it adds is new: the unit word "syncmers" and the
  regime of the authors' 3%. I do not re-file the attribution point.

**Web search.** One search on 2026-11-12, for "Rateless Invertible Bloom Lookup Tables set reconciliation Yang Gilad
Alizadeh SIGCOMM 2024". It found Yang, Gilad & Alizadeh, "Practical Rateless Set Reconciliation", SIGCOMM 2024,
DOI 10.1145/3651890.3672219, arXiv 2402.02668. The library does not hold this paper.

**Steering.** The claim list confines C2's straw-man test to the resend protocol "km-peeler supports", and C13's to
"IBLTs in the literature (or km-peeler)". The C2 wording steers the attack toward the implementation and away from the
literature, which is where the protocol the paper should have compared against lives (F1).

## F1 -- a held, uncited paper extends a failed IBLT for no extra bytes, which removes C2's 38.5 and C13's "resent whole"

**Claims:** C2 and C13, plus the Introduction's "an invertible Bloom lookup table that fails must be resent whole".

**Evidence.**
- **The held paper.** Mitzenmacher & Pagh 2017 (`../library/papers_md/mitzenmacher_2017_multiparty-reconciliation.md:310-322`)
  is in `references.bib` but cited nowhere in the paper (grep of `ATTACKED-paper.md`: 0 hits). Without an upper bound on
  d, they say, the prior work [17, 28] "correspond[s] to repeated doubling", and:

  > "Another option is to double the IBLT size by adding one additional lower-order bit to each hash value. Then it
  > suffices to send every odd-numbered cell of the IBLT arrays, since the even-numbered cells can be found by
  > subtraction from the old IBLT. In this way, the total number of cells transmitted is the same as if the final IBLT
  > had been transmitted initially."

- **The premise holds here.** IBLT cells are linear (counts add, keysums XOR), so a size-m cell is the sum of its two
  size-2m children. The difference protocol `T_A − T_B` stays linear under that refinement. The receiver keeps its own
  old table.
- **At D/8 the protocol's total is one table at 1.05·D.** It starts at n = 1.05·D/8 and doubles three times to
  n = 1.05·D, which succeeded 30 of 30 (`results/kp10_vs_shibuya_2026-11-11.md:108`).
  - Total ≈ 20.53 B/diff in 4 rounds, against the sketch's 20.76 B/diff in 5 rounds (`results/...-11-11.md:142`).
  - The paper's 38.5 is 20.53 × (1/8 + 1/4 + 1/2 + 1) = 38.49: every failed table resent in full.
- **The results file concedes the gap.** At `:166` it lists "a smarter IBLT protocol ... not what km-peeler
  implements" as **not tested**. The paper drops that caveat.
- **The paper's own protocol is published, uncited and resends nothing whole.** Rateless IBLT (Yang, Gilad & Alizadeh,
  SIGCOMM 2024; web, 2026-11-12) streams extra coded symbols until decoding succeeds and needs no estimate of the
  difference.

**Consequence.**
- **C2's comparison changes.** "20.8, against 38.5 for the table" compares the sketch with a naive protocol. With the
  doubling protocol of a paper the project holds, the expected figures at D/8 are about 20.8 against about 20.5. The
  sketch keeps its lead only at D/4 (19.4) and D/2 (16.7).
- **C13's general sentence is false as written.** "An invertible Bloom lookup table that fails must be resent whole" is
  true of km-peeler only.

**Correction.**
- **Introduction:** "...whereas km-peeler's table, if it fails, is resent whole (an IBLT can instead be doubled by
  sending only the new cells [mitzenmacher_2017_multiparty-reconciliation], or made rateless [Yang 2024])."
- **Abstract C2:**
  - Either compare against the doubling protocol (owed: one laila job; km-peeler's `ibflib.c:121` partition sizing
    needs power-of-two chunks, so this is a small code change), or
  - say "against 38.5 for km-peeler's table resent at double size, a protocol the IBLT literature improves on".
- **"Because an overflowing bucket is resent alone"** cannot stand as the reason the sketch beats IBLTs. With the
  doubling protocol, an IBLT also sends no wasted cells at D/8.

**Verdict: C2 REFUTED, C13 REFUTED.**

## F2 -- the power-sum vs IBLT size gap is a known 2011 result, and C4's "comes from the second round" holds only for m = 64

**Claims:** C1, C4 and the Discussion's "Because an overfull bucket is resent on its own, the sketch spent 11.3 ... where
a one-round table spent 20.5".

**Evidence.**
- **Eppstein et al. 2011 (cited for another purpose).** `eppstein_2011_set-difference.md:1470-1475`: the IBF grows at
  "24 Bytes ... per key in the difference", while "CPISync requires less memory in practice at approximately 10 bytes
  per element in the difference". Lines 1606-1608 add that CPISync needs "only 10d bytes of transmission". That is one
  round, single bucket, with no second round.
- **Boral et al. 2014 (held, uncited).** `boral_2014_multiparty-reconciliation.md:73-76`: "using characteristic
  polynomials generally requires a constant factor less communication than using IBLTs, and is almost optimal".
- **So the sketch's size advantage over an IBLT is the expected, published result.** The paper's ratio of 1.82×
  (11.3/20.5) is in line with Eppstein's 24/10. The paper presents it as a finding, attributes it to its two-round
  design, and cites neither comparison.
- **The one-round 21.7 is a property of m = 64 with 2^17 buckets, not of power sums.**
  - Arithmetic: 3,238,101 / 131,072 = 24.7 mean differences per bucket against capacity 64, and 67 × 8 / 24.7 = 21.7.
    The headroom is the Poisson tail over many small buckets.
  - Minsky, Trachtenberg & Zippel 2002 (`minsky_2002_practical-reconciliation.md:282-330`, cited) frame exactly this
    as the trade between partition bound m and rounds or computation.
  - Back-of-envelope (not run; owed if the authors want a number): with m = 1,024, a mean load near 850 sits about
    5 SD below capacity for ~3,800 buckets. That gives about 1,027 × 8 / 850 ≈ 9.7 B/diff in one round, at higher
    decode cost.

**Correction.**
- **C4:** "Sized for a single round with m = 64, the sketch needed 21.7, so at this bucket capacity its smaller size
  comes from the second round; with larger buckets a single round costs fewer bytes and more decoding
  [eppstein_2011_set-difference; minsky_2002_practical-reconciliation]."
- **C1 and the Results:** add that a gap of this order between characteristic-polynomial and IBLT reconciliation was
  measured by Eppstein et al. (about 10 against 24 B/diff) and stated by Boral et al.
- **Discussion:** drop the "Because" causal clause (see F1).

**Verdict: C4 NARROWED, C1 NARROWED** (the numbers stand; the sentence's implied novelty and causal reading do not).

## F3 -- "k-mer table at 1.05 times the difference" is 1.283 cells per difference

**Claims:** C7, and the Table 4 row label shared by C1 and C6.

**Evidence.**
- **km-peeler sizes per its paper's eq. (1).** Shibuya et al. (`shibuya_2022_genomic-reconciliation.md:212-219`) and
  Goodrich & Mitzenmacher (`goodrich_2011_iblt.md:482`) give m = n(c_r + ε) with c_3 ≈ 1.222.
- **km-peeler's code does the same:** `bench/data/kmer_poly_work/km-peeler/ibflib.c:16` (`ck_table ... 1.222`) and
  `:121`, `chunk_size = ceil((ck_table[r] + epsilon) * n / r + 1)`, at commit b8ec39b.
- **The bytes confirm it.** 66,477,002 B / 16 B per cell / 3,238,101 = **1.2831** cells per difference, which equals
  1.05 × 1.222.
- **Why the label matters.** A reader who knows the held IBLT literature reads "table at 1.05 times the difference" as
  1.05 cells per key, below the peeling threshold. A 30-of-30 success there would be impossible. The label also hides
  where the table's bytes go: 1.283 cells × 16 B.

**Correction.** Table 4 row: "km-peeler, k-mer table for n = 1.05 × the difference (1.28 cells per difference, its
default c_3 = 1.222, ε = 0)". Use the same wording in the Methods sentence of C7.

**Verdict: C7 NARROWED. C6 SURVIVED** on its numbers, with the same label fix.

## F4 -- the extended-syncmer "61 percent against about 3 percent" compares unlike units and unlike regimes

**Claim:** C8.

**Prior round.** Round 1 adjudication-L1 F7 and L4 F7 settled the arithmetic (61% of D, 38% of the returned set) and
endorsed "spurious k-mers numbering 61 percent of the true difference". What follows is new.

**Evidence.**
1. **The unit is wrong.** The paper says "spurious **syncmers** numbering 61 percent". The source says 1,986,807
   spurious **k-mers** (`results/kp10_vs_shibuya_2026-11-10_prior.md:104-105`; ledger row 111). Shibuya's "about 3%"
   is also spurious k-mers (`shibuya_2022:744-747, 782-783`). As printed, the sentence sets a syncmer count against a
   k-mer difference, which is not what was measured.
2. **The authors' 3% comes from another regime, and their own text predicts the cannabis result.**
   - The 3% (`shibuya_2022:672-678`, Table 1) is from SARS-CoV-2 genomes (mean diff 325) and 30-kb random sequences
     with point mutations at p_m = 0.001 and 0.01.
   - Their Section 3.5 (`:377-383`) says the output is a superset "due to the fact that an extended syncmer can contain
     a k-mer which belongs to another extended syncmer common to both datasets". It is exact only "where each k-mer
     occurs once".
   - The GigaDB difference is dropped unplaced scaffolds that "share much sequence with the chromosomes"
     (`..._prior.md:106-107`), in a repetitive ~800-Mb plant genome. That is exactly the case the authors flag.
   - So the high spurious rate is the method's documented behaviour on repeats, not evidence of a harness error. The
     round-1 open question is largely answered by the cited paper. The 7 missed differences remain untraced and are the
     harness question that is still open.
   - The bare "against about 3 percent in the authors' tests" invites the reader to suspect either the harness or the
     method, and omits the authors' own condition.

**Correction.** "...returned spurious k-mers numbering 61 percent of the true difference. Its authors report about 3
percent on viral and simulated 30-kb genomes, and note that k-mers shared with common extended syncmers are returned
when k-mers repeat, as the unplaced scaffolds' do here [shibuya_2022_genomic-reconciliation]. It also missed 7 differences,
which we did not trace."

**Verdict: C8 NARROWED.**

## F5 -- "needs no canonicalization" is a property of any strand-symmetric hash, which the held literature already supplies

**Claim:** C3 (abstract claim 6).

**Evidence.**
- **ntHash already combines its strand hashes.** ntHash takes the smaller of its two strand hashes and ntHash2 adds them
  (`kazemi_2022_nthash2` summary, per round-1 adjudication-L4 F11; paper Introduction). Either is orientation-invariant
  without canonicalizing the k-mer.
- **The same power-sum sketch over ntHash2 values reduced into F_q** would also "need no canonicalization". It would
  recover hash values that each side maps back by rescanning, exactly as KP-10 does (`results/...-11-11.md:65-70`).
- **What S adds is the collision bound, not the absence of canonicalization.** The paper's Discussion concedes the
  symmetry point ("this shows that the hash supplies the symmetry, not that it beats canonicalization"). The abstract
  sentence instead attributes the property to "these hashes".
- **The other parts hold.** The as-written k = 31 and k = 63 results and "linear" (Dodis 2008, Minsky 2003; settled in
  round 1) are not contested here.

**Correction.** Abstract: "Power sums of these hashes form a linear set sketch that, like any strand-symmetric hash,
needs no canonicalization, and here with a collision bound."

**Verdict: C3 NARROWED.**

## F6 -- the literature locates the cost in decoding, but the paper's cost sentence counts encoding

**Claim:** C12.

**Evidence.**
- **The literature puts the cost in decoding.**
  - Eppstein 2011 (`:1602-1608`): "the cubic decoding complexity for CPISync dramatically increases its latency", with
    "d^3 memory operations to solve its matrix" against "3d × hash_count" for IBF.
  - Boral 2014 (`:76-84`): about O(d^3) to recover, and O(|S_i| d) to encode.
- **The paper's sentence counts the wrong phase.** "Recovery costs more computation than a table listing. The sketch
  needs about m + c = 66 field multiplications per element" gives a per-element **insertion** count as the measure of
  **recovery**.
- **The measured recovery gap is not reported.** The results file has KP-10's decode "about 30× slower than the IBLT's
  diff and list" (`results/...-11-11.md:19-20`). The paper leaves it out.
- **Bounded buckets (m = 64) are what keep decoding tractable,** which is Minsky 2002's point. That should be said.

**Correction.** "Building the sketch costs about m + c = 66 field multiplications per element, against 3 cell updates for
the table, and decoding was about 30 times slower than the table's listing; bounding each bucket at m = 64 keeps the
cubic decoding of characteristic-polynomial methods [eppstein_2011_set-difference] per bucket."

**Verdict: C12 NARROWED.**

## Not reached through this lens

C5, C9, C10, C11, C14 and C15. I found nothing concrete in the literature that changes them.
- **C14:** "a canonical packed code and a universal hash give strand symmetry with a proof" is standard, and I found no
  held paper that contradicts it.
- **C15:** the ledger rows I opened (105, 109-112, 149-155, 199) say what their sources say. The F3 label issue is in the
  paper, not the ledger.

| claim | verdict |
|---|---|
| C1 | NARROWED |
| C2 | REFUTED |
| C3 | NARROWED |
| C4 | NARROWED |
| C5 | UNCHALLENGED |
| C6 | SURVIVED |
| C7 | NARROWED |
| C8 | NARROWED |
| C9 | UNCHALLENGED |
| C10 | UNCHALLENGED |
| C11 | UNCHALLENGED |
| C12 | NARROWED |
| C13 | REFUTED |
| C14 | UNCHALLENGED |
| C15 | UNCHALLENGED |
