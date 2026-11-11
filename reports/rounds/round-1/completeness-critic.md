# Round 1 completeness critic: kmer-polynomials manuscript

Target: `reports/rounds/round-1/ATTACKED-paper.md` (draft d049b77).

**What I read before asserting any absence.** I read all of these in full:
- `CRITIC-BRIEF.md` and `round.json`;
- the seven attacks, `attack-L1-correctness.md` through `attack-L7-reproduction.md`;
- the seven adjudications, `adjudication-L1.md` through `adjudication-L7.md`.

Every "nobody examined X" below was also checked with a keyword grep over those 14 files. I searched for other names for each gap too: "canonical code", "kmerset", "mix64", "universal", "Carter", "min-wise", "FracMinHash", "sourmash", "O3", "march", "vector instr", "error analysis", "in seconds", "Voss", "white noise", "about 3", "unplaced", "2.5" and "colliding pairs". I quote the grep result next to each gap.

`history` withholds round 1, and this is the first round on the repository, so there are no closed rounds to search.

## 1. What filed, per lens

All seven declared lenses filed one attack, and each attack was adjudicated. `aborted`, `unbriefed` and `lens_drops` are empty, so **no lens is UNRUN**.

| lens | findings | claims given a verdict other than UNCHALLENGED |
|---|---|---|
| L1 correctness | 12 | 17 |
| L2 sampling | 7 | 5 |
| L3 statistics | 9 | 12 |
| L4 provenance | 15 | 21 |
| L5 inference | 11 | 17 |
| L6 prior-art | 8 | 9 |
| L7 reproduction | 11 | 18 |

Lenses that engaged each claim:

| claim | lenses | claim | lenses | claim | lenses |
|---|---|---|---|---|---|
| C1 | 5 | C8 | 4 | C15 | 7 |
| C2 | 4 | C9 | 6 | C16 | 6 |
| C3 | 4 | C10 | 6 | C17 | 5 |
| C4 | 5 | C11 | 6 | C18 | 5 |
| C5 | 5 (priority: L6 only) | C12 | 5 | C19 | 3 |
| C6 | 3 (priority: L6 only) | C13 | 6 | C20 | 3 |
| C7 | **2** (L4, L6) | C14 | 6 | C21 | 3 |

**Where effort went.**
- **Concentrated.** C16, the claim that overflow detection is "deterministic", drew a finding from six lenses. All six make the same point. The 2.7×10⁻¹⁷ arithmetic in C1 drew five findings, and the 1/20 break-even (C11) drew six.
- **Thin.** The two priority claims were each examined by one lens.
  - C5: L6 ran a single WebSearch, and adjudication-L6 F8 calls it "thin".
  - C6: L6 ran six searches.
  - L4 F15 only noted that the paper carries no dated search marker for either claim.
- **Two lenses only.** C7 was read against the PDFs by L4 alone. L6 says it "did not open the Das PDF myself".

**UNCHALLENGED is not the same as UNRUN, condition by condition.** Each claim's own "Refuted if" list names its conditions. Four of those conditions were run by nobody:
- **C10, "if the comparison is unfair (flags …)".** No attack checks the compile flags (grep "O3", "march": 0 hits in all 14 files). I closed it myself, and it is **not a defect**. `hash/build_nthash.sh` builds ntHash c26bd45 from source with `OPT="-O3 -march=native"  # must match hash/Makefile OPT`. `hash/Makefile:4` uses the same flags.
- **C10, the scalar code matching AVX-512 "without any vector instructions".** No attack checked this (grep "vector instr", "auto-vec": 0 hits). I spot-checked the local build only. In `objdump -d hash/kp7e.o`, the ymm instructions in `kp7e_sym_seq` are table-initialisation stores, not the rolling loop. This is **low risk**. The laila build was not inspected.
- **C21, "if a 'did not pay off' verdict does not follow".** L4, L6 and L7 checked numbers and attribution. L5, the inference lens, marked C21 UNCHALLENGED. See gap 8.
- **C19, the journal version of Marçais.** L4 noted that only the preprint was read. Nobody read the journal version. This is a minor gap.

## 2. Gaps that matter

### Gap 1 (framing, strongest). In the sketch experiment, orientation invariance comes from canonicalising before hashing, not from the strand-symmetric hash

The paper attributes the draft-pair result to S several times:
- Introduction: "that is where strand symmetry pays off … two assemblies of the same k-mers give the same sketch whatever the orientation of their contigs";
- Results: "Our sketch recovered the 7,931 exactly … **because** a strand-symmetric hash makes the sketch blind to orientation";
- Abstract and Discussion: the orientation result is presented as what the hash "gives".

The code that produced those numbers does something else:
- `recon/kp10.py:14-17` says: "Inputs are the partitioned canonical k-mer sets of recon/kmerset.py … The bucket key is mix64(code)".
- `recon/kmerset.py:11` says: "Canonical = min(forward code, reverse-complement code)".
- `kp10.py:50-51` builds each side as `SideE.build(codes, sym_hash(codes), mix64(codes))`. Here `codes` are already canonical.

So every element is canonicalised by packed code before S is evaluated, and the bucket key is a heuristic mix of the canonical code. Any hash of these inputs, strand-symmetric or not, would give the same orientation invariance. The measured 137 KB result therefore cannot support "because a strand-symmetric hash". The Discussion concedes that a canonical k-mer table is orientation-blind, but it never says that its own sketch *was* built on a canonical k-mer table.

A second consequence follows. Bucketing by the packed code restricts the sketch, as built, to k ≤ 32. That is exactly the regime where canonicalisation is free.

**Nobody examined this.** Grep for "kmerset", "mix64", "canonical code" and "canonicaliz" finds:
- L7's orientation test, which uses mix64 for its own sampling;
- L4, L5 and adjudication-L7 on canonicalising km-peeler's syncmer strings;
- nothing about how the KP-10 inputs were built.

Adjudication-L4 argued from "the canonical k-mer difference is only 7,931" to the orientation cause in km-peeler's run. It did not turn that observation on KP-10's own inputs.

**Why it matters.** The abstract's orientation sentence and the Introduction's motivating paragraph assert a causal role for S that the experiment did not test.

**What would close it.** Either:
- state that the inputs were canonicalised and drop "because"; or
- rerun KP-10 on forward-strand, as-written k-mers with a strand-symmetric bucket key (for example, bits of S), and show the draft pair still gives 7,931.

### Gap 2 (framing). No lens asked whether a simpler provable strand-symmetric hash already exists for k ≤ 32, which covers every experiment in the paper

The Introduction asks "what can be proven instead" of ntHash's measured quality, and the abstract builds its contrast on that question. For k ≤ 32 there is a simple alternative: take the canonical 2-bit code min(f, r), which the paper's own sampler already rolls ("exact for k ≤ 32"), and apply any seeded universal or polynomial hash. That hash is strand-symmetric by construction, carries a proven pairwise bound of order 1/q, and costs about one multiplication on top of the packed roll.

Every experiment in the paper uses k = 31 (or 7 and 8). So the paper never measures the regime where S is the only provable option, k > 32. It also never states S's advantage over this baseline: field-linear algebra, and spaced seeds.

**Nobody examined this.** Grep for "universal", "Carter", "Wegman", "min-wise", "k > 32" and "canonical code" finds no hit on this point. L6's searches looked for published strand-symmetric *polynomial* hashes only.

**Why it matters.** It changes what the "Second contribution" and the abstract's first two sentences can claim. With no baseline, the paper's speed comparison (24.7 against 11.4 cycles) also never asks what a provable hash costs at k ≤ 32.

**What would close it.** Add a row to the speed table for canonical code plus a seeded universal hash. Then scope the provability argument to k > 32, to spaced seeds and to the algebra.

### Gap 3. The sampler's speed advantage comes from the heuristic packed-code filter, and its comparison baseline was not matched

The sampler decides which k-mers to keep with mix(f) + mix(r) on 2-bit packed codes. The paper says this filter is heuristic. The sampler computes S only for the kept k-mers. Batched ntHash is timed hashing every k-mer.

The sampler's intercept is 9.96 cycles (L2, L3, L5 and L7 all recompute it), below ntHash's 11.37. So the advantage at sparse rates is the cheap filter. The provable hash only adds cost, 28.5 cycles per kept k-mer.

An ntHash user could apply the same packed-code filter and compute ntHash only for kept k-mers. Nobody measured that baseline. There is also a conceptual gap. In FracMinHash, the selection hash is the one whose uniformity makes containment estimates unbiased. Here, selection is by the heuristic mix. "Every kept k-mer carries the provable hash" is true, but the proof does not cover the sampling decision. C13 *measured* uniformity, but it did not prove it. The paper's contrast with ntHash, "measured, not proven", therefore applies to this pipeline too.

**Nobody examined this.** All six C11 findings are about the precision of the 1/20 crossing and the inline-versus-collected modes. Grep for "heuristic" finds it only in the C16 bound disputes.

**What would close it.**
- Time "packed filter + ntHash on kept" in the same job.
- State that the selection is heuristic, so the provable guarantee covers collisions among kept k-mers, not which k-mers are kept.

### Gap 4. "FracMinHash sketches commonly sample far more sparsely than that" has no citation and no ledger row

This is the sentence that makes the 1/20 result matter. Grep finds no Irber, sourmash, FracMinHash or "scaled" entry in `manuscript/references.bib`. `manuscript/ledger.md` has no row for it. No attack mentions it (grep "FracMinHash", "sourmash": 0 hits).

It is probably true, since sourmash's default is scaled = 1000. It is still an unsourced empirical claim about practice.

**What would close it.** Cite the FracMinHash and sourmash papers with their default scaled values.

### Gap 5. The competitor's extended-syncmer row was never examined for harness error: 61% of D spurious against the authors' about 3%, and 7 differences missed by a superset method

`results/kp10_vs_shibuya_2026-11-11.md:105` says "1,986,807 spurious k-mers (61% of D; **the paper saw about 3%**)" and "**missed 7**" true differences in all 5 trials. The Shibuya library summary (lines 27-29) gives "within about 3% of the true set difference size" and describes the output as "a superset". The run therefore reports a spurious rate about 20 times the published one, and it misses elements that the method's design should return.

The results file offers one explanation for the spurious excess: unplaced scaffolds share sequence with the chromosomes. Nothing explains the 7 misses. Both could come from a harness defect, for example N or soft-mask handling, or canonicalisation during scoring in `recon/kmpeeler_run.py:169-174`, which canonicalises the listed syncmers' k-mers for scoring. A defect there would inflate the competitor's error, and the paper uses that error to dismiss the competitor ("2.7 times smaller … but returned a superset").

The paper does not tell the reader about the 3% against 61% discrepancy.

**Nobody examined this.** L1, L4 and L7 checked only whether 61% should be 38%. Adjudication-L1 mentions "Shibuya's 'about 3%'" only as a reason to keep the per-D basis. Grep "missed 7": L4 F7 and L7 quote it, and neither asks why.

**What would close it.**
- List the 7 missed k-mers and find their cause.
- Report the published 3% next to the observed 61%, with the regime difference, in the paper.

### Gap 6 (a claim nobody could attack as worded). "Reconciliation sketches can run on the hash values with one error analysis"

The Discussion sells the field structure as "one error analysis" covering collisions and sketches alike. The round's own adjudications have already ruled that no such analysis exists:
- adjudication-L5 F6: the per-pair bound 2(k−1)/q "is not the governing quantity";
- adjudication-L6 F1: a proven bound would need "a union over candidate subsets … not obviously small";
- adjudications L3 and L7 overturned q^{-c} as unproven.

Even so, no finding names the "one error analysis" sentence, and no correction touches it (grep "error analysis": 0 hits). The sentence is vague enough that no attack keyed to it. It now conflicts with what the round established.

**What would close it.** Delete it, or replace it with "the collision bound covers the hash; the sketch's overflow check has no derived bound".

### Gap 7. Unmeasured predictions stated as facts

- "A C implementation would hash a genome in seconds, but the power sums would still dominate" (Results, sketch cost). Nothing was measured, and no lens examined it (grep "in seconds": 0 hits).
- "So it pays only for hash-only graphs with k > 32" (Supplement, de Bruijn). This is not measured for k > 32 beyond the laila table at 47 and 63, which shows a slowdown there too.

**What would close it.** Mark both as estimates with their arithmetic, or measure them.

### Gap 8. The inference leg of C21 is UNRUN, and one Supplement citation was never read

- No lens tested whether the Supplement's "did not pay off" verdicts follow from the data: "faster only for arms longer than 32 bases", "a reconciliation sketch is as large as the set", "pays only for k > 32".
- Nobody checked the reading behind "Windows of k-mer length behave like white noise [@voss_1992_fractal-dna]". The only mention of Voss in the 14 files is adjudication-L6, which notes that the bullet cites prior art. The library summary (`voss_1992_fractal-dna.md:14,23-24`) reports 1/f^β long-range correlations, a white-noise region at high frequency, and a period-3 peak inside it. The paper's claim may hold at 31-base scale, but the cited paper's headline finding is the opposite of white noise. For an eight-coefficient bound in a 31-base window, the period-3 peak falls near coefficient 10.

This gap carries low stakes, because these are negative results in a supplement.

**What would close it.** Have one inference pass over S1, and quote Voss's white-noise-region sentence.

## 3. Framing of the round

- **The runner wrote the claim list and is the paper's author** (`round.json`: `runner_is_author: true`, `claim_list_author: runner`). Twenty of the 21 claims are "unclassified".
- **Attackers already flagged items missing from the list.** Five attackers wrote steering notes about items the list left out. Those items surfaced only because individual lenses read the text:
  - L1: x and λ derived from one seed;
  - L4: the spaced-seed timing and the CM31 bound;
  - L5: the two-set decoder;
  - L6: the Introduction's ntHash sentence, the uncited (k−1)/q, and MAFFT.
- **No lens was assigned to baseline adequacy**, meaning "is the comparison the right one?". Gaps 1 to 3 all sit there. The declared lenses check whether each number is right. None asks whether S is what produced the advantage being claimed.
- **Gap I considered and did not file.** I checked whether the per-genome collision allowance (about 2.5 pairs, from L3 F2 and L5 F4) breaks the sketch's "exact" claim. It does not. Buckets are keyed independently of S, so only same-bucket pairs matter. That gives about 3.2×10⁶ × 6,700 × 2.6×10⁻¹⁷ ≈ 6×10⁻⁷ expected per run.

## Gap list

| # | gap | why it matters | what would close it |
|---|---|---|---|
| 1 | KP-10 inputs canonicalised by min(code) before S; bucket key mix64(canonical code) (`recon/kp10.py:14-17`, `kmerset.py:11`) | "because a strand-symmetric hash" (Results, Abstract, Introduction) is untested; any hash would give the 7,931 | say inputs were canonical, or rerun on as-written k-mers with an S-derived bucket key |
| 2 | no canonical-code + universal-hash baseline at k ≤ 32 (all experiments are k = 31) | "what can be proven instead" has a simple answer at k ≤ 32; S's regime (k > 32) is never measured | add the baseline row; scope provability to k > 32, spaced seeds and algebra |
| 3 | sampler speed comes from the heuristic packed filter (intercept 9.96 < 11.37); no "filter + ntHash on kept" baseline; selection is not covered by the proof | the 1/20 claim compares unmatched pipelines; "measured, not proven" applies to this sampler too | time the matched baseline in the same job; state that selection is heuristic |
| 4 | "FracMinHash commonly samples far more sparsely": uncited, no ledger row | it is the sentence that gives C11 its relevance | cite FracMinHash/sourmash defaults |
| 5 | extended-syncmer run: 61% of D spurious against Shibuya's about 3%, and 7 misses by a superset method, unexplained | a harness defect would inflate the competitor's error that the paper uses against it | inspect the 7 misses; report 3% against 61% |
| 6 | "one error analysis" (Discussion) | contradicted by the round's own adjudications (no overflow bound exists); no finding named it | delete or reword |
| 7 | "C implementation … in seconds"; "pays only for k > 32" | unmeasured predictions stated as fact | label as estimates or measure |
| 8 | C21 inference leg UNRUN; Voss "white noise" reading unchecked | low stakes; the cited paper's headline is 1/f | one inference pass over S1; quote Voss |
| — | UNRUN but closed by me: C10 flags (same `-O3 -march=native`); "no vector instructions" (local spot check clean) | not defects | none |
| — | thin coverage: C5 priority (one WebSearch), C6 priority (one lens), C7 (two lenses) | priority claims rest on one lens | a second, independent dated search for C5 and C6 |
