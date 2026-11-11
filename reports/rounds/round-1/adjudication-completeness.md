# Round 1 adjudication of the completeness critic

Target: `reports/rounds/round-1/completeness-critic.md`, read before the claims as the brief requires. The filing is unsealed, so I quote the text I rule on in each section.

`scripts/round history "orientation" --kind adjudication` reports "no prior mention ... in 1 round directories", and round 1 is withheld because it is still open. No earlier ruling bears on these gaps. No method card in `refs/methods/` governs them, because they concern coverage and framing, not a statistical method.

Reproduction script: `/home/agent/.claude/jobs/0d0eac3e/tmp/ADJ-crit/g1.py`.

## Section 1: conditions the critic says nobody ran, and closed itself

- **C10 compile flags.** The critic wrote "`hash/build_nthash.sh` ... `OPT="-O3 -march=native"` ... `hash/Makefile:4` uses the same flags". I confirmed both files, so this is closed and not a defect.
- **C10 "without any vector instructions".** The critic checked only the local build and said so ("The laila build was not inspected"). I accept that as stated: low risk, and only partly closed.
- **C19 journal version.** A minor gap and correctly scoped.

These are coverage facts, not findings against a claim, so they need no table rows.

## Gap 1: KP-10's orientation invariance comes from canonicalization. UPHELD on the diagnosis; NARROWED on "untested"

The critic quoted: "`recon/kp10.py:14-17` ... 'Inputs are the partitioned canonical k-mer sets' ... `kmerset.py:11` 'Canonical = min(forward code, reverse-complement code)' ... `SideE.build(codes, sym_hash(codes), mix64(codes))`."

**Premise verified in the code.**
- `kmerset.build` calls `canonical_kmers(encode(seq), k)`. Its docstring in `prototype/kmers.py:31` reads: "min of forward and reverse-complement 2-bit codes".
- `kp10.side` passes those canonical codes to both `sym_hash` and `mix64`.
- `kp10.run` scores against `np.setdiff1d` of the same canonical sets, so the 7,931 target is itself a canonical-set difference.

**Re-run (g1.py, k=31, a random 20 kb sequence and its reverse complement):**
- canonical k-mer sets are equal across orientation: `True`;
- **the forward-only polynomial F (not strand-symmetric) on canonical inputs gives equal multisets: `True`**, so any element hash would have produced the orientation invariance;
- S(as-written) equals S(canonical) element by element: `True`;
- S value sets on the as-written strands are equal across orientation: `True`;
- **mix64 keys on as-written codes are equal across orientation: `False`.** The bucket key is the component that actually needs canonical input.

**What survives.** The sentence at paper line 130, "because a strand-symmetric hash makes the sketch blind to orientation", gives S a causal role that this experiment cannot isolate, since canonicalization had already made the inputs orientation-free. The results file states the same attribution ("KP-7's symmetric hash makes KP-10 independent of ..."). Introduction line 23 ("When the hash is strand-symmetric, so is the sketch, and ... whatever the orientation") is true as mathematics, but it is illustrated by an experiment that does not exercise it. Methods line 59 ("by a mix of its packed code") does not say which code is used, and nowhere does the paper say that the sketch's own inputs were canonical k-mers. That makes the Discussion's concession that "a table of canonical k-mers is orientation-blind too" incomplete: the paper's own sketch is built on such a table.

**Narrowing.** "Untested" applies to the causal "because" inside the experiment, not to S's strand symmetry. That symmetry is tested exhaustively at k=7/8 and by `kp7e_hash.self_test` (`sym_hash(codes) == sym_hash(rc_codes(codes))`), and I re-confirmed it above. The orientation claim is therefore a correct mathematical consequence and an unsupported experimental attribution.

**No false gap.** I grepped all 14 round files for canonical, kmerset, kp10.py and orientation. L4 F-C17, L5 F8, L7 F10 and adjudication-L4 all examine the cause of km-peeler's 48 M, and none examines KP-10's inputs. The critic is the only agent that saw this. I record it as a correct minority (rule 8).

**Critic's own premise.** "Bucketing by the packed code restricts the sketch, as built, to k ≤ 32." UPHELD: the codes are `uint64` holding 2 bits per base.

**Repairs.**
- (a) "State that inputs were canonicalised and drop 'because'": **ENDORSED; the fix is free.** Suggested Results wording: "Both sides were canonical k-mer sets, so the sketch was orientation-blind by construction; a strand-symmetric hash makes canonicalization unnecessary, which this experiment did not test." In the Discussion, say the comparison is against the extended-syncmer method only, "as for any canonical k-mer table, including ours". Apply the same change to the results file's "Cause" line.
- (b) "Rerun on as-written k-mers with an S-derived bucket key": **endorsed as the only test of the causal claim, but not free.** It needs a new KP-10 run, about an hour per pair in numpy. It also needs care: buckets keyed on bits of S make the bucket depend on the element value, so the decode-failure accounting must be rechecked. Without the rerun, (a) must be applied.

## Gap 2: no canonical-code plus universal-hash baseline at k ≤ 32. NARROWED

**The diagnosis holds.** For k ≤ 32, min(f, r) is an injective 64-bit representative of the double-stranded class. Any seeded universal or polynomial hash of it into F_q is strand-symmetric and carries a provable ~1/q pairwise bound. Gap 1 shows that this is exactly the pipeline KP-10 ran, with S as the element hash. The grep for universal, Wegman, k > 32 and min(f, r) found no attack that raised this.

**Scope narrowed.** "Every experiment in the paper uses k = 31 (or 7 and 8)" is false for the Supplement:
- de Bruijn at k = 47 and 63 (`results/ideas_eisenstein_2026-11-10.md:82-83`);
- palindrome arms longer than 32 bases (paper line 151).

The claim is true for every main-text experiment. The critic also omits that the paper's main-text spaced-seed result (C14) is a regime where canonical-code shortcuts need care. That does not remove the gap.

**Repairs.** A scoping sentence that names what S adds over canonical code plus universal hash (k > 32, rolling field algebra, spaced seeds, values usable directly as sketch elements) is **endorsed and free**. The speed-table row is endorsed but costs one more laila benchmark job. The scope sentence is owed now. The row is optional.

## Gap 3: the sampler's advantage is the heuristic filter, and its baseline is unmatched. NARROWED

**Partly a false gap.** The critic prescribes "state that the selection is heuristic". The paper already says this at line 67: "This filter is strand-symmetric, and its collisions are heuristic." That repair has already been made.

**What survives.** The intercept of 9.96 against 11.37 (recomputed by four lenses) shows that sparse-rate speed comes from the packed filter. "Faster than batched ntHash" (line 106 and the abstract) compares filter-plus-S-on-kept against ntHash on every k-mer. The matched alternative is filter-plus-ntHash-on-kept, which would be faster still. No lens framed it this way: all six C11 findings concern the 1/20 crossing.

**Repairs.**
- One clause is **endorsed and free**: "the guarantee covers collisions among kept k-mers; which k-mers are kept is decided by the heuristic filter, whose uniformity we measured (C13), not proved."
- The matched timing is endorsed but costs a benchmark run.
- The paper's claim is literally true as worded, so the timing is not owed for correctness. The clause is owed.

## Gap 4: "FracMinHash sketches commonly sample far more sparsely than that" is uncited. UPHELD

I confirmed the premise: `grep -i "irber|sourmash|fracmin" manuscript/references.bib` and the ledger return nothing. The repair is endorsed: cite FracMinHash and sourmash and give the default scaled value. The cost is one bib entry and one ledger row, with no re-measurement.

## Gap 5: 61% spurious against "about 3%", and 7 misses. NARROWED

**Partly a false gap.** The critic writes "Nothing explains the 7 misses ... a superset method". The same results line it quotes (`kp10_vs_shibuya_2026-11-11.md:104-105`) says "the listing returned 507,219 of 507,221 every trial". Two unlisted extended syncmers are enough to drop k-mers. The Shibuya summary (line 33) also says peeling "can fail" and the output is a superset only "in general". The misses are explained by peeling, not by a harness defect, although listing the 7 would confirm this cheaply.

**What survives.** The 61%-against-3% discrepancy is real. The results file explains it by regime ("unplaced scaffolds share much sequence with the chromosomes"), but the paper does not report it, and the paper uses that error against the competitor. Reporting both figures is **endorsed and free**. Inspecting the 7 misses is optional and cheap.

## Gap 6: "one error analysis" (Discussion, line 138). UPHELD

The sentence is present: "set operations and reconciliation sketches can run on the hash values with one error analysis". The critic is right that adjudications L3, L5, L6 and L7 established that no derived bound covers the sketch's overflow check, and that no finding named this sentence (grep "error analysis": 0). The delete-or-reword repair is **endorsed; the fix is free**. Use the critic's wording or simply delete the clause.

## Gap 7: unmeasured predictions. NARROWED

- **"A C implementation would hash a genome in seconds"** (line 132). UPHELD. No measurement or ledger source supports it (grep "in seconds" in results and the ledger: 0). Labelling it an estimate is free.
- **"pays only for hash-only graphs with k > 32"** (line 148). Narrowed. As worded this is a necessary condition, not a claim that it does pay. Even so, the laila table (`ideas_eisenstein:82-83`) shows 4.3× and 4.2× slowdowns at k = 47 and 63 against 128-bit packed codes. "Pays" should become "could pay only where packed codes are unavailable". That fix is free.

## Gap 8: C21's inference leg was UNRUN, and the Voss reading. NARROWED

- **"UNRUN".** This is a correct coverage fact (L5 marked C21 UNCHALLENGED). Per the brief, it says nothing about the claim.
- **Voss.** The critic says "the cited paper's headline finding is the opposite of white noise". That overstates the case. The summary (`voss_1992_fractal-dna.md:23`) itself describes "the white-noise region" at high frequency, which is the scale of a 31-base window, so the citation is roughly apt. The real caveat is the period-3 peak inside that region. Quoting Voss's white-noise-region sentence and noting the period-3 peak is **endorsed and free**.

The stakes are low, as the critic says.

## Thin coverage (C5, C6, C7)

This is a correct fact about coverage, and adjudication-L6 already called the C5 search "thin". A second, independent dated search costs one agent pass. It is not free, but it is owed before the priority sentences reach an abstract.

## Section 3 framing

The critic notes "no lens was assigned to baseline adequacy". I endorse this as a structural observation: gaps 1 to 3 all fall there, and the only agent that found gap 1 was the critic. The "gap considered and not filed" (about 6×10⁻⁷ same-bucket collisions per run) is sound reasoning, and I accept it.

| finding | verdict | why |
|---|---|---|
| Gap 1 | UPHELD | KP-10 inputs are min(f,r) canonical and keyed by mix64(canonical); re-run shows a non-symmetric F gives the same invariance, so "because a strand-symmetric hash" is untested by the experiment (S's symmetry itself is tested elsewhere); free fix: disclose canonical inputs, drop "because"; rerun with S-keyed buckets is the costly alternative; correct minority |
| Gap 2 | NARROWED | real baseline/framing gap, but Supplement does run k=47/63; free scoping sentence owed, speed row costs a benchmark |
| Gap 3 | NARROWED | paper already states the filter is heuristic (line 67); unmatched baseline is real; free clause that selection is not covered by the proof; matched timing optional, costs a run |
| Gap 4 | UPHELD | no FracMinHash/sourmash entry in bib or ledger; add citation with default scaled (cheap, no re-measurement) |
| Gap 5 | NARROWED | the 7 misses are explained by 507,219/507,221 listed (peeling is not a guaranteed superset); 61% vs ~3% unreported; free fix: report both with regime |
| Gap 6 | UPHELD | "one error analysis" contradicts the round's rulings that no overflow bound exists; free delete/reword |
| Gap 7 | NARROWED | "in seconds" unmeasured (free: label estimate); "pays only ... k>32" is conditional, but k=47/63 show 4.2-4.3x slowdown; free rewording |
| Gap 8 | NARROWED | UNRUN is a coverage fact; Voss does describe a high-frequency white-noise region, so the citation is roughly apt; the period-3 caveat is the free fix |
