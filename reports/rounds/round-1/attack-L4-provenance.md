# Round 1 attack, slot L4, lens: provenance

Attacked text: `reports/rounds/round-1/ATTACKED-paper.md` (draft d049b77). For each claim, the lens asks one thing: does the number or reading trace to the artifact the ledger names, and does that artifact say what the paper says it says?

**Where things live.** The brief has no appendix, so I found the sources myself:
- Literature: summaries in `bench/repos/library/papers_md/summaries/`, full-text conversions in `papers_md/`, PDFs in `bench/data/papers/`.
- Method cards: `bench/repos/library/methods/`. Two apply here, `substitution_scoring.md` and `reproducibility.md`.
- Prior rounds: `scripts/round history "strand-symmetric"` reports "no prior mention ... in 1 round directories". This is the first round, so nothing here re-files an adjudicated point.
- Scratch scripts are in `/home/agent/.claude/jobs/0d0eac3e/tmp/L4/`.

**Steering note.** The claim list is reasonable, but it leaves out three load-bearing points:
1. Methods line 73 says every timing in the paper was run twice, each in its own counted process, with the load cost subtracted. The spaced-seed timing in C14 meets none of these conditions (F3).
2. The CM31 bound "1.4 × 10⁻¹⁷" in the Implementation section inherits an arithmetic slip from its source (F1).
3. The "61 percent ... spurious" sentence states a ratio against the wrong base (F7).

## F1 -- the quoted bound 2.7×10⁻¹⁷ is 2k/q, not the 2(k−1)/q the body states, and 1.4×10⁻¹⁷ is wrong under either formula

**Claim:** C1. **Verdict:** NARROWED.

**What the paper says.** Methods line 55: "they collide with probability at most $2(k-1)/q$ ... At $k = 31$ that is about $2.7 \times 10^{-17}$."

**What the arithmetic gives.** With q = 2⁶¹ − 1 = 2,305,843,009,213,693,951:
- 2(k−1)/q = 60/q = **2.60 × 10⁻¹⁷**.
- 2k/q = 62/q = 2.69 × 10⁻¹⁷.

The 2.7 comes from `results/kp7_eisenstein_2026-11-10.md:29-30` ("≤ about 2k/q ≈ 2.7 × 10⁻¹⁷"). The abstract keeps the source's "about 2k/q". The body changed the formula to the correct total-degree bound but kept the number that belongs to the other formula.

**The same slip in the CM31 bound.** Implementation line 65 says "about $1.4 \times 10^{-17}$". Its source, `results/kp7_m31_simd_2026-11-10.md:13`, says "≤ about 2k/p² ≈ 1.4e-17". With p = 2³¹ − 1:
- 2k/p² = **1.34 × 10⁻¹⁷**.
- 2(k−1)/p² = **1.30 × 10⁻¹⁷**.

Neither rounds to 1.4. Only 64/p² does, which would be k = 32. The ledger copies the source and the source is wrong. The same 1.4 appears in `results/kp10_vs_shibuya_2026-11-11.md:43`.

**What survives.**
- The total degree 2(k−1) is right: the λ·(F+R) term has degree k and the F·R term has degree 2(k−1).
- The exhaustive checks reproduce from an independent implementation. My `eis.py` uses a different random x and λ over F_q with ω = g^((q−1)/3). At k = 7 and k = 8, S takes 8,192 and 32,896 values, and F takes all 16,384 and 65,536. This matches `hash/mutant_e.out` and `results/kp7_eisenstein_2026-11-10.md:33-34`.

**Correction.**
- Body: "at most $2(k-1)/q$ ... about $2.6 \times 10^{-17}$". Alternatively, write 2k/q in both places.
- Implementation: "about $1.3 \times 10^{-17}$".
- Fix the source notes and the ledger rows.

## F2 -- F·R merges more than a k-mer and its reverse; 16,362 is not the count the stated mechanism predicts

**Claim:** C2. **Verdict:** NARROWED.

**What reproduces.** The table's counts come from `prototype/kp7_symmetric_check.py`: fourth-root map, F_{q²} with q = 2³¹ − 1, x complex, seed 1. Its output `results/kp7_symmetric_check_2026-11-10.txt:1-10` gives 729 / 4,160 / 8,192 and 6,561 / 16,362 / 32,896. The caption's encoding and field are right.

**The F·R count does not match the stated mechanism.** The text says F·R "cannot tell a k-mer from its plain reverse". If that were the only merge, the count would be the number of orbits under {identity, reverse, complement, reverse-complement}. By Burnside:
- k = 7: (16,384 + 256)/4 = **4,160**. This matches.
- k = 8: (65,536 + 256 + 256)/4 = **16,512**. The table gives 16,362.

So at k = 8, F·R merges 150 more values than reversal explains. These are structural: random collisions over F_{q²} are about 10⁻⁹ expected. The excess also depends on the encoding. Under the Eisenstein map I get **16,500** (`eis.py`: `8 6561 16500 32896 65536`).

**The abstract's 729 has no committed Eisenstein artifact.** The abstract places "the sum F+R takes 729 values" in the Eisenstein context, but the committed output is fourth-root. I reproduced 729 under Eisenstein as well, so the number holds. It still needs its own artifact.

**Correction.**
- Write: "$F \cdot R$ is an autocorrelation, so it cannot tell a k-mer from its plain reverse (or from other k-mers with the same autocorrelation)".
- Commit an Eisenstein run of the three functions.

## F3 -- the timings were not taken at 2.0 GHz, and the spaced-seed timing breaks the paper's own measurement statement

**Claims:** C10 and C14 (timing part). **Verdict for C10:** NARROWED.

**The clock was not 2.0 GHz.** `results/paper_bench_2026-11-11/launcher.log` records cycles and task-clock for every run. Cycles divided by task-clock:
- The `__load_only__` runs (run_25, run_49, run_46, run_22) ran at **2.22–2.37 GHz**.
- Compute runs ran at up to **2.98–3.00 GHz** (run_4, run_28, run_51).
- The headline runs: run_2 at 2.57 GHz, run_6 at 2.85, run_30 at 2.76.

Turbo was on. Two texts are therefore wrong:
- The abstract's "On a 2.0 GHz Xeon".
- `results/kp7_eisenstein_2026-11-10.md:59`: "The counter cycles are core clock cycles at the 2.0 GHz base frequency".

The cycle counts themselves look sound. Load-only cycles stay at 6.96–7.01 × 10⁹ while frequency moves between 2.22 and 2.47 GHz, so subtracting the load cost is not biased by frequency. I reran `manuscript/paper_bench_table.py` on the raw directory and the output is byte-identical to `table.md`. Every cycle number in C10 traces.

The problem is the label. A reader who turns "24.7 cycles on a 2.0 GHz Xeon" into 12.3 ns gets a wall time about 40% too slow.

**The spaced-seed timing breaks Methods line 73.** Methods says: "Every timing in the paper comes from one job, in which each method ran in its own counted process ... twice, in interleaved order. The cost of loading ... subtracted." The 0.76× in Results line 114 comes from somewhere else. `results/paper_bench_2026-11-11.md:54` describes it as "run_51; bench_spaced_e's own wall-clock rates, first 60,000,000 bases of cs10, 3 rounds". `run.log` line 344 confirms it ran once, with both methods in one process, timed by wall clock, with no load subtraction.

**Correction.**
- Drop "2.0 GHz", or say "turbo enabled; cycles are core cycles".
- Change Methods to "Every timing except the spaced-seed comparison, which is a wall-clock rate from one run of `bench_spaced_e` in the same job, ...".

## F4 -- the bit-for-bit check is a checksum plus a 2-million-k-mer prefix, it is not in the paper's job, and no output was committed

**Claim:** C12. **Verdict:** NARROWED.

**The check is not in the paper's job.** `hash/bench_e.cpp` is the binary that produced every ntHash number in the paper. Its identity pairs (lines 164-166) have **no NT_BATCH vs NT_LIB pair**. Every run in the job used `only=...`, and such runs return at line 162 before any identity check.

**The check that exists.** It is in `hash/bench_pack.cpp:276-295`. It compares the XOR of all hashes plus the count over the whole input. It compares position and hash element by element only for the first 2,000,000 k-mers of the first record.

**No output was committed.** No output of this check exists in `results/`, in `results/laila_counters_*`, or in `bench/data/laila_perf/` (grep for "positions and hashes identical": no hits).

**I reran it locally.** On `cs10_head6m.fa` (5,882,341 bases), `bench_pack` prints "identical to ntHash canonical (fwd+rev): yes" and "identical on the first 2000000 k-mers of NC_044371.1: yes". So the claim is probably true. As written, though, it rests on a check that was neither run in the job nor recorded.

"About twice": 23.53/11.37 = 2.07. This holds.

**Correction.**
- Commit the `bench_pack` identity output on full cs10, at ntHash c26bd45.
- Describe the check as "the same XOR checksum and count over all 736.5 M k-mers, and identical (position, hash) pairs on the first 2 M".

## F5 -- at exactly 1/20 the paper's own model says the sampler is slower

**Claim:** C11. **Verdict:** NARROWED.

`results/paper_bench_2026-11-11.md:47` gives the batched model: slope 28.5 and intercept 9.96. At a rate of 1/20 this predicts 9.96 + 28.5/20 = **11.385 cycles**, against batched ntHash's **11.37**. The break-even is 1/20.3, which is sparser than 1/20.

So the abstract's "faster than batched ntHash at sampling rates of 1/20 and below" includes 1/20, where the model predicts slower. The margin is 0.13%, well inside the 0.7–0.8% repeat spread of the two rows it is computed from (lines 22 and 36). The line fits the measured batched rows at 1/100 (10.22 measured, 10.245 predicted) and 1/1000 (9.95 measured, 9.99 predicted). The break-even itself is an extrapolation past the 1/16 endpoint into a region where the two methods cannot be told apart.

**Correction.** Abstract: "at sampling rates below about 1/20", which matches Results line 106. Better: "at 1/100 and sparser (measured, 1.14×)".

## F6 -- the 0.747 vs 0.818 test measures no strand consistency; it shows our hash turns a weight-20 seed into a weight-29 seed

**Claim:** C14. **Verdict:** REFUTED, for the comparison sentence and its framing. The algebraic symmetry holds.

**The mask and its mirror.** The non-palindromic mask is `1101100111010111011100101110101` (`ideas/eisenstein/bench_spaced_e.cpp:37`). It has weight 20. Its union with its mirror image has 29 positions.

**Both numbers are exact powers of 0.99.**
- 0.99²⁰ = 0.81791. SeedNtHash's measured value is 0.81794 (`results/paper_bench_2026-11-11.md:62`).
- 0.99²⁹ = 0.74717. Ours is 0.74722.

So SeedNtHash behaves exactly as the weight-20 seed it was given. Our strand-symmetric product depends on every position in mask ∪ mirror(mask), so it behaves as a 29-position seed.

**The "exact double-stranded value" is defined to match our construction.** `bench_spaced_e.cpp:149-155` defines "exact" as "care bases of s and of rc(s) unchanged as a set", which is the union. It agrees with our hash by construction.

**No reverse complements are compared.** The test compares s with a mutated s on the same strand (lines 139-147). It cannot show the inconsistency ntHash2 warns about, which is hash(s) ≠ hash(rc(s)). The warning itself is real: ntHash c26bd45 `src/seed.cpp:96-101`, "is not symmetric, reverse-complement hashing will be inconsistent".

**Effect.** Placed directly after that warning, the paper's sentence presents a sensitivity loss of 0.818 → 0.747 as a consistency win. For a spaced seed, the fraction of windows that survive substitutions is the figure of merit.

**Correction.**
- State the cost: "for a non-palindromic mask the symmetric hash depends on the mask and its mirror (29 of 31 positions here), so at 1% substitutions 0.747 of windows survive, against 0.818 for ntHash2's weight-20 seed, which is not strand-consistent".
- Test strand consistency directly: hash(s) = hash(rc(s)) for both methods.
- See F3 for the timing's provenance.

## F7 -- "61 percent of it spurious" uses the wrong base, and "equal reliability" compares a two-round protocol with a one-shot one

**Claim:** C15. **Verdict:** REFUTED, for these two statements. The other numbers trace.

**The 61%.** Results line 130 says the extended-syncmer method "returned a superset, 61 percent of it spurious on this pair". The source says "1,986,807 spurious k-mers (61% of D)" (`results/kp10_vs_shibuya_2026-11-11.md:105`, and `full_kmp_ext_giga.json`: missing 7, spurious 1,986,807 in all 5 trials).
- The returned set holds 3,238,094 true differences plus 1,986,807 spurious ones, 5,224,901 in total.
- The spurious share of what it returned is **38.0%**. 61% is spurious ÷ true difference.
- It is not strictly a superset either: 7 true differences are missing. The table says so; the prose does not.

**Equal reliability.**
- `full_kp10_eis_giga.json` (and the Gaussian file) records `rounds_max 2` and `splits 1275`. The 36,494,096 bytes equal (65,536 + 2 × 1,275) × 67 × 8, so the resend is counted.
- So at 11.3 B/diff, 1,275 of 65,536 buckets (1.9%) overflowed. The sketch reached exactness only through a second round.
- The k-mer IBLT's 20.5 B/diff is one-shot, without retry, at 30/30 trials, measured on partitions 0–3 only (`full_kmp_kmers_giga_r1.05.json`: parts [0,1,2,3], D = 202,557). It was also sized from the oracle true D (`recon/kmpeeler_run.py`, n = ratio × D).
- With the same allowance, one retry on failure, the IBLT at n = D (19.6 B/diff, 15/30) is the like-for-like comparison.

The byte accounting itself is consistent. The reliability regimes are not equal.

**Correction.**
- "returned 1,986,807 spurious k-mers, 38 percent of what it returned (61 percent of the true difference), and missed 7".
- Replace "of equal reliability" with "succeeding in 30 of 30 one-shot trials on 4 of 64 partitions; our sketch needed a second round for 1.9 percent of buckets". Change the abstract and Discussion ("0.55 times") the same way.

## F8 -- "detected deterministically" contradicts the cited source; the check is probabilistic over the hash

**Claim:** C16. **Verdict:** REFUTED.

**What the decoder accepts.** The decoder (`recon/recon61.py:10-11`) accepts a bucket when two things hold:
- the flagged roots reproduce the c extra power sums;
- the flagged counts match |A_b| − |B_b|.

**How a wrong answer can pass.** Take a bucket with n > m differences. It passes if some set of at most m genuine elements has its first m + c power sums equal to those of the true n-element difference. That is a nonzero polynomial identity in the random (x, λ). It fails with probability bounded by the degree over q. It does not fail always.

**The cited source says the same.** The paper cites Minsky and colleagues for this step (Methods line 59). The library summary (`minsky_2003_set-reconciliation.md:23-25`) says equality "is tested probabilistically, with error per check at most (D − 1)/|E|". Hash collisions inside a set also defeat any power-sum check. The Methods sentence ("is detected rather than resolved wrongly") and the Discussion ("detected its own overflow deterministically") both overstate.

**Correction.** "detected except with probability about (m + c)·2k/q per bucket, over the random hash".

## F9 -- the Ts/Tv range mixes subsets, the alignment set is mis-described, and the source alignments violate the substitution-scoring card

**Claim:** C9. **Verdict:** NARROWED.

**The arithmetic is right.** `results/ideas_eisenstein/kp1_tstv_NC_044378.1.json` gives Ts/Tv of 1.585, 1.584, 1.669 and 1.847. The transition shares, 61.3–64.9%, follow from them (1.585/2.585 and 1.847/2.847).

**The description is not.** "Pairwise alignments of cannabis assemblies" here means three comparisons:
- Two queries (KOMPa and KOMPb) aligned against cs10. They come from a WF-I panel that holds dozens of labels; `ls .../NC_044378.1/bare/` shows 79X, ABC, ACBD, and so on.
- One KOMPa vs KOMPb comparison. It is not a pairwise alignment. It is projected through cs10 coordinates (`kp1_tstv.py` docstring; `ideas_eisenstein_2026-11-10.md:20-21`).

The 65% / 1.85 endpoint comes only from the **uppercase-only** subset of that projected haplotype comparison. Over all columns the range is 61.3–62.5% (Ts/Tv 1.58–1.67).

**The method card.** `library/methods/substitution_scoring.md` decision 3: "An estimate taken from alignments produced under matrix M is not independent evidence about M." The MAFs are lastz output: `bench/data/wfi_laila/NC_044378.1/bare/KOMPa.maf` line 1 reads `##maf version=1 scoring=blastz`. A blastz-family matrix scores transitions far above transversions. The paper then uses this Ts/Tv to argue that its transition-favouring metric fits real substitutions. The project's own card says this is partly circular. The circularity is probably small at near-identity, but the card requires saying so.

**Correction.** "Between cs10 and two other assemblies on NC_044378.1, 61 percent of aligned substitutions are transitions (Ts/Tv 1.58), rising to 65 percent (1.85) between two haplotypes outside soft-masked repeats; the alignments were made with a transition-favouring scoring matrix."

## F10 -- the Marçais reading hardens "can negate" and erases a table that moves both ways

**Claim:** C19. **Verdict:** NARROWED.

**"Can negate", not "negates".** Discussion line 140 says canonicalization "negates the window guarantee of syncmer-style schemes". The source says it "can negate" (`marcais_2024_knonical.md`, quoting lines 89-90) and "may not preserve" (line 84).

**Fractional sampling moves both ways.** The paper says there is "no visible effect on fractional sampling". The library summary's rebuild of Table 2 from the PDF shows fractional deserts moving in both directions at k = 31:
- chromosome 1: 63.20 → 38.14 Mb;
- chromosome 10: 24.68 → 57.52 Mb, a 2.3× growth.

The source also reports a standard deviation of 34.7 Mb. "No visible effect" is the repo note's gloss (`notes/kp7_strand_symmetric_hash.md:45`) of the caption's "every method except fractional". It is not a measured null.

**Version.** The summary was read from the bioRxiv preprint. The paper cites the Bioinformatics DOI, which "was not the copy read".

**Correction.** "show that folding ... can negate the window guarantee of syncmer-style schemes, while fractional sampling shows no systematic change (deserts shrink on some chromosomes and grow on others, with high variance)". Cite the preprint, or check the journal version.

## F11 -- ntHash's canonical value is not GF(2)-linear, and combining the two strands predates ntHash2

**Claim:** C18. **Verdict:** NARROWED.

**What holds.**
- Fixed seeds: confirmed. See `hash/ntbatch.h` SEED_A..T, copied from ntHash, and `MULTISEED` in ntHash `src/internal.hpp:95`.
- The 8-mer test exists: `results/kp7_rolling_vs_nthash_2026-11-10.raw.txt:4`, "32896 distinct values over 65536 k-mers".

**Two defects.**
- **Linearity.** ntHash's canonical value is `fwd + rev`, an integer sum mod 2⁶⁴ (ntHash c26bd45 `src/internal.hpp:27-29`; ntbatch.h line 75). Integer addition is not GF(2)-linear. Only the per-strand XOR-of-rotations values are. The properties the Discussion lists (constant-time substitution, concatenation) hold for ntHash by working on the (fwd, rev) pair and then adding, not because "ntHash ... is linear over GF(2)".
- **Attribution.** "Strand symmetry by combining the two strand hashes is ntHash2's idea." ntHash 2016 already combined them, by minimum. ntHash2 says it replaced "the minimum of the hashes, which was used in the previous version" with the sum (`kazemi_2022_nthash2.md:13`). ntHash2's contribution is the sum.

**Correction.**
- "whose per-strand hashes are linear over GF(2)".
- "Combining the two strand hashes is ntHash's idea, and the symmetric sum is ntHash2's."

## F12 -- two Supplement numbers are described more broadly than their sources

**Claim:** C21. **Verdict:** NARROWED.

**The de Bruijn test was not on a chromosome.** The supplement says the test was exact "on 13.6 million edges of a cs10 chromosome's graph". The source ran `dbg_e{31,47,63}` "on the first 20 Mb of cs10" (`results/ideas_eisenstein_2026-11-10.md:77`). That is part of one chromosome, not a chromosome's graph.

A local run disagrees on the slowdown. `results/ideas_eisenstein/dbg_e31.txt:2-4` gives 570.7 vs 971.4 ns/node, a 1.7× slowdown. The laila table gives 5.9× (204.9 vs 1,210.7). The "four to six times" depends on the machine, and the paper reports only the laila value.

**The 42% is not a share of 31-mers.** "Two haplotypes ... differ in 42 percent of their chromosome 7 31-mers" uses D/|A| (`kp9_kompa_kompb_2026-11-10.md:41`, D = 16,316,892 = 0.42 |A|), and D counts both sides. The share of the union's 31-mers that differ is 16.32/(30.42 + 16.32) = **35%**.

**What traces.** 7,760; 1,370× (`ideas_kp5_kp11:74`, 0.81 s vs 1,108 s); 0.210 / 0.766; 0.971–0.991 / 0.635–0.576.

**Correction.** "exact on 13.6 million edges from the first 20 Mb of cs10"; "their symmetric difference is 42 percent of one haplotype's 31-mers (35 percent of the union)".

## F13 -- the readings of the literature check out against the PDFs, but two paraphrases tidy their sources and one citation is unread

**Claim:** C7. **Verdict:** NARROWED.

**Verified against full text and PDF:**
- Cristea, equations (3)–(5): `papers_md/cristea_2002_genomic-signals.md` lines 770-835. Equation (5) is the fourth-root map up to sign. Equation (4) scaled by 1/√2 equals the Das family at θ = π/4.
- Das 2017, equation (2) with "the ș value is taken as pi/3" (conversion lines 122-128). Their map is ours times ω, so the two are related by a rotation.
- The printed 97.36% = 1218/1251, their Eq. 8 (specificity).
- Three ranges for the second exon. In `pdftotext -layout` output: text 2540-2867 (line 187), Table 1 2599-2906 (line 412), Table 3 2649-2906 (line 446).
- Anastassiou: t = a*, g = c*, so complement is conjugation.
- Mendizabal Table 6: real and paired-numeric spectra are identical for complement and reverse complement.

**Defects:**
- **(a) Das's rationale as printed is false for the bases it names.** Das writes "Purines -G and pyrimidines-C have identical imaginary but opposite real parts" (conversion lines 131-132). Under their Eq. 2, G = ½ + j√3/2 and C = −½ − j√3/2, whose imaginary parts are opposite. The paper's paraphrase ("purines share imaginary parts, as do pyrimidines") silently corrects this into Cristea's true statement, then reports that Das "restates Cristea's".
- **(b) Cristea's purine/pyrimidine sentence belongs to a different map.** It is attached to equation (3), the conjugation map (lines 777-783). Equation (3) is identical to Anastassiou's equation (2). The sentence is not attached to the negation maps (4) and (5).
- **(c) Cristea did not use θ = π/4.** The "θ = π/4 point" needs a 1/√2 scaling that Cristea applied only on the way to (5), not to (4).
- **(d) "Reviewed by Das and colleagues" is unsupported for Cristea.** `das_2021_dna-encoding` never mentions Cristea (grep of summary and conversion: no hits). It reviews the trigonometric family and attributes it to Das 2017.
- **(e) An unread paper is cited.** The paper cites `[@das_2019_kaiser-exon]`, which it says it "has not read". The "later comparison" that cites it is `gulocak_2022_mapping-effect`, which is held in the library and not cited (`literature/eisenstein_wantlist.md:93-96`).

**Correction.**
- Quote Das's sentence as printed and note that it pairs G with C.
- Say that Cristea's purine/pyrimidine remark is made for his conjugate map (3).
- Attach the das_2021 citation to "the trigonometric family" explicitly.
- Cite Gulocak 2022 for the 2019 attribution rather than citing the unread paper as a source.

## F14 -- "tens of long gaps per seed" fails for one ntHash seed

**Claim:** C13. **Verdict:** NARROWED.

**What traces.**
- 1,774 = 341 + 363 + 371 + 332 + 367 (`results/kp7_stats_seeds_eisenstein_2026-11-11.md:19-23`).
- 367 long gaps for the Eisenstein hash, 0 outside repeats.
- 0.00090 vs 0.00088.
- The 0.14 reproduces: 736.5 M × (2/21) × (19/21)²⁰⁰ = 0.144.

**What does not.** "Every method left tens of threshold-sampling gaps longer than 200 k-mers per seed" fails for ntHash seed 1, which left **8** (`results/kp7_next_steps_2026-11-10.md`, table: 8, 87, 56, 94, 96).

**Definition.** "In soft-masked repeats" means that at least half of the gap's spanned bases are soft-masked (`kp7_next_steps:60`). The genome is 47.7% soft-masked, so the paper should state this definition.

**Correction.** "Every method left up to tens of such gaps per seed (8 to 96)". Define "in a repeat".

## F15 -- the absence claims rest on dated searches that the paper never mentions

**Claims:** C5 and C6. **Verdict:** NARROWED for both.

**The searches exist but are not in the paper.** The repository holds dated negative searches:
- `literature/wantlist.md:69-83`, OpenAlex and web, 2026-11-10, including "polynomial rolling hash reverse complement DNA";
- `literature/eisenstein_wantlist.md:57-83`, including "k-mer hashing finite field reverse complement symmetric" and "DNA numerical representation transitions transversions complex".

The paper states both novelty claims ("Our contribution here is narrower...", "We contribute four things") with no dated qualifier. Under this project's absence-claim rule, such a claim is unmeasured.

**Held works the paper does not discuss.**
- The novelty search found `benbassat_2014_incremental-hashing`, a Karp–Rabin polynomial rolling hash over DNA, held in the library.
- `chen_2026_hash-functions`, the current survey of nucleotide hashing (canonical encodings, ntHash), is also held.

Neither is cited. I do not claim either pre-empts C5 or C6; that is the prior-art lens's question. The provenance defect is that the paper's absence claims do not point at the searches that support them.

**Correction.** "As far as searches of OpenAlex and the web on 2026-11-10 found, ..." for both contributions. Cite or dismiss chen_2026.

## Brief notes on the remaining claims

**C3: SURVIVED.** I derived R_{j+1} directly. R_j = −Σ_{u=0}^{k−1} s_{j+u} x^{k−1−u}, so x(R_j + s_j x^{k−1}) − s_{j+k} = −Σ_{u=1}^{k} s_{j+u} x^{k−u} = R_{j+1}. The F rule is immediate. The ROLL mutant is killed (`kp7_eisenstein_2026-11-10.md:49-50`).

**C4: NARROWED.**
- The congruences hold: 2⁶¹ ≡ 2 (mod 3), so q ≡ 1; q ≡ −1 ≡ 3 (mod 4).
- The cycle ratio is 47.37/24.69 = 1.92.
- The paper says "about twice as many *multiplications*". The project's own count (`kp7_eisenstein_2026-11-10.md:19-21`) is 3 for Eisenstein against "4 plus a full F_{q²} product", which is 7 or 8, a ratio of 2.3–2.7×. Instructions per k-mer are 179.3 / 78.5 = 2.3×. Say "about twice the cycles".

**C8: SURVIVED.** |1−ω|² = 3, |1+ω|² = |ω²|² = 1, |2|² = 4. The sentence and the Discussion state that the field hash keeps no metric.

**C17: NARROWED.**
- 137 KB = 2⁸ × 67 × 8 = 137,216 bytes. This traces (`full_kp10_eis_draft.json`).
- No artifact measures the draft contigs' orientation. The cause "because the draft's contigs lie in other orientations" is inferred.
- Canonicalizing the extended-syncmer strings still leaves 33,446,895 differences (`kp10_vs_shibuya:112`). That shows the method's dependence on orientation runs through syncmer selection, not just string orientation.
- I checked whether soft-mask case could explain the 48M instead: `kmpeeler_run.py` applies `tr a-z A-Z`, so it cannot.
- Either measure orientation (count the draft contigs that map reversed) or write "consistent with".

**C20: NARROWED.**
- "Block sizes from 1 to 4,096" means the discrete sets {1, 3, 5, 64, 4096} or {1, 3, 64, 4096} (`hash/test_kp7e.c:149,197`). CM31 uses 11 sizes.
- The mutants are killed (`hash/mutant_e.out`: "FAILED: 64072 check(s)").
- One-sided recon tests exist in `tests/test_kp10.py`.
- Write "block sizes from 1 to 4,096 (1, 3, 5, 64, 4,096)".

## Verdicts

| claim | verdict |
|---|---|
| C1 | NARROWED |
| C2 | NARROWED |
| C3 | SURVIVED |
| C4 | NARROWED |
| C5 | NARROWED |
| C6 | NARROWED |
| C7 | NARROWED |
| C8 | SURVIVED |
| C9 | NARROWED |
| C10 | NARROWED |
| C11 | NARROWED |
| C12 | NARROWED |
| C13 | NARROWED |
| C14 | REFUTED |
| C15 | REFUTED |
| C16 | REFUTED |
| C17 | NARROWED |
| C18 | NARROWED |
| C19 | NARROWED |
| C20 | NARROWED |
| C21 | NARROWED |
