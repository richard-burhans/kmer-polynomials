# Round 1 adjudication of L7 (reproduction)

Adjudicator: independent skeptic, L7 slot. The attack file `attack-L7-reproduction.md` hashes to
`15205881f37207eaea...`, which matches the sealed digest `15205881f372`, so I adjudicated the sealed attack.

**History (rule 5).** `scripts/round history "overflow" --kind adjudication` and `... "kmer-polynomials"` both return
a searched absence. Round 1 is the first round, so none of these findings re-files a settled point.
**Method cards (rule 6).** `refs/methods/` does not exist in this repository, so no method card governs any dispute
below. Where a statistical point arises (F1, F6), I rule from the raw data.

What I re-ran myself, rather than taking the attacker's word for it:

- `tmp/L7/be.py` on `results/paper_bench_2026-11-11`: break-even pairs 1/19.1 to 1/21.7, mean 1/20.33, model at
  r = 1/20 gives 11.388 against ntHash's 11.365. All of this matches the attack. I also checked the fit against the
  third batched row, which the attacker did not use. 9.961 + 28.53/100 = 10.246 is predicted, and 10.204 / 10.237
  was measured.
- `tmp/L7/overflow.py` (20,000 trials each): q=31,c=1 gave 62 wrong answers accepted; q=31,c=2 gave 3; q=101,c=1
  gave 17. These are the same concrete counterexamples the attack prints. I read `recon/recon61.py::decode` (lines
  100-136), and it re-checks only S_{m+1}..S_{m+c} and the count.
- `round1_L7/nteq` on cs10: 221 records, 736,547,799 k-mers on both sides, 0 mismatches.
- `round1_L7/gf2`: the canonical value matched 0/10000, forward 10000/10000 and reverse 10000/10000.
- `kp1_tstv.py` on the ACDC and ANC-1 labels: cs10 vs ACDC 1.511 (60.2%), cs10 vs ANC-1 1.546 (60.7%). Both match
  the attack.
- Raw JSONs of `results/kp10_vs_shibuya_2026-11-11/`: `rounds_max: 2`, `splits: 1275`, 4,561,762 elements =
  65,536 x 67 + 1,275 x 134. The IBLT's 30/30 comes from `parts: [0,1,2,3]`. All of this is confirmed.
- Arithmetic: 60/q = 2.60e-17, 62/q = 2.69e-17, 62/p^2 = 1.34e-17, 60/p^2 = 1.30e-17, 64/p^2 = 1.39e-17. Also
  (19/21)^200 x 70.1 M = 0.142, and q mod 3 = 1, q mod 4 = 3.
- `run.log` line 346 (spaced seeds), `results/kp7_stats_seeds_2026-11-10.eisenstein.tsv` (ntHash seed 0: 8 long
  gaps), `hash/test_kp7e.c:149,197` and `test_kp7pack.c` (block-size sets), `manuscript/run_paper_bench.sh:35`
  (build line, no tests).

## F1 (C11, the 1/20 break-even): NARROWED

**Diagnosis upheld.** The paper's own fit puts the break-even at 1/20.3, so at exactly 1/20 the model predicts the
sampler slightly *slower*: 11.388 against 11.365 cycles. The individual repeats spread the break-even over 1/19.1 to
1/21.7. Treating 1/20 as a sharp boundary ("1/20 and below") is not supported.

**Premise partly overturned.** The fit is not a fragile extrapolation. Both fitted rows come from the same batched
mode. (The claim list's worry about "DIFFERENT evaluation modes" does not apply to the fit.) The line predicts the
unused 1/100 batched row within 0.3%. The linear model is therefore validated by a third point. The only real issue
is that 0.2% of predicted difference at 1/20 sits inside the 1.3% repeat spread.

**Repair overturned.** The attacker's fix ("faster at 1/100 and sparser (measured)") throws away a validated model
and installs a claim weaker than the evidence supports. The Results already says "about 1/20 and below". The
correct repair is to make the abstract match that: "at sampling rates of about 1/20 and below". Optionally add
"(interpolated; 1/19 to 1/22 across repeats)" in Results. **The fix is free.**

## F2 (C10, "every timing ... one job ... counters"): NARROWED

**Defect upheld.** Methods says every timing comes from performance counters in one job, one counted process per
method. Several timings in the paper break that:

- KP-3 "four to six times slower" comes from job af1eb52f65c6d471.
- KP-11 "1,370 times faster" comes from job 9f793f18b3112141.
- The spaced-seed 0.76x is `bench_spaced_e`'s wall clock, with three methods in one process.

The 0.76 is also a time ratio over unequal hash counts. `run.log:346` gives 42,088,079 hashes in 1.953 s against
42,667,269 in 1.482 s, which is 0.7485 per hash, i.e. 0.75x.

**Scope narrowed on two legs:**

- The spaced-seed run IS in job 08653d9e6a7453a2 (run_51 in that job's `run.log`). The paper's "In the same job" is
  true. Only "counters" and "own process" fail for it, so it needs a method label, not a different job label.
- The attacker invokes the ledger's no-cross-job-comparison rule against the Supplement timings, but that rule is
  misapplied. Each Supplement ratio (5.9x, 1,370x) is a within-job ratio. Nothing is compared across jobs. The
  defect is only the false "every timing" sentence.

**Repair endorsed with that correction**, and it is free:

- Methods: "Every timing in the speed table comes from one job ...".
- Spaced seeds: "wall-clock rates in the same job".
- 0.76 -> 0.75 (per hash).
- Supplement items get their own job ids.

## F3 (C12, bit-for-bit batched ntHash): UPHELD (claim survives; evidence gap real)

I re-ran `nteq`: 0 mismatches over all 736,547,799 k-mers, positions and values, against the pinned ntHash. The
ratio 23.53/11.37 = 2.07 supports "about twice". The attacker's secondary point also holds:

- The repository's only check (`hash/bench_pack.cpp:274-296`) is XOR plus count genome-wide, and position-by-position
  on the first 2 M k-mers of one record.
- No results file records that check's output.
- The paper-bench job does not run it.

The claim is true, and the paper should point at a committed full-genome check. **Cost:** commit `nteq.cpp` and its
output, about 20 s to run. No text change is needed.

## F4 (C9, Ts/Tv 1.58-1.85): NARROWED

**Defect upheld.** The range takes its minimum from the all-columns rows and its maximum from the uppercase-only
row. The all-columns range alone is 1.58-1.67 (61.3-62.5%). I reproduced two of the four extra pairs (cs10 vs ACDC
1.511, cs10 vs ANC-1 1.546). These fall outside the stated range, so the range is a property of the alignment set
chosen, not of "cannabis assemblies".

**Attacker's framing narrowed.** "One plant's two haplotypes" overstates. Two of the three rows are cs10 against
KOMPa and against KOMPb. Those are genuine pairwise alignments between two different plants' assemblies. Only the
third row is haplotype-vs-haplotype, projected through cs10.

**Repair.** Either of the attacker's two options is accurate. The cheaper honest one is free: "In pairwise
alignments of cs10 to the two haplotypes of one other cultivar on NC_044378.1, 61 to 65 percent ... (Ts/Tv 1.58
with all columns to 1.85 outside soft-masked repeats)". The broadened 60-66% / 1.5-1.9 version needs the four extra
runs committed. They take seconds; I timed one at 2 s.

## F5 (C16, deterministic overflow detection): UPHELD on the defect; NARROWED on the repair

**Diagnosis upheld.** This is the strongest finding of the attack. I reproduced the counterexample with the
repository's own `decode` logic on toy fields: q = 31, m = 2, c = 1, true {18,29,30}/{13,22,23} is accepted as
{}/{}. `decode` checks only the c extra sums and the count. Whether a wrong answer passes therefore depends on the
hash values, and "detected deterministically" (Methods "detected rather than resolved wrongly", Discussion,
`kp10_vs_shibuya_2026-11-11.md:97,124`) is false as a mathematical statement. The practical claim survives at
q = 2^61 - 1. No run is in doubt: every recovery was compared element by element (C20).

**Two parts of the prescribed repair are overturned.**

1. *"Comparing against all m+c sums would cost nothing"* implies the full check would make detection deterministic.
   It would not.
   - A wrong answer E (|E| <= m) passes a full check exactly when the signed set D minus E has vanishing power sums
     for j = 0..m+c.
   - By the Vandermonde argument, that needs at least m+c+2 distinct elements in D xor E. Any overflow with
     |D| >= c+2 can reach that.
   - So the full check closes small overflows deterministically, but not large ones. It is a worthwhile hardening,
     not a cure.
2. *"Except with probability about 1/q^c"* installs an unproven, stronger bound.
   - It treats the hash values as independent uniform field elements. They are not: they are polynomials in only two
     random parameters (x, lambda).
   - With two random parameters, Schwartz-Zippel gives a bound of the form degree/q. It does not give 1/q^c in
     general.
   - The attacker's own toy table does not fit 1/q^c cleanly either. At q = 31, c = 1 it gives 0.0033, against
     1/q = 0.032.

**Endorsed wording:** "is detected except with small probability over the random choice of x and lambda". Drop
"deterministically" in the Discussion and the results-table row. **This fix is free.** Proving a specific bound is
not free, and the paper should not quote one without the derivation.

## F6 (C15, "of equal reliability"): NARROWED (attack upheld in substance)

Every number reproduces, and the bytes are accounted the same way. KP-10's 11.27 includes its round-2 resends,
which is to its credit. But the two methods do not have equal reliability:

- KP-10 is a two-round protocol: 1,275 of 65,536 buckets were split and resent.
- The IBLT is one round, 30/30 on partitions 0-3 only, and 26/30 at the same 1.05 x D on the draft pair.
- Both were sized from the true difference.

So "equal reliability" compares unequal regimes. The attacker's point that an IBLT could also be resent on failure
is fair.

**Repair endorsed, shortened**, and it is free: "11.3 bytes per difference over two rounds, against 20.5 for a
one-round invertible Bloom lookup table that succeeded in 30 of 30 trials on 4 of 64 partitions; both were sized from
the true difference". Make the same change to the Discussion's "0.55 times ... of equal reliability".

## F7 (C18, "linear over GF(2)"): NARROWED (attack upheld)

I re-ran `gf2`: forward and reverse each satisfy the affine test 10000/10000, and the canonical f + r (with integer
carries) satisfies it 0/10000. The paper's conclusion ("not advantages") stands, because the operations do hold
through ntHash's two strand states. The stated reason is wrong for the value compared.

**Repair endorsed:** "... hold for ntHash too, through its two strand hashes, each linear over GF(2)". **The fix is
free.** The 8-mer test exists, which the attacker confirmed. I did not re-run it, and nothing in the finding depends
on it.

## F8 (C1, bound arithmetic): NARROWED (attack upheld)

The arithmetic is verified. Methods writes 2(k-1)/q, which is 2.60e-17, but quotes "about 2.7e-17", which is 2k/q.
The Implementation section's "1.4e-17" matches neither 2k/p^2 = 1.34e-17 nor 2(k-1)/p^2 = 1.30e-17.

The whole set (rule 3) needs checking, because the same 1.4 appears in several places:

- `results/kp7_m31_simd_2026-11-10.md:13`
- `manuscript/ledger.md:17`
- `results/kp10_vs_shibuya_2026-11-11.md:43,129` (including the D table row "Gaussian 1.4 x 10^-17, Eisenstein
  2.7 x 10^-17")

All of them must change together. The abstract's "at most about 2k/q, 2.7e-17" is a looser, true bound and can stay.
The rest of C1 survives: the injectivity argument, the degree 2(k-1), and the exhaustive k = 7/8 checks.

**Repair:** Methods "about 2.6e-17", Implementation "about 1.3e-17", plus the three source files. **The fix is
free.**

**Minority note (rule 8):** attack L3 independently flagged the same 1.4e-17 slip. Both attackers come from one model
family, so I record this as agreement, not corroboration. The finding stands on my own arithmetic.

## F9 (C20, testing): UPHELD (claim survives; one wording fix)

All tests and mutants are as the attack reports. I did not re-run the mutants and rely on the attacker's run for
that leg.

**Wording leg upheld.** "Block sizes from 1 to 4,096" reads as a range, but the tests use the discrete sets {1, 3, 5,
64, 4096} and {1, 3, 64, (600,) 4096}. Repair: "block sizes 1, 3, 5, 64 and 4,096". **The fix is free.**

**Second leg overturned.** "Every variant is checked" does not assert that the check happened inside the benchmark
job. The AVX-512 path was tested on laila from the same source (`results/kp7_m31_simd_2026-11-10.md:32`). To make
this leg hold, the attacker should have shown that the benchmarked source differed from the tested source.

## F10 (C17, orientation): UPHELD (claim survives; units caveat)

137,216 bytes = 2^8 x 67 x 8, one round, reproduces from `full_kp10_eis_draft.json`. I did not re-run `orient.py`.
The cause does not rest on it alone, though. `kp10_vs_shibuya_2026-11-11.md` records that canonicalizing the
extended-syncmer strings still leaves 33.4 M differences, because syncmer selection and extension depend on
orientation. That supports the orientation cause independently.

The units caveat is right: 48 M counts 51-mer extended syncmers, not k-mers. Repair: "where a strand-specific method
saw 48 million differing extended syncmers". **The fix is free.**

## F11 (bundle of smaller notes): NARROWED

The bundle as a whole carries real but small defects, so the governing verdict is NARROWED.

- **C2: upheld.** The script computes over Z[i]/(q), q = 2^31 - 1, which is F_{q^2}. The caption reuses the letter
  q, which the paper defines as 2^61 - 1. Repair: "F_{p^2}, p = 2^31 - 1" (free). The Eisenstein k = 8 count of
  16,500 for F.R is offered as context. The paper's count is stated for the fourth-root encoding and is correct for
  it, so this is not a defect.
- **C3, C4, C8: upheld as survived.**
- **C13: upheld.** ntHash seed 0 has 8 long gaps (TSV verified). Repair: "all but one method-seed left tens of ..."
  (free). The 0.14 is verified.
- **C14: upheld as survived.** The suggestion to report the direct 0/199,970 SeedNtHash(s) vs SeedNtHash(rc s)
  result is good: it is what tells a reader that 0.818 > 0.747 is not "better". It costs one committed test.
- **C21: narrowed.**
  - "A cs10 chromosome's graph" is the first 20 Mb of cs10. Repair: "13.6 million edges from the first 20 Mb of
    cs10" (free).
  - On "four to six", the attacker's local 6.8x at k = 31 comes from a different host. Laila's 4.2-5.9 across
    k = 31/47/63 supports the printed range, so that sub-point is overturned.
  - 1,370x is a single run on one host, and the attacker's local 409-478x shows the ratio depends on the host.
    Repair: add the qualifier "on one run on the Xeon" (free). Do not adopt "400 to 1,400", which mixes hosts that
    the paper otherwise keeps apart.
- **C5: survived. C6, C7, C19: unchallenged by this lens.** I make no ruling on them here.

## Verdict table

| finding | verdict | why |
|---|---|---|
| F1 | NARROWED | break-even 1/20.3 makes "1/20 and below" too sharp; the fit is validated by the 1/100 row, so the fix is "about 1/20" (free), not "1/100 measured" |
| F2 | NARROWED | "every timing ... counters, one job" is false for the Supplement and the spaced seeds; spaced seeds ARE the same job, and no cross-job comparison is made; 0.76 -> 0.75 |
| F3 | UPHELD | bit-for-bit holds (0 of 736,547,799 mismatches, re-run); the repository has no committed full check |
| F4 | NARROWED | range mixes all-columns and repeat-masked rows and does not generalise; but it is cs10 vs two haplotypes, not only "one plant" |
| F5 | UPHELD | detection is probabilistic (counterexample reproduced); repair: "except with small probability" (free); the attacker's 1/q^c and "free full check makes it deterministic" are overturned |
| F6 | NARROWED | "equal reliability" compares a two-round, true-D-sized sketch with a one-round IBLT on 4 of 64 partitions |
| F7 | NARROWED | canonical ntHash is not GF(2)-linear (0/10000); only its strand hashes are; conclusion stands (free fix) |
| F8 | NARROWED | 2(k-1)/q = 2.6e-17 and the CM31 bound is 1.3e-17, not 1.4; fix across paper, ledger and two results files (free) |
| F9 | UPHELD | tests and mutants reproduce; "1 to 4,096" should list the discrete sizes; the AVX-512-not-in-job leg is overturned |
| F10 | UPHELD | 48 M is an orientation effect (also shown by the 33.4 M after canonicalisation); state the units: extended syncmers |
| F11 | NARROWED | C2 caption field, C13 "tens" (ntHash seed 0 = 8), C21 "chromosome" and the single-run 1,370x need small free fixes; C21 "four to six" stands on laila's data |
