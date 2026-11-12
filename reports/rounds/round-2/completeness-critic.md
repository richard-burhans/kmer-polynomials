# Round 2 completeness critic

Critic: completeness. Frozen text: `ATTACKED-paper.md` (AP) and `ATTACKED-ledger.md` (AL) in this directory.

**What I read before asserting any absence.** All seven `attack-L*.md` and all seven `adjudication-L*.md` in this
directory, verdict tables in full and finding bodies where a term hit. I grepped all 14 files for each gap's subject
under other names; the terms are listed in each gap. `scripts/round history` (closed rounds only; round 2 is withheld)
was run for "int64 counter", "16-byte cells", "packed tighter" and "collision". My local checks were arithmetic only
(`python3`, under `nice -n 19`), plus reads of `recon/kp10c.c`, `recon/kp10_vs_iblt_job.sh`, the pinned km-peeler
source (`bench/data/kmer_poly_work/km-peeler`, the clone L1 named), `results/kp10_vs_shibuya_2026-11-11.md` and
`results/round1_owed_2026-11-11/draft_orientation.json`. No genome was read, and no job was submitted.

## 1. What filed, per lens (UNCHALLENGED vs UNRUN)

| lens | findings | claims left without a ruling of substance |
|---|---|---|
| L1 correctness | 9 | none (C14 timings deferred "to another lens") |
| L2 sampling | 7 | C11, C12 "not reached under this lens" (UNRUN for this lens); C14 UNCHALLENGED with a reason |
| L3 statistics | 5 | C8, C10, C11, C12, C14 are blank rows, with no reason given, so I cannot tell them from UNRUN |
| L4 provenance | 9 | C13 (declared out of lens) |
| L5 inference | 12 | none |
| L6 prior-art | 6 | C5, C9, C10, C11, C14, C15 "not reached"; C15 had a partial row check |
| L7 reproduction | 12 | none |

Sixty findings filed. Every lens filed and was adjudicated, so no lens is UNRUN. **L3 and L6 are thin.** L3 has five
blank rows. L6 left six claims unreached, which is expected for a literature lens. Every claim still has at least four
substantive rulings.

**Effort converged on a few defects.** The x/λ sentence was filed by 4 lenses (L1 F5, L4 F1, L5 F9, L7 F2). The
km-peeler consecutive seeds were filed by 5 (L1 F4, L2 F3, L3 F4, L4 F3, L7 F6), and the k = 63 provenance by 5 (L1 F6,
L2 F5, L3 F5, L4 F2, L5 F6.3, L7 F3). The power-of-two account of "comes from the second round" was filed by 5 (L1 F3,
L2 F2, L3 F2, L5 F1, L6 F2). Meanwhile one refutation condition that C1 names in so many words, **"entry and cell
sizes"**, drew no finding at all (Gap 1). That is the framing miss. C1 bundles five refutation conditions, and the
lenses took the ones that have a code or data trail to follow.

## 2. Gaps

### Gap 1: the 20.5 is km-peeler's storage format, half of it a 64-bit counter. Packed, the 1.8x gap disappears and the draft-pair comparison reverses (C1, C6, C13, abstract claim 4). STRONG

**Searched:** "counter", "int64", "keysum", "16-byte", "16 bytes", "cell size", "bytes per cell", "compress", "gzip",
"entropy" and "packed" across all 14 files, plus `history` for "int64 counter", "16-byte cells" and "packed tighter".
- L5 F4 (attack l.134) mentions the `int64 counter` only as a source of a size estimate.
- adjudication-L6 l.90 states "km-peeler has 16-byte cells" as a fact, while comparing with Eppstein.
- No file asks what those 16 bytes are, or whether the comparison is like for like in encoding.
- Closed rounds have no mention.

**The artifact.**
- `km-peeler/ibflib.h:71-79` and `ibflib.c:81-83` (b8ec39b): a serialized cell is an `int64_t counter` written as 8
  bytes plus `keysum[WSIZE]`. With `STORE_SEQUENCES 31` (`kp10_vs_iblt_job.sh:52`), `constants.h:23-25` gives
  WSIZE = 8 and sets GLEN and RPOS, so there is no length or position field. The cell is 16 bytes, and 8 of them are
  the counter.
- `iblt_bytes` is `stat -c %s a.iblt` (`kp10_vs_iblt_job.sh:140`), so the 20.5 is the size of the shipped file.
- The results file already concedes the point (`results/kp10_vs_shibuya_2026-11-11.md:203-204`): "Bytes are as the
  implementations store them ... Each could be packed tighter: 61-bit entries, a 62-bit key plus a few counter bits."
- `grep 'packed tighter|as the implementations store'` finds nothing in `manuscript/paper.md`. The caveat never reached
  the paper, and nobody priced it.

**Why it matters (arithmetic, not a run).**
- At 1.2831 cells per difference (L6 F3, adjudicated), each cell holds on average 3/1.2831 = 2.34 keys.
- A Poisson(2.34) counter has an entropy of 2.59 bits, and a 1-byte signed counter holds it.
- With an 8-byte keysum and a 1-byte counter, the table costs 1.2831 × 9 = **11.55 B/diff**, against the sketch's
  **11.27**, or **10.74** with its 61-bit entries packed. Entropy-coded, the table costs 10.36.
- On the GigaDB pair, the abstract's 11.3 against 20.5 (1.82x) therefore becomes about 1.0-1.08x.
- On the draft pair (C6), the table has the same cells per difference (162,938 / 16 / 7,931 = 1.284). Packed, it costs
  about 11.6 against the sketch's 17.3 in one round. **C6 reverses.** L1, L3, L4, L6 and L7 all ruled C6 SURVIVED or
  NARROWED on its numbers alone.
- The D/8 row (20.8 against 38.5) scales the same way: the table's side roughly halves, to about 21.7.
- The sketch's own per-bucket 8-byte count is 1 entry in 67, so it barely moves.

**What would close it.**
- (a) Free: state that the bytes are the shipped encodings, and give packed figures for both methods beside them. The
  results file has the sentence already.
- (b) Owed, one laila run: keep the `a.iblt` files from the 1.05·D trials and measure the maximum |counter| and the
  counter entropy. That replaces my Poisson model.
- Until (a) is done, abstract claim 4's "against 20.5 for a one-round invertible Bloom lookup table" compares one
  implementation's padding, not the two data structures.

### Gap 2: a hash collision can drop a difference without any check firing, and the run saw a collision that the paper omits (C5). MODERATE

**Searched:** "collision", "birthday", "silent", "undetect" and "false accept". The hits are L1 F9 (the AL:159 row-count
quibble), L6 F5 (the bound as a selling point), adjudication-L1 l.214 and adjudication-L2 l.112 (km-peeler partition
collision patterns). None examined collisions as a failure mode of the sketch. Round 1's L3 F2 (NARROWED, "about 2.5
collisions") ruled on the hash's claims, not on the set sketch, so this is a new angle and not a re-filing.

**The artifact.**
- AP:146: "Neither method returned a wrong difference, but the sketch detects an overfull bucket only with high
  probability". That sentence names overflow as the only probabilistic failure.
- `kp10c.c:10-11` deduplicates S values. If a differing class collides with any class in A ∪ B, its value cancels or
  merges, and no c-check or flag count can see it.
- That event is bounded by D·|A ∪ B|·2(k−1)/q = **0.037 per run** at the bound, or 6 × 10^-4 for a random function.
- This very job saw a collision: `results/kp10_vs_shibuya_2026-11-11.md:127-131` and `tables/collisions.tsv:2` record
  1 run in 37, which is 31 distinct draws (adjudication-L1 l.214). It was harmless only because both classes were
  shared. It is in the ledger (AL:159) and absent from the paper.

**Why it matters.** The numbers stand, so this is not a refutation. But "neither returned a wrong difference, but ..."
presents the sketch's only exposure as overflow, when a second, undetectable exposure exists and occurred here.

**What would close it.** Free: one sentence, "a collision of S with a differing class would drop it undetected
(≤ 0.04 per run by the bound); one of 31 draws had a collision, between shared classes."

### Gap 3: Figure 4A was never in any claim, and its caption drops the hedge (C9 scope). WEAK

**Searched:** "4A" and "Figure 4A". The only hit is L1 l.90 ("Figure 4A's ... may stay as the mechanism"). C9 covers 4B
only, and the ledger's 4A row (AL:213) cites the paper's own Methods.

**The artifact.**
- AP:148 (A): "A bucket with more than m differences fails its check and is split". That is a deterministic statement.
- Methods (AP:63) and Results (AP:146) both say detection fails "with small probability". Round 1's adjudication-L1 F4
  UPHELD a 140-difference bucket that passes, and results :126 says that "stands".
- 4A also says "the power sums of the difference", where Methods correctly says *signed* power sums.

**Why it matters.** The caption restates as a certainty the one property round 1 made the paper hedge.

**What would close it.** Free: "fails its check, except with small probability", and "signed power sums".

### Gap 4: the forward-only control comes from the superseded numpy prototype, not "the same sketch", and nothing independent confirms its count (C3). WEAK

**Searched:** "forward-only", "forward only", "398" and "not rerun". L5 F6.2 attacked the verb ("saw"). Adjudication-L5
l.177-180 UPHELD that and endorsed "faced". L3 l.179 and adjudication-L3 l.192 only quote the number. **L4, the
provenance lens, never mentions the control.**

**The artifact.**
- The control is `run_..._k31_fwd_b8_r1.json` (L5 F6.1): numpy, 2^8 buckets.
- The ledger's carried-items row (AL:194) says "The forward-only control was not rerun", while the 7,931 beside it
  comes from the C run.
- AP:150 nonetheless says "in the same sketch". AL:142 says so too.
- 398,217,248 is `D_true_values` from that prototype. Unlike every other difference in the section, it has no KMC
  confirmation: KMC checked canonical sets only.

**My check.**
- 398,217,248 / (2 × 442,159,825) = 0.450.
- `draft_orientation.json` gives 45.36% of bases flipped, plus 5.5% on mixed contigs.
- The count is therefore plausible. The defect is the attribution, not the number. "The draft writes 45 percent ..."
  (AP:150) is otherwise unattacked, and I found it consistent. I do not file it.

**What would close it.**
- Free: "in the numpy prototype of the sketch".
- Owed: a forward-only run in `kp10c`, or a KMC non-canonical (`-b`) difference, on laila.

### Gap 5: two published sketch numbers have no ledger row, contrary to L1 F9's "every number has a row" (C15). WEAK

**Searched:** AL for "2^16", "65,536", "m = 64" and "c = 2". There are zero rows.
- AP:144 publishes "started from 2^16 buckets".
- AP:63 publishes m = 64, c = 2.
- L1 F9 asserted "Every number in AP:117-136 has a row". Adjudication-L1 F9 did not test that sentence.

**My check.** (65,536 + 2,526) × 67 × 8 = 36,481,232 B, exactly AL:103's C-run size. So 2^16, m + c = 66, the 8-byte
entries and L1's 2,526 resent halves all reproduce. The gap is the missing rows, not the values.

**What would close it.** Free: add rows citing `results/kp10_vs_shibuya_2026-11-11.md` (section 2), and this identity.

## 3. A claim no evidence could settle

"Except with small probability" (Methods) and "only with high probability over the random hash" (C5) carry no number.
Round 1 (adjudication-L1 F4) already noted that "'negligible' is unquantified". The C decoder adds checks the numpy one
lacked (own-root flag counts, `kp10c.c:294-296`), and no lens could attack a probability the paper never states. This
is a finding about the claim list, not about the attackers.

**Owed:** a stated bound, or a measured false-accept rate on synthetic overfull buckets. The latter is a small local or
laila job.

## 4. Summary list

| # | gap | strength | why it matters | what closes it |
|---|---|---|---|---|
| 1 | IBLT bytes are a 16-byte cell with an 8-byte counter | STRONG | packed, 11.3 vs 20.5 becomes about 1.0-1.08x, and C6 (17.3 vs 20.5) reverses to 17.3 vs about 11.6 | free: report packed bytes for both; owed: measure the counters from kept `a.iblt` files |
| 2 | S collisions give a silent-miss mode, observed once | MODERATE | C5's sentence names overflow as the only probabilistic failure | free: one sentence with the 0.037 bound and the 1-of-31 observation |
| 3 | Figure 4A caption says "fails its check" deterministically | WEAK | contradicts the hedge round 1 imposed | free: add the hedge and "signed" |
| 4 | the forward-only control is the numpy run, called "the same sketch" | WEAK | provenance misattributed; no independent count | free: rewording; owed: kp10c forward-only run or KMC `-b` |
| 5 | "2^16", "m = 64, c = 2" have no ledger rows | WEAK | C15's own test | free: add the rows (values verified) |
| — | "small probability" is unquantified | framing | no lens could test it | owed: a bound or a false-accept measurement |
