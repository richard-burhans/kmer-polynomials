# Round 2 adjudication of the completeness critic

Adjudicator: skeptic on `completeness-critic.md`. Read in the fixed order: the critic's filing, then the claims (C1-C15)
as quoted in `ADJ-BRIEF-completeness.md` and `ATTACKED-paper.md` (AP) / `ATTACKED-ledger.md` (AL).

**What I checked, and how.** For every gap I first tested whether it is FALSE (someone examined it under another name):
grep of all seven `attack-L*.md` and `adjudication-L*.md` here and every file in `reports/rounds/round-1/` for the
critic's terms plus "encoding", "like for like", "serializ", "padding", "collide", "silent", "undetect", "fwd", "2¹⁶",
"65,536", "m = 64". `scripts/round history "16-byte" --kind adjudication`: a searched absence (round 2 withheld). Code
read: km-peeler at b8ec39b (`bench/data/kmer_poly_work/km-peeler/ibflib.h:73-81`, `ibflib.c:72-100`, `constants.h:20-25`),
`recon/kp10_vs_iblt_job.sh:115-130`, `recon/kp10c.c:1-30`. Data read: `results/kp10_vs_shibuya_2026-11-11.md:120-131,
195-210`, `results/round1_owed_2026-11-11.md:80-115` and its `run_..._k31_fwd_b8_r1.json`. Local work was arithmetic
only (`nice -n 19 python3`). No genome read, no job run.

---

## Gap 1 -- "the 20.5 is km-peeler's storage format, half of it a 64-bit counter" (C1, C2, C6, C13)

The critic: *"a serialized cell is an `int64_t counter` written as 8 bytes plus `keysum[WSIZE]` ... The cell is 16
bytes, and 8 of them are the counter"*, and *"On the GigaDB pair ... 1.82x therefore becomes about 1.0-1.08x ...
C6 reverses."*

**Is it a false gap? No.** The only round-2 mentions are adjudication-L6:90 ("km-peeler has 16-byte cells", stated as
a fact in an Eppstein comparison), attack-L5:134 (the counter as an information source) and attack-L1:269 / adj-L7:42 /
round-1 adj-L5 F5, whose "like for like" rulings are about **one side's message and rounds**, not encoding. C1 names
"entry and cell sizes" as a refutation condition, and no lens ruled on it. The gap is real, and it is not a re-filing.

**Premise, re-run.**
- `ibflib.h:73-81`: `int64_t counter; uint8_t keysum[WSIZE];` and, with `STORE_SEQUENCES 31` (job.sh:52), `GLEN` and
  `RPOS` are defined (`constants.h:20-25`), so no length or position field. `ibf_bucket_store` writes 8 + WSIZE = 8 + 8
  bytes. `iblt_bytes = stat -c %s a.iblt` (job.sh:123). CONFIRMED.
- `results/kp10_vs_shibuya_2026-11-11.md:203-204` already says "Bytes are as the implementations store them ... Each
  could be packed tighter: 61-bit entries, a 62-bit key plus a few counter bits." Nothing like it is in
  `manuscript/paper.md` (grep). CONFIRMED.

**The critic's model is applied to the wrong table.** It sizes the counter from Poisson(3/1.2831 = 2.34). That is the
distribution in the **subtracted** table `d.iblt`. The bytes counted are `a.iblt`, **one party's own table**, which
`ibltseq build` fills with all of A (job.sh:119). Its counters average 3 × 442,159,825 / (1.2831 × 3,238,101) = **319**,
entropy about 6.2 bits, maximum above 255. So:
- "a 1-byte signed counter holds it" is false of the shipped message as argued;
- "Entropy-coded, the table costs 10.36" is false: it prices the difference table's counters, which nobody sends.

**The conclusion survives anyway, by a route the critic did not give.** The receiver needs only the difference's
counter, so the counter can be sent modulo 2^b. Over the ~4.15 M cells, Poisson(2.34) predicts 0.79 cells with
|count| ≥ 14 and 0.12 with ≥ 15, so b = 5 to 6 bits suffices (a prediction, not measured). With a 62-bit key:
- GigaDB: 1.2831 × (62 + 5..6)/8 = **10.75-10.91 B/diff**, against the sketch's 11.27 × 61/64 = **10.74**. The 1.8x
  gap does vanish (ratio 1.00-1.02).
- D/8 (C2): 38.5 × 67/128 ≈ **20.2**, against 20.8 × 61/64 = **19.8**. Claim 5's gap vanishes too.
- Draft (C6): 162,938/16/7,931 = 1.284 cells/diff, so ≈ **10.8**, against 17.3 × 61/64 = **16.5**. Under packing the
  ordering does reverse.
- Every one of these is arithmetic about encodings **neither of which was built or run**.

**Diagnosis: UPHELD.** C1's, C2's and C6's comparisons are of shipped file formats, half of km-peeler's being a 64-bit
counter, and the paper does not say so. Abstract claim 4 and claim 5 read as comparisons of data structures, which the
measurement does not support.

**Repair (a), "give packed figures for both methods beside them": NARROWED.** Printing packed byte counts (11.55,
10.36, "about 11.6", "about 21.7") would install counterfactual numbers for a table that was never built, one of them
(10.36) wrong. **Endorsed instead, free** (Results, after Table 4, and the same qualifier in the table caption):

> "Bytes are as each implementation stores them: 8-byte field elements for the sketch, and for km-peeler 16-byte
> cells, an 8-byte key sum and an 8-byte counter. Packed, 61-bit entries and a 62-bit key with a few counter bits would
> make the two about the same size on the GigaDB pair, and the table the smaller on the draft pair; we built neither
> packing."

and in the abstract, claim 4: "... against 20.5 for a one-round invertible Bloom lookup table **as km-peeler stores
it, half of each cell a 64-bit counter**". Whether claims 4 and 5 stay in the abstract at all, once the gap is
attributed to encoding, is the user's decision. The free sentence makes the abstract true; it does not keep it
interesting. The Discussion's "the sketch spent 11.3 bytes per difference ... where a one-round table spent 20.5" and
C13's "advantage ... due to resending only overfull buckets" need the same qualifier (free).

**Repair (b), "keep the `a.iblt` files ... measure the maximum |counter| and the counter entropy": NARROWED.** Measuring
`a.iblt` measures the 319-mean counters, the wrong instrument. If any packed number is ever printed, what is **owed**
(one laila run, the job deletes `$d`) is the maximum |counter| of **`d.iblt`** across the 1.05·D trials, which is what a
modular counter must hold. Not owed for the endorsed wording.

**Correct minority.** Only the critic saw this. No attacker in seven lenses did, though C1 named the condition.

## Gap 2 -- "a hash collision can drop a difference without any check firing" (C5)

The critic: *"If a differing class collides with any class in A ∪ B, its value cancels or merges, and no c-check or
flag count can see it ... bounded by D·|A ∪ B|·2(k−1)/q = 0.037 per run."*

**False gap? No.** Round 2's collision hits (attack-L1:249, adj-L1:214) rule on AL:159's row count; L3, L5, L7 confirm
"neither method returned a wrong difference" from `undetected_failures 0`. Nobody asked whether a collision is a silent
failure mode. Round 1 adjudication-L5 F4 (its "genome-scale reading", l.83-87) ruled on the per-assembly 2.5-pair allowance of the hash, not on the sketch.

**Premise, re-run.** `kp10c.c:8-11` deduplicates the distinct values. The arithmetic is right: 3.24 × 10^6 × 4.42 ×
10^8 × 60 / (2^61 − 1) = 0.037; ideal 61-bit hash 6 × 10^-4. The observed collision (results :127-131, 1 of 37 runs,
1 of 31 OS-seed draws per adj-L1:214) was between shared classes.

**One overstatement.** "Any class in A ∪ B" is too wide. A differing class colliding with another class on its **own**
side is not lost: the merged value is still a difference, and the owner's rescan emits every k-mer with that value. The
silent loss needs a collision with a class the other party also holds (a shared class or one of the other side's
differences). The bound is unchanged as an upper bound.

**Is it a defect?** NARROWED. Nothing published is false: the recovery was checked against KMC. But AP:146's "the
sketch detects an overfull bucket only with high probability" offers overflow as the sketch's probabilistic exposure,
and the collision exposure, which occurred in this job, is omitted from the paper (it is in the ledger only).

**Repair: UPHELD with a tightened wording, free:**

> "A differing k-mer whose hash equals that of a k-mer the other party holds would cancel without any check firing;
> the bound puts this below 0.04 per run on the GigaDB pair, and the one collision seen, in 1 of 31 seed draws, was
> between shared k-mers."

## Gap 3 -- Figure 4A caption "fails its check" deterministically (C9 scope)

The critic: *AP:148 "A bucket with more than m differences fails its check and is split" ... "the power sums of the
difference", where Methods correctly says signed.*

**False gap? No.** C9 is 4B only, and the sole mention of 4A (attack-L1:90) is about its mechanism. Round 1 already
**ruled** that detection is probabilistic (adj-L1 F4 UPHELD; adj-L7 F5 endorsed "except with small probability"), and
Methods (AP:63) carries the hedge. The caption contradicts a settled ruling, so this is a missed application of it, not
a new dispute.

**Verdict: UPHELD, repair endorsed, free:** "A bucket with more than $m$ differences fails its check, except with small
probability, and is split ..." and "leaves the signed power sums of the difference". A schematic caption is a weak
venue for the error, hence WEAK, but the ruling it breaks is binding.

## Gap 4 -- the forward-only control is "the numpy prototype, not 'the same sketch'" (C3)

The critic: *"The control is `run_..._k31_fwd_b8_r1.json` ... numpy, 2^8 buckets ... AP:150 nonetheless says 'in the
same sketch'."*

**False gap? Partly.** L5 F6.2 attacked this exact sentence and adj-L5:177-180 endorsed "faced". The *provenance* leg
was not raised by L4, so that leg is new.

**Premise, re-run.** The fwd JSON shows `recon/aswritten.py` (numpy) with hashes from `hash/kp7e_stream.c`: bits 8,
m 64, c 2, bucket key "top bits of the value", 137,216 bytes (`results/round1_owed_2026-11-11.md:82-94`). That is the
**same sketch design and size** as the published as-written S run. It is a different implementation from `kp10c`.
"Same sketch" is therefore true of the design and misleading only about the code. The critic's "superseded" is also
wrong: the paper still uses the numpy prototype for k = 63 (AP:63).

**Diagnosis: NARROWED.** A small provenance imprecision, not a misattribution of the result. **Repair, free, merged
with adj-L5's endorsed verb so the two do not collide:**

> "a forward-only polynomial hash, in the numpy prototype of the same sketch (2^8 buckets, m = 64, c = 2), faced
> 398,217,248 differences and overflowed every bucket."

AL:142 takes the same words. **Owed item (a kp10c forward-only run or KMC `-b`): OVERTURNED.** The control's only job
is to exceed the sketch's capacity of 2^8 × 64 = 16,384 differences, which a count wrong by four orders of magnitude
would still do, and the critic's own orientation check (0.450 against 45.4 % flipped bases) finds the count
plausible. A laila run to confirm a number that cannot change the conclusion is not owed.

## Gap 5 -- "2^16" and "m = 64, c = 2" have no ledger rows (C15)

The critic: *"AL for '2^16', '65,536', 'm = 64' and 'c = 2'. There are zero rows."*

**Premise, re-run.** "2^16" / "2¹⁶" / "65,536": zero rows. CONFIRMED. But **AL:113 "m + c = 66 | field
multiplications per element"** exists. The critic's search terms missed it. That row cites the superseded
`kp10_vs_shibuya_2026-11-10_prior.md:8` for parameters that the current C run states at
`results/kp10_vs_shibuya_2026-11-11.md:64` ("m = 64, c = 2"). The identity reproduces: (65,536 + 2,526) × 67 × 8 =
36,481,232 B, AL:103's C-run size.

**Verdict: NARROWED.** 2^16 lacks a row. m and c have a row, but by their sum and from the prior file. **Repair, free:**
add "| 2^16 | KP-10 first-round buckets, GigaDB, C run ((65,536 + 2,526 halves) × 67 × 8 = 36,481,232 B) |
results/kp10_vs_shibuya_2026-11-11.md:104 |", and repoint AL:113 to `results/kp10_vs_shibuya_2026-11-11.md:64`, with
m = 64, c = 2 in its meaning.

## §3 -- "small probability" is unquantified (framing)

The critic: *"no lens could attack a probability the paper never states ... Owed: a stated bound, or a measured
false-accept rate."*

**This is a settled point.** Round 1 adj-L1 F4 noted "'negligible' is unquantified", and adj-L7 F5 endorsed the
unquantified hedge "except with small probability" while **overturning** the attacker's attempted 1/q^c bound. The
paper carries that endorsed hedge. Re-filing it without new evidence is a weaker finding. The critic's own framing,
"a finding about the claim list, not about the attackers", concedes that no published claim is false. **OVERTURNED as
a gap.** The false-accept measurement is a nice-to-have, not owed. Should the paper ever quantify the hedge, a
measurement on synthetic overfull buckets through the C decoder (with its flag-count checks, `kp10c.c:294-296`) is the
right instrument, not a bound derived from the numpy decoder.

---

## What the critic should have checked (rule 7)

- Gap 1: which file `iblt_bytes` measures (`a.iblt`, one party's full table), before modelling its counters as the
  difference table's.
- Gap 2: whether a same-side collision is in fact lost, given that the owner rescans its genome by value.
- Gap 4: the fwd JSON's parameters against the S run's, before calling it a different sketch, and AP:63 before calling
  the prototype superseded.
- Gap 5: AL for "m + c", not only "m = 64".
- §3: round 1's adj-L7 F5, which already ruled on this hedge.

## Free repairs endorsed

1. Gap 1: the Results/caption sentence on stored encodings, the abstract qualifier "as km-peeler stores it, half of
   each cell a 64-bit counter", and the same qualifier in the Discussion. No packed numbers.
2. Gap 2: the collision sentence as worded above.
3. Gap 3: "except with small probability" and "signed" in the Figure 4A caption.
4. Gap 4: "in the numpy prototype of the same sketch (2^8 buckets, m = 64, c = 2), faced ..." (paper and AL:142).
5. Gap 5: a 2^16 row, and AL:113 repointed to the current results file.

**Owed only if packed figures are ever published:** max |counter| of `d.iblt` over the 1.05·D trials (one laila run).

| finding | verdict | why |
|---|---|---|
| Gap 1 | UPHELD | 16-byte cells with an 8-byte counter, undisclosed in the paper; packed, the GigaDB and D/8 gaps vanish. The critic's counter model fits the wrong table (shipped `a.iblt` counters average 319) and its 10.36 is wrong; packed figures OVERTURNED as a repair; the stored-encoding disclosure and abstract qualifier are endorsed, free |
| Gap 2 | NARROWED | the silent-cancel mode is real and was observed (harmless); it needs a collision with a class the other party holds, not with "any class in A ∪ B"; the tightened sentence is free |
| Gap 3 | UPHELD | the caption contradicts round 1's binding hedge; "except with small probability" and "signed" are free |
| Gap 4 | NARROWED | same sketch design and size, a different (numpy) implementation, not superseded; rewording merged with adj-L5's "faced" is free; the owed rerun is OVERTURNED |
| Gap 5 | NARROWED | 2^16 has no row, but m + c = 66 does (AL:113, from the prior file); add the row and repoint, free |
| §3 framing | OVERTURNED | settled in round 1 (adj-L1 F4, adj-L7 F5 endorsed the unquantified hedge); no published claim is false; the measurement is optional, not owed |
