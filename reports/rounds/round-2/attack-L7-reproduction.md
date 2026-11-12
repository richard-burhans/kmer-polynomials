# Round 2 attack, lens L7: reproduction

Attacker: reproduction lens. Frozen text: `reports/rounds/round-2/ATTACKED-paper.md`, `ATTACKED-ledger.md`.

**What I re-ran.** I re-ran `recon/kp10_vs_iblt_summarize.py` on the committed `run.log` into my scratch
directory under `nice -n 19 ionice -c3`. All eight tables (`collisions`, `iblt_protocol`, `iblt_trials`, `kp10_dest`,
`kp10_seeds`, `sets`, `timed`, `truth`) came out **byte-identical** to the committed `tables/*.tsv`, so the tables are a
faithful projection of `run.log`. I recomputed the derived numbers from the raw JSON lines with python. I read
`recon/kp10c.c`, `recon/aswritten.py`, `recon/recon61.py`, `hash/kp7e_stream.c` and `recon/kp10_vs_iblt_job.sh`. Nothing
heavy was run, and no job was submitted.

**Method cards.** None of the claims here is a regression, ICC or brand analysis. The project's card on power (cited in
round 1's L5 adjudication as `power.md`) governs C5, and C5 complies with it (F9).

**Literature.** The project's literature is held in `../library/papers_md/`. I read
`shibuya_2022_genomic-reconciliation.md` for C8 and C13.

**Prior rounds.** I searched with `scripts/round history` for "one side", "both sides", "oracle", "message", "own
differing" and "listing gives". The one prior ruling on my main point is round 1's L5 adjudication F5. I disclose it in
F1, together with what has changed since.

---

## F1 -- "recovered every difference with 11.3 bytes" counts one of the two messages that recovering every difference took

**Claims:** C1 (REFUTED as worded), C2 and C6 (NARROWED).

**The defect.** The 11.27 B/diff is `sketch_bytes / D`, and `sketch_bytes` is ONE party's sketch:
`recon/kp10c.c:14-15` says "ONE side's message". The run, however, recovered both sides' k-mers
(`run.log:28`: `kmers_recovered_A_minus_B 3237326`, `kmers_recovered_B_minus_A 775`). In the code that produced the
number, that needs both sketches to cross:

- **Each side needs the other's sketch.** `decode()` starts from `S = SA - SB` (`kp10c.c:259`), so each side must
  hold the other's sketch before it can decode at all.
- **Acceptance needs both sides' checks.** A bucket is accepted only when side A's check AND side B's check pass
  (`kp10c.c:21-24`, `:300-304`).
- **Each party gets only its own k-mers.** The paper says so itself: "Each party checks its buckets with its own data
  only and finds its own differing k-mers", and "Each party of the sketch recovers only its own differing k-mers,
  while one listing of the table gives both sides'".

So "every difference recovered" took A's sketch to B **and** B's sketch to A. km-peeler's 20.53 is one file
(`kp10_vs_iblt_job.sh:119-126`: `iblt_bytes` = the size of `a.iblt`), and one `ibltseq list` returned both sides
(`run.log:30`: listed 3,237,326 + 775). Counted like for like, as the total bytes that recovered every difference,
using the same raw outputs:

| quantity | sketch, both messages | table, one message | source |
|---|---|---|---|
| GigaDB at the true D | 2 × 36,481,232 / 3,238,101 = **22.53** (30-seed range 22.49–22.57) | 20.53 | `kp10_dest.tsv:6`, `kp10_seeds.tsv:2`, `iblt_trials.tsv:9` |
| GigaDB at D/8 (C2) | 2 × 67,212,256 / 3,238,101 = **41.51** | 38.49 | `kp10_dest.tsv:2`, `iblt_protocol.tsv:6` |
| draft at the true D (C6) | 2 × 137,216 / 7,931 = **34.60** | 20.54 | `kp10_seeds.tsv:3`, `iblt_trials.tsv:3` |

Under this accounting, the sketch's headline advantage reverses on every row.

**Premise test.** Could one message suffice? With B's sketch alone, A decodes P (its own) and Q, but it gets B's
differences only as roots of the unfactored Q. Those are hash values, not k-mers, and root-finding is "not
implemented" (`results/kp10_vs_shibuya_2026-11-11.md` §6). A one-message protocol plus root-finding plus a reply of
775 values is plausible. It was not built or run, and acceptance in the code that ran requires both sides' checks.
In the draft pair all 7,931 differences are on B's side, so a one-way A→B message might suffice there. That variant
was not run either.

**What changed since round 1.** Round 1's adjudicator ruled that "both byte counts are one message's worth … the
accounting is like for like" and that "the sketch's byte count does not depend on who decodes" (adjudication-L5.md
F5, NARROWED). At that time the decoder was the one-process `recon61.py`, and the adjudicator OVERTURNED a sentence
about a two-party decoder because none had been built. One has now been built: `kp10c` is a two-party decoder. It
shows that each party decodes only from the difference of the two sketches and returns only its own k-mers. The
paper now says this in Methods and Results, yet it still prices "every difference" at one message.

**Correction.**
- **Abstract.** Write "…each party recovered its own differing k-mers in 30 of 30 runs from 11.3 bytes per difference
  sent in each direction over two rounds (22.5 in total), against 20.5 for a one-round invertible Bloom lookup table
  from which one party lists both sides", or restrict the comparison explicitly to a one-way message.
- **C2 and C6.** Give the same both-ways totals (41.5 vs 38.5; 34.6 vs 20.5), or say that the bytes are per direction.
- **Check owed (cheap, local).** Run kp10c's acceptance with side A's check alone, on a tiny pair via
  `kp10_vs_iblt_launch.sh dry`, before claiming any one-way figure.

## F2 -- "the GigaDB sketch ran before x and λ were separated" is false of the run the paper's GigaDB numbers come from

**Claim:** C11 (REFUTED, second sentence).

**Evidence.**
- **The run used separate seeds.** `run.log:28` (the 11.27 run) records `"seed_x": "0x5eed", "seed_lam": "0x1a3d"`.
  These are two separate seeds, set by `kp10_vs_iblt_job.sh:25-26`. The 30 seed runs behind "30 of 30" and "11.3"
  used `"seed_x": "os"` with λ also from the OS (60 records in `run.log`).
- **The separation predates the code.** Separation is commit `8c34af9` ("Draw x and lam independently", 2026-11-11
  01:14). The run's code is `e76fad1` (13:47 the same day). `git merge-base --is-ancestor 8c34af9 e76fad1` succeeds.
- **The sentence is stale.** It describes the superseded 2026-11-10 numpy run, which the paper no longer reports.

The first sentence of C11 holds: `kp10c` hashes as written, via `kp7e_sym_seq` (`kp10c.c:149-160`), and keys buckets on
the top bits of S (`kp10c.c:171-176`). The ledger has no row for the x/λ sentence.

**Correction.** Delete ", and the GigaDB sketch ran before $x$ and $\lambda$ were separated". If a caveat is wanted,
write: "The GigaDB and draft comparisons drew x and λ independently, from fixed seeds for the timed and sweep runs and
from the operating system for the 30 repeated runs."

## F3 -- the k = 63 provenance sentence describes the wrong code and the wrong sizing rule

**Claim:** C3 (NARROWED).

The paper says: "The $k = 63$ run used a numpy prototype, whose hashes match the C code on 20,000 sampled k-mers, with
buckets sized from the true difference."

- **The hashes are C, not numpy.** The k = 63 values come from `hash/kp7e_stream.c`, which is C (`kp7e_sym_seq`).
  `recon/aswritten.py:12,47` streams it. numpy does only the reconciliation (`recon61.py`). The k = 63 truth check
  hashes the KMC strings with the C direct evaluator `kp7e_direct` (`aswritten.py:222-230`).
- **The 20,000-k-mer check is about something else.** It comes from `kp10_vs_shibuya_2026-11-10_prior.md:32-33`: the
  numpy ports `recon/kp7e_hash.py` of the k = 31 packed-code pipeline. Those ports cannot evaluate a k = 63 k-mer, so
  that check says nothing about the k = 63 run. Ledger row 115 carries the same mis-scoping.
- **The sizing rule differs from the one Methods states.** The run used `bits_for`, D/2^b ≤ **m/2** = 32
  (`aswritten.py:136-142`; run JSON `"bits": 10, "bits_rule": "auto: bits_for(exact value D, m)"`). Methods states
  "at most 0.8m". The 0.8m rule would give 2^9 buckets (21,382/512 = 41.8 ≤ 51.2) and 274,432 B, half the 548,864 B
  reported. Whether 2^9 finishes in one round was not run. So "549 KB and one round" was bought at a more generous
  load than the rule the paper states.
- **The decoder was the one-process check.** `recon61.py:122-128` checks `psum(fa) - psum(fb)` with both sides' data,
  which is the decoder round 1 flagged (adjudication-L5 F5). The paper's "Each party checks its buckets with its own
  data only" therefore does not cover the k = 31 as-written run of `round1_owed` or the k = 63 run. The k = 31 result
  is independently reproduced by kp10c's two-party decode (137,216 B, 1 round, 30/30 OS seeds; `kp10_seeds.tsv:3`).
  k = 63 has no such reproduction.

**Correction.** "At k = 63 the hashes came from the C roller and the reconciliation from a numpy prototype that
checked both sides' data in one process, with $2^{10}$ buckets from the true difference (at most m/2 per bucket)."
Drop the 20,000-k-mer clause. The owed check is a laila job, and I did not submit it: kp10c generalised to k = 63
(it rejects k > 32, `kp10c.c:515`) or aswritten.py run at 2^9, to test the stated rule.

## F4 -- "resent about 1,250 of them" halves what round 2 sent

**Claims:** C4 and C15 (NARROWED).

Table 4's caption says "The sketch started from $2^{16}$ buckets and resent about 1,250 of them in its second round."

- **What round 2 sent.** It sent each failed bucket's two halves, each a full 67-entry sketch:
  `run.log:27`, "round 2: bits 17, **2526** buckets sent" for 1,263 failures. Over 30 seeds, the mean is 1,254
  failures and about 2,508 buckets sent.
- **The bytes are right.** 2,526 × 67 × 8 = 1,353,936 B, which equals 36,481,232 − 35,127,296. Only the sentence is
  wrong.
- **The ledger misreads its source.** Ledger row 178 ("about 1,250 | buckets resent in round 2") cites
  `kp10_vs_shibuya_2026-11-11.md:118`, which says the buckets "were split in two and resent".

**On "comes from the second round" (C4).** It reproduces as arithmetic: round 1 alone is 10.85 B/diff
(`kp10_seeds.tsv:2`), and the one-round run is 21.70 at 2^17 (`kp10_dest.tsv:7`). But 2^17 is simply the next size
available. Buckets are top-bit ranges, so the code can only double. The 21.70 is therefore the cost of a power-of-two
grid, not a measured minimum one-round size. One fixed-seed run, too, as the paper says.

**Correction.** "…resent the two halves of about 1,250 failed buckets (about 2,500 buckets) in its second round."
Change ledger row 178 to match.

## F5 -- the extended-syncmer row is from a different, superseded run, and the "7 missed" has a partial trace the paper omits

**Claims:** C8 and C7 (NARROWED).

- **A different run.** Table 4 sets the extended-syncmer row (4.2, "no, 1,986,807 spurious and 7 missed", 1 round)
  beside rows from laila job b69c5c76. Its numbers come from the 2026-11-10 run
  (`kp10_vs_shibuya_2026-11-10_prior.md:99-112`) with **5 trials**, not 30, in a different job. The results file says
  this method "was not rerun" (§6). The paper says neither, and its Methods sentence "each method ran 30 times with
  independent random seeds on one host" reads as covering it. Ledger rows 109-112 do cite the prior file.
- **The 61% and the authors' 3% compare like quantities.** Their Table 1 reports err/diff, spurious over the true
  difference: 11.03/325.29 = 3.4%, 60.03/1708.78 = 3.5% and 357.57/15342.81 = 2.3%
  (`../library/papers_md/shibuya_2022_genomic-reconciliation.md:747-782`, "within about 3% of the true set
  difference size"). 1,986,807 / 3,238,101 = 61.4%. That part SURVIVES. The ledger's "38% of the returned set" is a
  different quantity, and the paper rightly does not compare it.
- **"Missed 7, which we did not trace" omits a partial trace.** The prior file says "the listing returned 507,219 of
  507,221 every trial" (`prior.md:104`). Two extended syncmers were never listed, identically in all 5 trials. A
  deterministic shortfall in a randomised IBLT points to the harness or input, for example the D_ext count or
  duplicate or N-containing syncmers, rather than to a peeling failure. It is also an IBLT listing that exited
  incomplete, which bears on "Neither method returned a wrong difference" if km-peeler exited 0. The exit status is
  not recorded in the prior file. Round 1's critic gap 5 (harness error not ruled out) therefore stays open.

**Correction.** In the Table 4 caption: "The extended-syncmer row is from an earlier run (5 trials, 2026-11-10); its
listing missed 2 of 507,221 differing extended syncmers in every trial." Owed: rerun Section 3.5 in the b69c5c76
harness. That is a laila job; I did not submit it.

## F6 -- km-peeler's seeds were fixed consecutive integers, not "independent random seeds"

**Claim:** C7 (NARROWED).

- **What the trials used.** Every `iblt_trial` in `run.log` used seeds **1001–1030**. The same 30 seeds served every
  ratio and both pairs; I tallied this from the 300 trial records. The protocol repetitions used other fixed seeds
  (e.g. `"seed": 111000`, `run.log`, iblt_protocol line).
- **What KP-10 used.** KP-10's 30 runs did use OS randomness.
- **What this changes.** The rows at 0.95, 1.0 and 1.05 D are paired on seeds, not independent. Whether consecutive
  seeds give km-peeler effectively independent hash functions depends on its seeding, which I could not check because
  the km-peeler source is not in the repo.

The rest of C7 holds against the code. The source is unmodified, and the build adds only `-fcommon`
(`kp10_vs_iblt_job.sh`). r = 3 (`-r 3`, `:119-120`). The input was the canonical sets (`kmtools text` of
`sets/*.u64`). And "each party checks with its own data only" holds for kp10c (`kp10c.c:300-304`), though not for the
numpy runs (F3).

**Correction.** "…km-peeler with 30 fixed seeds (1001 to 1030) and the sketch with 30 seeds drawn from the operating
system."

## F7 -- "about 66 field multiplications per element" omits the decode, which costs about as much again on the GigaDB pair

**Claim:** C12 (NARROWED).

- **66 counts only the sketch.** m + c = 66 is the sketch cost: `psums()` does M = 66 `mulq` per element
  (`kp10c.c:200-221`).
- **The decode adds about 50 more on GigaDB.** It then runs Horner on P at every element of every bucket with a
  difference (`flag_roots`, `kp10c.c:228-237`), at about deg P + 1 multiplications per element. On GigaDB **every**
  bucket had a difference in round 1 (`run.log:28`, `buckets_with_difference_round1: 65536`). At 3,237,326 / 65,536 =
  49.4 differences per bucket on A's side, that adds about 50 multiplications per element of A, roughly 116 in all.
- **The measured CPU agrees.** Decode took 110.3 CPU-s against 97.7 for both sides' power sums (`timed.tsv:2`;
  `run.log:28`), and the results file puts the decode at about 30× the IBLT's diff and list (§5).

**Correction.** "The sketch needs m + c = 66 field multiplications per element to build, and decoding evaluates a
polynomial of degree up to m at each element of each bucket with a difference (about 50 more per element on the
GigaDB pair), against 3 cell updates."

## F8 -- the table's failed listings were large and correct, and the D/8 protocol threw them away

**Claims:** C2 and C13 (NARROWED).

C2 credits the gap at D/8 (20.8 vs 38.5) to the sketch resending only overfull buckets. The Introduction says "an
invertible Bloom lookup table that fails must be resent whole". The raw outputs show what a failed listing contained.
From `run.log` `iblt_trial` lines with `list_exit != 0`:

| pair, n | failed trials | mean correctly listed | spurious (max) |
|---|---|---|---|
| giga, 0.95 D | 30 | 1,353,288 (41.8% of D) | 0 |
| giga, 1.0 D | 12 | 2,021,206 (62.4% of D) | 0 |
| draft, 1.0 D | 22 | 4,701 (59% of D) | 0 |

A failed peel leaves a large, error-free partial listing. The table is linear too, so subtracting the listed items
and sending a smaller table for the remainder is the IBLT analogue of resending a bucket. The harness discards the
partial listing and resends at double size. The results file concedes "A smarter IBLT protocol … is not what
km-peeler implements" (§4). "Must" is therefore a property of the harness's protocol, not of the table. Even
without F1, the 38.5 is an upper bound on the table's cost, not its cost.

**Correction.** "…while the table, resent whole at double size after each failed listing (km-peeler's partial listings,
which were correct, were not reused), cost 38.5." In the Introduction, replace "must be resent whole" with "is, in
km-peeler, resent whole". An IBLT subtract-and-resend protocol is owed before any causal "because". It is a laila job,
and I did not submit it.

## F9 -- the 9.5% bound and "neither method returned a wrong difference" reproduce

**Claim:** C5 (SURVIVED).

- **The bound.** 1 − 0.05^(1/30) = 0.09503, which matches `failure_upper95` = 0.095 in both `kp10_seeds.tsv` and
  `iblt_trials.tsv`.
- **No wrong difference, KP-10.** `missing_total 0`, `spurious_total 0` over 30 + 30 runs.
- **No wrong difference, IBLT.** `undetected_failures 0` in every `iblt_trials` row, and every failed listing above
  had 0 spurious.
- **Independence.** The KP-10 seeds are independent OS draws. The IBLT caveat is in F6.

## F10 -- the "almost all from unplaced scaffolds" clause has no measurement and no ledger row

**Claim:** C10 (NARROWED).

- **The counts reproduce.** `truth.tsv:2-3` gives 3,238,101 and 7,931, with KMC 3.2.4 identical element for element
  (`run.log:19-20`, `same_amb`/`same_bma`). 3,238,101 / 442,151,894 = 0.732%.
- **The clause is unmeasured.** The only support is the prose of `kp10_vs_shibuya_2026-11-11.md:42` and `:193`, with
  no source line. No artifact assigns the 3,237,326 A∖B k-mers to Cannbio-2 records. Record and base counts
  (147 vs 10 records, 10,317,051 fewer bases; `sets.tsv`) make the clause plausible but do not measure it, and the
  775 B∖A k-mers show that the placed sequence differs too.
- **KMC is independent.** KMC is a separate counter run on the same FASTA, so its agreement is independent of kmtools.

**Owed (cheap on laila).** Rescan Cannbio-2 for the A∖B set and tally hits by record.

## F11 -- "2.4 times as fast as our hash" divides across two programs the paper says not to compare

**Claim:** C14 (NARROWED).

- **Where 2.4 comes from.** It is 24.67 / 10.24. The 24.67 is from `bench_e` (`paper_bench_2026-11-12.md:33`), and
  the 10.24 from `bench_kept` (`:51`). Ledger row 57 says "across the two programs".
- **The paper warns against exactly this.** Its Figure 3 caption says the panels "are not compared with each other,
  because the two programs give the same sampler code 28.7 and 22.5 cycles per kept k-mer", a 28% disagreement.
- **The same sampler differs between programs.** At 1/4 it takes 17.12 in `bench_e` and 15.43 in `bench_kept`
  (1.11×).
- **No within-program number exists.** `bench_kept` has no every-k-mer Eisenstein row, so the ratio cannot be formed
  inside one program.
- **The 1.12× is within one program and holds.** 11.43 / 10.24 = 1.116.

**Correction.** Either time the Eisenstein every-k-mer path in `bench_kept`, or write "about twice as fast (2.4 times
across our two benchmark programs, which differ by up to 28 percent on identical code)".

## F12 -- Figure 4B matches the tables, and plotting one sketch run against the table's mean does not mislead

**Claim:** C9 (SURVIVED).

- **The plotted values come straight from the tables.** `make_figures.py:249-262` reads `kp10_dest.tsv` and
  `iblt_protocol.tsv` directly. It asserts that every repetition ends exact with constant rounds, and both tables
  reproduce byte for byte from `run.log`.
- **The table's mean is every repetition's value.** Its protocol cost is deterministic given the rounds: the sizes
  are fixed by n, and the rounds are identical in all 10 repetitions (`mean_rounds = max_rounds`). The mean is
  therefore each repetition's value.
- **One sketch run is enough at the true D.** At the true D the sketch's 30-seed byte range is ±0.18%
  (36,415,840–36,547,696). The single run at D/8 has no spread measured, which the caption discloses.
- **The caption's protocol description matches the job.** It doubles n from 1.05 D_est after each non-zero
  `ibltseq list` exit. F1 and F8 bear on the figure's interpretation, not its values.

---

## Verdicts

| C1 | C2 | C3 | C4 | C5 | C6 | C7 | C8 | C9 | C10 | C11 | C12 | C13 | C14 | C15 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| REFUTED | NARROWED | NARROWED | NARROWED | SURVIVED | NARROWED | NARROWED | NARROWED | SURVIVED | NARROWED | REFUTED | NARROWED | NARROWED | NARROWED | NARROWED |

**Steering note (per the brief).** The claim list asks whether the byte accounting is "one side's message" but never
asks how many messages "recovered every difference" took. In the code it took two (F1), and that is the weakest
point of claims 4 to 6.
