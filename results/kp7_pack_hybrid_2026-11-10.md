# KP-7 hybrid v2: packed-code gate, provable hash for kept k-mers (2026-11-10)

**Verdict: at sketch rates the hybrid beats ntHash. At 1/1000 sampling it costs 11.6 cycles per input k-mer against
ntHash's 23.4 (2.0× faster), and every kept k-mer carries the provable KP-7 symmetric hash (F + λ)(R + λ),
bit-identical to kp7_roll64(SYM). The gate alone, a fast heuristic symmetric hash for every k-mer, is 1.26× faster
than ntHash. At high sampling rates (1/16) direct evaluation dominates, and hybrid v1's rolling is better.**

Idea (from KMC, refresh-bio/KMC; ideas only, KMC is GPL-3 and no code was copied): for k ≤ 32 a k-mer packs exactly
into 64 bits, 2 bits per base, and both strands roll with shifts and masks:
- fwd' = ((fwd << 2) | c) & mask
- rc' = (rc >> 2) | ((3 − c) << 2(k−1))

with codes A0 C1 G2 T3, so complement is 3 − c. The packed codes are exact. The gate is mix(fwd) + mix(rc) mod 2^64,
strand-symmetric because reverse complementing swaps fwd and rc, with heuristic collisions as for ntHash's fwd + rev.
A k-mer with gate < threshold gets the provable hash by DIRECT evaluation (kp7_direct + kp7_sym), so no field
multiply runs on the rolling path at all. Code: hash/kp7pack.{h,c} (75f76d5), bench: hash/bench_pack.cpp,
hash/run_bench_pack.sh (9ceb922).

Tests (`make test_pack`, `make mutation_pack`):
- the gate equals the reference from directly built codes at every window
- gate(s) == gate(rc s) for k = 1..32
- exactly one gate value per double-stranded class at k = 7 (8,192) and k = 8 (32,896)
- the kept set == {gate < thr}, with provable hashes and positions equal to kp7_roll64(SYM), for x real and x complex,
  Karatsuba on and off, every N layout and block size of test_kp7
- the mutants (an asymmetric gate, a wrong rc shift) are killed

## Laila hardware counters

Xeon Gold 6330, via cannabis-genome/scripts/laila_perf.py (job 061106f07b427759). One process per method:
`bench_pack cs10.fa 31 1 throughput only=<method>`, built `-O3 -march=native` in the job, ntHash pinned at c26bd45.
The same job's `__load_only__` run (7.74e9 cycles, 20.73e9 instructions) was subtracted, then divided by 736,547,799
input k-mers. Raw logs: `laila_counters_pack_2026-11-10/`.

| method | cycles/k-mer | instructions/k-mer | IPC | vs ntHash |
|---|---|---|---|---|
| ntHash canonical | 23.4 | 110.1 | 4.71 | 1.00× |
| **pack gate**, every k-mer (heuristic symmetric) | **18.6** | 58.9 | 3.17 | **1.26×** |
| pack sample 1/16 (provable for kept) | 116.3 | 320.1 | 2.75 | 0.20× |
| **pack sample 1/1000** (provable for kept) | **11.6** | 46.6 | 4.03 | **2.02×** |
| kp7_sym_seq, 2 lanes (provable, every k-mer) | 47.6 | 179.3 | 3.77 | 0.49× |

Earlier, for comparison: hybrid v1 (rolling field state, gate in the field), sample 1/1000 = 33.4 cycles/k-mer and
sample 1/16 = 37.6 (kp7_laila_counters, the hybrid run).

Reading:
- **Direct evaluation costs about 1,700 cycles per kept k-mer:** (116.3 − 11.6) / (1/16 − 1/1000). That is 2k general
  F_{q²} multiplies through kp7_direct. Against hybrid v1's flat ~34–38 cycles, v2 wins below a sampling rate of
  about 1/65 and v1 wins above it. Typical sketch rates (FracMinHash 1/100 to 1/2000, minimizer-based sketches) are in
  v2's range.
- **The gate's density and bit uniformity match ntHash.** A local smoke run on 3 Mb of cs10 gave minimizer density
  ratios to 2/(w+1) of 1.0009, 1.0034 and 1.0040 at w = 10, 20 and 50 (ntHash: 1.0019, 1.0016 and 1.0037), and a
  maximum per-bit deviation of 0.00131 (ntHash 0.00135). Full-cs10 stats: `run_bench_pack.sh 31 3 "" <label>`.
- **The load baseline drifts between jobs** (6.9e9 to 7.7e9 cycles), so only same-job differences are reliable.
  ntHash's per-k-mer figure has been stable across jobs: 23.3, 23.6 and 23.4.

Next:
- Make direct evaluation cheaper for x real. kp7_direct uses general F_{q²} multiplies, while a real x needs half,
  and Horner over the packed code could use a precomputed 4^m-entry table per m-base chunk (another KMC-style idea).
  That moves the v1/v2 break-even upward.
- For k > 32, two-word packed codes.
