# Revisiting the ideas under the Eisenstein encoding (2026-11-10)

**Status: done. Encoding A = 1, C = ω, G = −ω, T = −1 (hash/kp7e.h). KP-8/9/10 are revisited separately in recon/.**

Code: `ideas/eisenstein/` (Makefile, dbg_e.c, spaced_e.{h,c}, test_spaced_e.c, bench_spaced_e.cpp, kp1_tstv.py,
kp13_compare.py, run_spectral.sh, run_laila.sh). KP-4 and KP-13 reuse `ideas/spectral/` with the new
`--encoding eisenstein` option. Raw output: `results/ideas_eisenstein/`; laila job af1eb52f65c6d471 (Xeon Gold 6330)
in `results/ideas_eisenstein/laila_20261110T234758Z/`.

In ℂ, Re(a·b̄) scores match +1, transition +½, A–C/G–T −½ and complement −1, and squared distance is 0 / 1 / 3 / 4.
Under the quarter map (A = 1, C = i, G = −i, T = −1) the scores are +1 / 0 / 0 / −1 and the distances 0 / 2 / 2 / 4.
The metric exists only in the complex embedding; the F_q hash keeps the algebra (complement = negation, ω²+ω+1 = 0).

## KP-1: the metric as a substitution score

**Verdict: transitions become the cheapest substitution, and they are the commonest real difference.**

`kp1_tstv.py` over the WF-I pairwise inputs for cs10 NC_044378.1 (`kp1_tstv_NC_044378.1.json`). KOMPa vs KOMPb is the
two haplotypes of one plant, compared only at cs10 positions both align.

| pair | columns | substitutions | Ts/Tv | transitions | A–C/G–T | complements | mean cost, quarter | mean cost, Eisenstein |
|---|---|---|---|---|---|---|---|---|
| cs10 vs KOMPa | 31.7 M | 933,307 | 1.585 | 61.3% | 20.2% | 18.5% | 2.37 | 1.96 |
| cs10 vs KOMPb | 32.4 M | 884,600 | 1.584 | 61.3% | 20.5% | 18.2% | 2.37 | 1.96 |
| KOMPa vs KOMPb | 29.1 M | 316,451 | 1.669 | 62.5% | 17.2% | 20.3% | 2.41 | 1.95 |
| same, uppercase only (outside soft-masked repeats) | 12.6 M | 144,065 | 1.847 | 64.9% | 17.3% | 17.8% | 2.36 | 1.88 |

Ts/Tv is 1.6–1.8, inside the 0.87–2 range in our notes; the 0.872 was raw HSPs, half noise. Under the quarter map a
transition costs the same as an A–C/G–T transversion (2). Under Eisenstein it costs half as much (1 against 3), so k-mer
distance ranks a neighbour that differs by the likelier substitution as closer. The mean cost of a real substitution
falls from about 2.4 to about 1.9 (complement = 4 on both scales). This holds only where distances are computed in ℂ
(KP-1, KP-4, KP-13), not for the F_q hashes.

## KP-4: truncated-DFT lower bound

**Verdict: unchanged. Still white noise; still not useful at k-mer scale.**

`kp4_dft_bound.py --encoding eisenstein` against `--encoding quarter`, same seed, cs10 NC_044378.1, k = 31.

- **Spectral energy:** for seeded neighbours (distinct k-mers sharing a 12-base prefix), the first 8 low-frequency
  coefficients carry 0.257 of the distance energy, against 0.258 for white noise (quarter 0.258). Random pairs: 0.296
  for both.
- **Hamming-safe pruning** (bound > 4d) on seeded neighbours, m = 8: d = 3 prunes 20.6% of 76.6% prunable (quarter
  21.0%); d = 5 prunes 3.9% of 70.7% (quarter 3.2%).
- **Weighted thresholds** (neighbour means squared distance ≤ W; new in this run, both encodings): the bound is tighter
  against the metric itself. W = 6 prunes 48.3% of 79.4% (quarter 49.8% of 80.1%); W = 10 prunes 28.9% of 73.6%.
  That still misses most prunable pairs, and exact distance on 2-bit codes (XOR + popcount, with a 16-entry class
  table for the weights) is cheaper than computing the bound.

## KP-13: FFT correlation as a diagonal finder

**Verdict: recall does not improve; it drops slightly. Not a lastz/yama pre-filter.**

`kp13_fft_diag.py --encoding eisenstein`, the same 40 windows as the quarter run (same `--rng`), pooled by
`kp13_compare.py`. Both encodings are pooled the same way here: hits/peaks, found/heavy clusters, found weight/heavy
weight. These pooled values differ from the per-window figures in results/ideas_spectral_2026-11-10.md, but each
row below compares like with like (`kp13_compare.txt`).

| variant | z | quarter precision / weighted recall | Eisenstein precision / weighted recall |
|---|---|---|---|
| complex | 6 | 0.685 / 0.703 | 0.703 / 0.673 |
| complex | 10 | 0.879 / 0.448 | 0.895 / 0.426 |
| complex, soft-mask zeroed | 6 | 0.577 / 0.783 | 0.581 / 0.762 |
| **complex, soft-mask zeroed** | **10** | **0.971 / 0.635** | **0.991 / 0.576** |
| 4 indicators (encoding-free; a sanity check) | 6 | 0.831 / 0.458 | 0.831 / 0.458 |

**Why:** with ~60% transitions, a diverged aligned column scores about the same mean under both maps (+½ transitions,
−½ and −1 transversions roughly cancel). But the per-base background variance rises from 0.5 to 0.625 (four score
levels, ±1 and ±½, instead of three), so z-scores fall a little. Precision at z ≥ 10 rises slightly (0.971 → 0.991) and
recall falls (0.635 → 0.576). The structural problem, short diagonals broken by indels, is untouched.

## KP-3: de Bruijn adjacency from hashes

**Verdict: still correct and exact; storage is now even with packed codes at k ≤ 32 and half at k > 32, but navigation
stays 4–6× slower.**

`dbg_e{31,47,63}` on the first 20 Mb of cs10, the same input as ideas/algebra/dbg.c. Laila timings:

| k | distinct k-mers | edges | false / missed | packed codes | hash only (Eisenstein F) | bytes/node packed vs hash |
|---|---|---|---|---|---|---|
| 31 | 13,530,525 | 13,608,457 | 0 / 0 | 204.9 ns/node | 1,210.7 ns/node (5.9×) | 8 vs 8 |
| 47 | 14,292,098 | 14,329,963 | 0 / 0 | 309.3 ns/node (128-bit) | 1,316.7 ns/node (4.3×) | 16 vs 8 |
| 63 | 14,712,782 | 14,734,711 | 0 / 0 | 312.2 ns/node (128-bit) | 1,325.5 ns/node (4.2×) | 16 vs 8 |

The non-branching node counts agree exactly at every k. The Gaussian version needed 16 bytes per node and broke even
on storage only for k > 64. The Eisenstein hash is 8 bytes, so hash-only graphs **halve** storage for 32 < k ≤ 64, as
predicted. The time cost is structural: 32 candidate lookups per node against 8 for packed codes. The provable
false-positive bound is ≤ 0.03–0.07 expected false edges over the whole set; zero were observed.

## KP-14: spaced seeds

**Verdict: twice as fast as the Gaussian version, 0.78× ntHash2's SeedNtHash, and still exactly strand-consistent
for any mask.**

`spaced_e.{h,c}`: one 61-bit multiply per strand per mask plus two table additions per block. Tests (`make test`):
- rolling equals direct evaluation for 6 masks with N, and strand symmetry holds;
- the all-ones mask is bit-identical to kp7e_sym over the kp7e reference roller;
- exhaustive k = 7 hashes equal the exact masked double-stranded classes (512, 8,192, 512, 2,080);
- `make mutation` kills the wrong-power mutant.

Laila, 58.8 Mb of cs10, 3 interleaved rounds, one hash per seed:

| masks | Gaussian | Eisenstein | SeedNtHash | Eisenstein / SeedNtHash |
|---|---|---|---|---|
| 1 | 11.0 M/s | 21.4 M/s | 27.8 M/s | **0.78×** (Gaussian 0.40×) |
| 2 | 10.1 M/s | 19.9 M/s | 28.2 M/s | **0.71×** (Gaussian 0.36×) |

Sensitivity (fraction of windows unchanged under substitutions at 0.001 / 0.01 / 0.05): with the non-palindromic mask,
Eisenstein gives 0.97150 / 0.74724 / 0.22572, exactly the exact double-stranded fraction, as Gaussian did. SeedNtHash
gives 0.98025 / 0.81795 / 0.35843 and warns that its reverse-complement hashing is inconsistent. With the palindromic
mask all are identical.

## KP-5 and KP-11

**KP-5 (RC palindromes): unchanged.** The test p(x) = −x^{k−1} p(1/x) relies only on complement = negation and the
reversal identity, and both hold under Eisenstein. The F_q version can now test it on a single field element.

**KP-11 (generating functions by k-mer code): unchanged.** It indexes k-mers by their base-4 code, not by the base
values, so the encoding never enters.

## KP-7 against ntHash: speed lines updated

The IDEAS.md comparison block now gives the Eisenstein figures against batched ntHash
(results/kp7_eisenstein_2026-11-10.md): provable every-k-mer hashing 0.46× (24.6 cycles/k-mer, scalar, 2 lanes); the
v2 sampler 1.21× at 1/1000 and 1.12× at 1/100; break-even near 1/19.
