# KP-4 and KP-13: the spectral ideas, tested (2026-11-10)

**Verdicts:**
- **KP-4** (truncated-DFT lower bound for mismatch search): **not useful at k-mer scale.**
- **KP-13** (FFT correlation as a diagonal finder): **not a lastz/yama pre-filter.** At best it is a cheap
  high-precision anchor finder for already-paired regions.

Code: `ideas/spectral/` (common.py, kp4_dft_bound.py, kp13_fft_diag.py). Raw output: `results/ideas_spectral/`.

## KP-4: Parseval lower bound vs exact Hamming

Run: `kp4_dft_bound.py cs10.NC_044378.1.fa` (k = 31, 60.6 M valid k-mers; 200,000 pairs per set; raw output
`kp4_cs10_NC_044378.1.json`). LB_m = (1/k)·Σ over m coefficients of |P_j − Q_j|², which is ≤ ‖p − q‖². A pair within
Hamming d has ‖p − q‖² ≤ 4d, so LB_m > 4d prunes it safely (asserted: no true neighbour was ever pruned).

**The Voss note holds.** The share of ‖p − q‖² carried by the m lowest-|frequency| coefficients is about white noise:

| m | random pairs | seed-sharing pairs | white noise m/k |
|---|---|---|---|
| 1 | 0.051 | 0.036 | 0.032 |
| 4 | 0.158 | 0.132 | 0.129 |
| 8 | 0.296 | 0.258 | 0.258 |

**Pruning (m = 8, low-frequency order); each cell is pruned / could be pruned:**

| candidates | d = 1 | d = 2 | d = 3 | d = 5 |
|---|---|---|---|---|
| random pairs (Hamming mean 22.4) | 0.999 / 1.000 | 0.973 / 1.000 | 0.843 / 1.000 | 0.362 / 1.000 |
| distinct k-mers sharing a 12-base prefix (mean 9.4) | 0.590 / 0.888 | 0.398 / 0.814 | 0.210 / 0.766 | 0.032 / 0.707 |

**Why not useful:**
- On the candidates that matter (seed-sharing near neighbours), the bound misses most of the prunable pairs, and almost
  all of them at d = 5.
- For k ≤ 32, exact Hamming on 2-bit codes is one XOR, one OR/AND and one popcount on a single 64-bit word: a few
  cycles and 8 bytes per k-mer. The m = 8 bound needs 64 bytes per k-mer and about 32 floating-point operations per
  pair, and it is only a bound. (numpy timing, 25 vs 30 ns per pair, reflects numpy overhead, not C costs.)
- For indexing (avoiding all-pairs comparison), pigeonhole seeds (split into d + 1 parts, exact-match lookup) are
  exact and standard; low-dimensional DFT features have nothing to add.

## KP-13: FFT diagonals vs the KOMPa chr7 chains

Run: `kp13_fft_diag.py` with target cs10 NC_044378.1 and query KOMPa.chr7, chains from the 12-of-12 rerun (1,136
chains; the main one is on '−' and spans nearly the whole chromosome). 40 windows of 100 kb (20 '+', 20 '−'), each
correlated against ±500 kb of the query around the chain's median diagonal, on the chain's strand. Truth: chain gapless
blocks in the window, clustered by diagonal (±50 bp). "Heavy" means ≥ 1,000 aligned bases. Precision is the share of
peaks on a chain diagonal (±50 bp). The negative control is the same window against a same-sized query segment with
none of the window's chain diagonals. Raw output: `kp13_cs10_vs_KOMPa_chr7.json` and `kp13_band_cs10_vs_KOMPa_chr7.json`.

| variant | z | peaks (40 windows) | precision | heavy-cluster recall | aligned-base-weighted recall | negative peaks per window |
|---|---|---|---|---|---|---|
| complex (match − complement) | 6 | 219 | 0.685 | 0.442 | 0.608 | 0.90 |
| complex | 10 | 66 | 0.879 | 0.197 | 0.373 | 0.07 |
| complex, soft-mask zeroed | 6 | 350 | 0.577 | 0.571 | 0.700 | 6.33 (max 121) |
| **complex, soft-mask zeroed** | **10** | 136 | **0.971** | 0.410 | **0.564** | 0.05 |
| 4 indicators (exact matches) | 6 | 77 | 0.831 | 0.219 | 0.382 | 0.10 |
| complex, summed over ±25 diagonals | 6 | 40 | 0.950 | 0.131 | 0.296 | 0.00 |
| masked, summed over ±25 diagonals | 6 | 134 | 0.590 | 0.251 | 0.401 | 0.82 |

**Readings:**
- **Peaks are weak.** The median maximum z per window is 14 unmasked and 22 masked, against about 5.8 in the negative
  control. Cannabis haplotypes differ by frequent small indels, so an alignment is a run of short gapless blocks on
  shifting diagonals, while a correlation peak collects only one exact diagonal.
- **Summing over adjacent diagonals does not rescue it.** It adds the noise of all 51 offsets while most signal sits on a
  few, so recall falls.
- **The 4-indicator exact-match count is worse than the complex score** at the same z.
- **The best setting** (masked, z ≥ 10) gives 97% precision but recovers only about 56% of aligned bases per window,
  and misses most short or fragmented diagonals (heavy-cluster recall 41%).
- **Cost:** 0.165 s per correlation (100 kb against 1.1 Mb, FFT of length 2²¹, numpy complex64, single thread), against
  2.8–8.5 s for lastz on the same five window/segment pairs. But lastz returns gapped alignments (352–1,573 per pair),
  so the outputs are not comparable.
- **Whole-chromosome search** (a 100 kb window against all 65 Mb of KOMPa.chr7, both strands, overlap-save in 4 Mb
  chunks) costs 31–34 s per window, so about 6.5 CPU-hours for chr7 against NC_044378.1. One of three windows (ws
  34.67 Mb) was flooded by repeat peaks: 50 at z ≥ 10 per strand, 2 on chain diagonals. The other two put every
  z ≥ 10 peak on a chain diagonal.

**Verdict:** as a pre-filter that must keep recall (lastz seeding, or choosing yama's band), FFT correlation loses
30–60% of the aligned diagonals at any useful precision, so lastz would still be needed. A plausible niche: with
soft-masking, z ≥ 10 peaks are reliable anchors (97%), so the method could cheaply confirm that two candidate regions
are homologous, or seed band placement where a long collinear diagonal exists.
