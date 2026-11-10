# KP-7 on laila: hardware counters (2026-11-10)

Laila: Intel Xeon Gold 6330 at 2.0 GHz, kernel 6.1.0-52, perf_event_paranoid 2. Measured with
`cannabis-genome/scripts/laila_perf.py run` (UDT laila-perf, its own perf_event_open reader, user space only). Each
method ran in its own counted process: `bench cs10.fa 31 1 throughput only=<method>`, kp7 at 4fb12df, ntHash at c26bd45,
both built `-O3 -march=native` in the job. Raw logs: `laila_counters_2026-11-10/`.

Each count includes loading the 893 MB FASTA. A load-only run (`only=__load_only__`, three times) gave
6.94e9 cycles and 20.73e9 instructions (the cycle counts varied by under 2%). That was subtracted, then divided by
736,547,799 k-mers.

| method | cycles/k-mer | instructions/k-mer | IPC | branch misses |
|---|---|---|---|---|
| ntHash canonical | **23.3** | **110.1** | **4.73** | 0.22% |
| original `kp7_roll64` symmetric, x real | 136.3 | 271.1 | 1.99 | 6.99% |
| `kp7_roll_sym` all three options | 56.0 | 224.3 | 4.00 | 0.19% |
| `kp7_sym_batch`, x real | 51.8 | 158.5 | 3.06 | 0.34% |
| `kp7_sym_seq` x real, 1 lane | 50.1 | 160.2 | 3.20 | 0.31% |
| **`kp7_sym_seq` x real, 2 lanes** | **47.4** | 179.3 | 3.78 | 0.32% |
| `kp7_sym_seq` x real, 3 lanes | 52.6 | 183.3 | 3.48 | 0.33% |
| `kp7_sym_seq` x real, 4 lanes | 51.9 | 186.1 | 3.59 | 0.33% |
| `kp7_sym_batch`, x complex | 63.0 | 211.5 | 3.36 | 0.23% |
| `kp7_sym_seq` x complex, 2 lanes | 59.1 | 223.2 | 3.78 | 0.22% |

(IPC and branch-miss rates are whole-process figures, load included.)

Reading:
- The best version, 2 lanes, is 0.49× ntHash on this CPU. 3 and 4 lanes are slower, probably from register pressure.
- The gap is about 1.6× more instructions and 1.25× lower IPC. With IPC at 3.8 the loop is no longer strongly
  latency-bound, so further gains must cut instructions.
- The original per-call code mispredicted 7% of branches on this CPU, from the per-k-mer mode switch and restart
  checks. The batch rewrite removed that.
- ntHash's roll, from its source: one split-rotation plus two table XORs per strand, and canonical = fwd + rev, with no
  multiplies. Ours has about 6 widening multiplies per k-mer. The symmetric product (F + λ)(R + λ) alone is a full
  F_{q²} multiply.
- Fairness: ntHash ran through its per-k-mer `roll()` API, and ours through a batch API.

Next options (ranked): a cheap symmetric combiner mix(F) + mix(R) (heuristic, like ntHash); a hybrid that gates with
the cheap combiner and computes the provable (F + λ)(R + λ) only for kept k-mers; Mersenne-31 AVX2 or AVX-512 IFMA
SIMD; x = 2^r rotation (multiply-free rolling, but it loses the collision bound).
