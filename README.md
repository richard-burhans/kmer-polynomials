# Provably strand-symmetric polynomial hashing of k-mers

Richard Burhans (ORCID [0009-0006-5949-0958](https://orcid.org/0009-0006-5949-0958))

**Status:** draft, private. The abstract is provisional until the KP-10 vs Shibuya comparison is in.

## Abstract

<!-- The text between the markers is copied from manuscript/paper.md by `make -C manuscript readme`. Edit it there. -->
<!-- ABSTRACT:BEGIN -->
**Provisional abstract, pending the KP-10 sketch comparison.** Comparing k-mer content across genomes needs a hash that gives a k-mer and its reverse complement one value. The widely used rolling hashes ntHash and ntHash2 do this with fixed seeds, so their collision behaviour is measured, not proven. We read a k-mer as a polynomial over the prime field $\mathbb{F}_q$, $q = 2^{61} - 1$, with bases mapped to sixth roots of unity—$A = 1$, $C = \omega$, $G = -\omega$, $T = -1$, the $\theta = \pi/3$ member of the trigonometric encoding of Das, Nanda and Das. The cube root of unity $\omega$ lies in $\mathbb{F}_q$ because $q \equiv 1 \pmod 3$. Complement becomes negation and reversal a reflection of the polynomial, so the product $(F + \lambda)(R + \lambda)$ of the forward and reverse-complement evaluations is the same for both strands. Its unique factorization in $\lambda$ separates every pair of distinct double-stranded k-mers: two collide with probability at most about $2k/q$, $2.7 \times 10^{-17}$ at $k = 31$, for any input. The obvious alternatives fail, and on the 8,192 double-stranded 7-mers the sum $F + R$ takes 729 values and the product $F \cdot R$ takes 4,160. The same encoding makes a transition the cheapest substitution, at squared distance 1 against 3 for the A–C and G–T transversions. On a 2.0 GHz Xeon, a rolling implementation hashes every k-mer of the Cannabis sativa cs10 reference at 24.6 cycles per k-mer, 0.46 times the speed of a batched ntHash. A sampler that computes the provable hash only for the k-mers a packed-code filter keeps runs faster than batched ntHash at sampling rates of 1/19 and below, 1.21 times as fast at 1/1000.
<!-- ABSTRACT:END -->

## The encoding in one paragraph

A k-mer $s_0 \ldots s_{k-1}$ is the polynomial $p(x) = \sum_t s_t x^t$. Complement is negation, so the complement of a k-mer
is $-p(x)$, its reverse is $x^{k-1} p(1/x)$, and its reverse complement is $-x^{k-1} p(1/x)$. Two encodings are in the
repository:
- the fourth-root map $A = 1$, $C = i$, $G = -i$, $T = -1$, which needs $\mathbb{F}_{q^2}$
- the Eisenstein map $A = 1$, $C = \omega$, $G = -\omega$, $T = -1$, which lives in $\mathbb{F}_q$ itself and is the
  current choice.

## Read more

- **[manuscript/](manuscript/):** the paper ([paper.md](manuscript/paper.md)), its numeric ledger, the outline, and the
  build. `make -C manuscript pdf` builds paper.pdf with pandoc in a pinned docker image, `make -C manuscript lint` runs
  the writing-repo linter and reference validator, and `make -C manuscript readme` / `check` keep this abstract in sync
  with the paper.
- **[IDEAS.md](IDEAS.md):** ideas KP-1 to KP-14 with their status, the comparison with ntHash/ntHash2 (what is
  different and what is shared), and the open questions.
- **Key results** (all in [results/](results/)):
  - the strand-symmetric hash and why $F + R$ and $F \cdot R$ fail:
    [kp7_symmetric_check_2026-11-10.txt](results/kp7_symmetric_check_2026-11-10.txt),
    [notes/kp7_strand_symmetric_hash.md](notes/kp7_strand_symmetric_hash.md)
  - the Eisenstein hash: [kp7_eisenstein_2026-11-10.md](results/kp7_eisenstein_2026-11-10.md)
  - samplers against batched ntHash: [kp7_next_steps_2026-11-10.md](results/kp7_next_steps_2026-11-10.md),
    [kp7_kept_eval_2026-11-10.md](results/kp7_kept_eval_2026-11-10.md),
    [kp7_pack_hybrid_2026-11-10.md](results/kp7_pack_hybrid_2026-11-10.md)
  - $\mathbb{F}_{p^2}$ SIMD (AVX2, AVX-512): [kp7_m31_simd_2026-11-10.md](results/kp7_m31_simd_2026-11-10.md)
  - set reconciliation (KP-9, with KP-10 pending): [kp9_kompa_kompb_2026-11-10.md](results/kp9_kompa_kompb_2026-11-10.md)
  - negative and niche results: [ideas_spectral_2026-11-10.md](results/ideas_spectral_2026-11-10.md),
    [ideas_algebra_2026-11-10.md](results/ideas_algebra_2026-11-10.md),
    [ideas_eisenstein_2026-11-10.md](results/ideas_eisenstein_2026-11-10.md),
    [ideas_kp5_kp11_2026-11-11.md](results/ideas_kp5_kp11_2026-11-11.md)
- **[literature/](literature/):** the wantlists, including the Eisenstein novelty check
  ([eisenstein_wantlist.md](literature/eisenstein_wantlist.md)).
- **[hash/](hash/):** the C code. `make test` runs the test suites, `make mutation` checks that deliberately broken
  versions fail them, and `run_bench.sh` and `run_bench_pack.sh` run the benchmarks against ntHash
  (`build_nthash.sh` builds the pinned ntHash).
- **Hardware counters on laila** come from `scripts/laila_perf.py` in the separate cannabis-genome repository, which runs
  a counted job through a Galaxy user-defined tool.
