# Provably strand-symmetric polynomial hashing of k-mers

Richard Burhans (ORCID [0009-0006-5949-0958](https://orcid.org/0009-0006-5949-0958))

**Status:** draft, private. The abstract is provisional until the KP-10 vs Shibuya comparison is in.

## Abstract

<!-- The text between the markers is copied from manuscript/paper.md by `make -C manuscript readme`. Edit it there. -->
<!-- ABSTRACT:BEGIN -->
Comparing k-mer content across genomes needs a hash that gives a k-mer and its reverse complement one value. A k-mer of up to 32 bases can be canonicalized as a machine word and then hashed, but longer k-mers, spaced seeds and sketches that add hash values need the symmetry inside the hash. The rolling hashes ntHash and ntHash2 provide it with fixed seeds and no proven collision bound. We read a k-mer as a polynomial over $`\mathbb{F}_q`$, $`q = 2^{61} - 1`$, with bases at sixth roots of unity, $`A = 1`$, $`C = \omega`$, $`G = -\omega`$, $`T = -1`$, a published encoding that fits $`\mathbb{F}_q`$ because $`q \equiv 1 \pmod 3`$. Complement becomes negation, so the product $`(F + \lambda)(R + \lambda)`$ of the two strand evaluations is the same for both strands. Distinct double-stranded k-mers collide with probability at most $`2(k-1)/q`$, $`2.6 \times 10^{-17}`$ at $`k = 31`$, for every input, while the sum $`F + R`$ takes only 729 values on 8,192 classes of 7-mers. Power sums of the hashes form a linear set sketch. Between two cannabis assemblies it needed 11.3 bytes per difference, against 20.5 for an invertible Bloom lookup table. On k-mers as written, it recovered every difference between an assembly and its differently oriented draft, at $`k = 31`$ and $`k = 63`$. The hash costs 24.7 cycles per k-mer, 0.46 times the speed of a batched ntHash.
<!-- ABSTRACT:END -->

## The encoding in one paragraph

A k-mer $`s_0 \ldots s_{k-1}`$ is the polynomial $`p(x) = \sum_t s_t x^t`$. Complement is negation, so the complement of a k-mer
is $`-p(x)`$, its reverse is $`x^{k-1} p(1/x)`$, and its reverse complement is $`-x^{k-1} p(1/x)`$. Two encodings are in the
repository:
- the fourth-root map $`A = 1`$, $`C = i`$, $`G = -i`$, $`T = -1`$, which needs $`\mathbb{F}_{q^2}`$
- the Eisenstein map $`A = 1`$, $`C = \omega`$, $`G = -\omega`$, $`T = -1`$, which lives in $`\mathbb{F}_q`$ itself and is the
  current choice.

## Read more

- **[manuscript/](manuscript/):** the paper ([paper.md](manuscript/paper.md), and the built [paper.pdf](manuscript/paper.pdf)),
  its numeric ledger, the outline, and the build. `make -C manuscript pdf` rebuilds paper.pdf with pandoc in a pinned
  docker image, `make -C manuscript lint` runs the writing-repo linter and reference validator, and
  `make -C manuscript readme` / `check` keep this abstract in sync with the paper and fail if the committed PDF is stale.
- **[IDEAS.md](IDEAS.md):** ideas KP-1 to KP-14 with their status, the comparison with ntHash/ntHash2 (what is
  different and what is shared), and the open questions.
- **Key results** (all in [results/](results/)):
  - the strand-symmetric hash and why $`F + R`$ and $`F \cdot R`$ fail:
    [kp7_symmetric_check_2026-11-10.txt](results/kp7_symmetric_check_2026-11-10.txt),
    [notes/kp7_strand_symmetric_hash.md](notes/kp7_strand_symmetric_hash.md)
  - the Eisenstein hash: [kp7_eisenstein_2026-11-10.md](results/kp7_eisenstein_2026-11-10.md)
  - samplers against batched ntHash: [kp7_next_steps_2026-11-10.md](results/kp7_next_steps_2026-11-10.md),
    [kp7_kept_eval_2026-11-10.md](results/kp7_kept_eval_2026-11-10.md),
    [kp7_pack_hybrid_2026-11-10.md](results/kp7_pack_hybrid_2026-11-10.md)
  - $`\mathbb{F}_{p^2}`$ SIMD (AVX2, AVX-512): [kp7_m31_simd_2026-11-10.md](results/kp7_m31_simd_2026-11-10.md)
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
