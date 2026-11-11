# Leading k-mer tools: where KP-7/KP-10 could help, and where they cannot (2026-11-11)

**Status:** survey only. Source was read and nothing was benchmarked. There are three small correctness tests,
programs in `/home/agent/.claude/jobs/0d0eac3e/tmp/survey/`: KmerStream, Metagraph and Bifrost's hash collisions, and
an exhaustive check of GGCAT's ntHash. The pangenome rows come from a delegated source read. Scope: the 15 requested
tools plus LexicMap, sourmash and 10 pangenome tools. Clones are in `bench/data/kmer_tool_survey/`; each citation gives the commit as `repo@hash` and
then `file:line`.

**Our properties, as used below.** The hash is S = (F+λ)(R+λ) over F_q, q = 2⁶¹−1, with the Eisenstein encoding.
- (B) **Bound:** the collision bound is per pair, about 2(k−1)/q over random x and λ.
- (Sy) **Symmetry:** strand symmetry holds by factorization, for any k.
- (L) **Linearity:** field values give linear, subtractable power-sum sketches with exact difference recovery: 11.3
  B/diff, orientation-invariant, and deterministic overflow detection by check sums. Per the round, that detection is
  probabilistic in the hash, so call it "detected with probability ≥ 1 − O(m/q)", not "deterministic".
- (O1) **Updates:** O(1) roll, substitution and concatenation (ntHash shares these).
- (SS) **Spaced seeds:** strand-consistent for any mask, at a sensitivity cost: the mirror is unioned in, 0.99^29
  against 0.99^20 for a 20-of-29 mask.

Cost: 24.7 cycles per k-mer on every k-mer, against 11.37 for batched ntHash. The kept-only sampler beats batched ntHash
below about 1/20.

**Honesty constraint.** S is a degree-2(k−1) polynomial in (x, λ). It is NOT k-wise (or even pairwise-uniform)
independent in the sense that HyperLogLog, MinHash, SetSketch or count-sketch analyses assume. That gap closes for none
of the tools here, but none of them meets the assumption either. Every tool surveyed uses a fixed or seeded mixer
(Wang, Murmur3, ntHash, a cyclic polynomial), and none of these is k-wise independent; tabulation-style hashing is the
standard way to get provable sketch guarantees (Pătraşcu–Thorup). So "S then a bijective mixer" keeps (B), because a
bijection preserves distinctness, and is exactly as heuristic as everyone else for sketch uniformity. Never claim more.

## Per-tool table

"Hash / canonical" is from source. "Bottleneck?" means a published profile or measurement of hashing's share. UNMEASURED
means none was found, and we did not guess.

| tool | hash, and RC canonicalization (source) | hashing used for | hashing a measured bottleneck? | fit of our properties |
|---|---|---|---|---|
| **Jellyfish** (gmarcais/Jellyfish@43b1ab2) | Random invertible GF(2) matrix applied to the 2-bit code: `include/jellyfish/large_hash_array.hpp:292,317`; built at `:959` via `randomize_pseudo_inverse`, product at `rectangular_binary_matrix.hpp:156`. Canonical = lexicographic min of code and RC: `mer_dna.hpp:416-431`, iterator `mer_iterator.hpp:51`. | Slot in a lock-free open-addressing table with quadratic reprobes. Because the matrix is invertible, only the key bits not implied by the slot are stored (quotienting). | UNMEASURED. The paper's case is lock-free CAS and memory, not hash cost. | **No fit.** The hash is a bijection on keys, so there are no collisions to bound. Quotienting needs invertibility, which S lacks. |
| **KMC 3** (refresh-bio/KMC@751ef36) | Counting uses **no hash**. The k-mer canonical is a lexicographic min: `kmc_core/kb_collector.h:74-95`. Canonical signatures (minimizers) are a lexicographic min over "allowed" m-mers: `kmc_api/mmer.h:40-60,81-89`, binned by a frequency-balanced map at `kmc_core/s_mapper.h:263`. Counting is a radix sort (raduls). One hash: ntHash v1, **min(f, r)**, at `kmc_core/libs/ntHash/ntHashWrapper.h:292-296`. | The ntHash feeds an ntCard-style histogram estimate (`ntHashWrapper.h:306-316`, `kmc.h:1119-1121,1162`), used to pick parameters and for estimate-only runs. | UNMEASURED for the hash. KMC's costs are I/O and sorting. | **No fit** for counting, which has no hash. The estimator is ntCard's case (below). |
| **meryl** (marbl/meryl@3b3866c + meryl-utility@d2110f4) | **No hash.** Canonical = min of codes in ACTG order: `src/meryl/merylOp-count.C:542`, `merylOp-countThreads.C:246`, order at `kmers-v1/kmers-tiny.H:132,151-167`. | Prefix bits of the canonical code pick the bucket/file; suffixes are sorted (`merylOp-countThreads.C:249-256`, `merylCountArray.H:31-37`). | UNMEASURED (no hash). | **No fit.** Exact and cheaper on packed codes. |
| **DSK** (GATB/dsk@3737094, engine gatb-core@28b240a, HEAD rather than DSK's pinned db4e7da) | Canonical = min(code, RC code): `gatb-core/src/gatb/kmer/impl/Model.hpp:294`. Partitions by canonical minimizer (`Model.hpp:340+`). Per-partition `Hash16` table using `oahash64` (`tools/math/NativeInt64.hpp:191-202`); `hash64` (`:175-188`) for Bloom filters. DSK 2013 split iterations and partitions by hash(k-mer) mod n. | Partitioning, in-partition counting, Bloom filters. | UNMEASURED. | **No fit.** Keys are stored and compared exactly, so a hash collision only costs a probe. |
| **ntCard** (bcgsc/ntCard@f6ceae1) | ntHash v1, canonical = **min(f, r)**: `vendor/ntHash/nthash.hpp:274-279`; multi-hash `NTMC64` at `:381-390` uses h[0] = canonical. | **Cardinality / histogram estimation.** The top s = 11 bits select 2 sampled sub-streams, and the next r = 27 bits index the counters (`ntcard.cpp:133-145`); the estimator is at `:258-270`. Every k-mer must be hashed. | UNMEASURED share. Hashing every k-mer is the algorithm's only per-k-mer work, so it is plausibly dominant, but no profile was found. | **Weak fit, speed only.** It keeps about 1/2¹⁰ of k-mers, so our kept-only sampler (1.21× batched ntHash at 1/1000) could feed it. But (B) is irrelevant to an estimator, which needs uniformity, and S is not k-wise independent. ntHash isn't either. |
| **KmerGenie** (official 1.7051 tarball, kmergenie.bx.psu.edu; no GitHub repo) | By default it runs ntCard per k (`kmergenie:188-207`). Legacy `specialk`: RanHash (Numerical Recipes) of the canonical min (`minia/Hashing.cpp:29-40`, `minia/Kmer.cpp:122`), sampled by `hash % reduction_factor` (`main.cpp:62-69`). | Histogram estimation across many k, by sampling. | UNMEASURED. | As ntCard. |
| **KmerStream** (pmelsted/KmerStream@7c33d96) | RepHash, a cyclic polynomial over 128-bit rotations, with random tables seeded from time() when seed = 0 (`RepHash.cpp:4-17`). Canonical = **forward XOR reverse-complement** hash: `RepHash.hpp:92-94`; roll at `:111-125`. | Streaming F0/f1 estimation. The trailing-zero level plus bits feed `StreamCounter.hpp:64-80`, called from `KmerStream.cpp:593`. | UNMEASURED. | **Concrete defect, our (Sy) fixes it.** The XOR canonical is the paper's F+R failure in GF(2) form. For odd k the middle term is rol^m(T[b] ⊕ T[comp b]), so a k-mer and the same k-mer with its middle base complemented collide **deterministically**, under every seed. **Verified 10,000/10,000 at k = 31** (`tmp/survey/rephash_mid.cpp`). No issue reports it (issues #1-#16 checked 2026-11-11). See idea 2. |
| **KmerLight** (nsivad/kmerlight@56d7125) | Canonical = lexicographic min: `src/KmerHash.h:123` (k ≤ 32), `:160` (k ≤ 64). Then MurmurHash3_x64_128 with a random seed (`CountSketchInstance.h:112`), level = trailing zeros (`:114`), then multiply-add-shift (`:120`). | Count-sketch-style per-level abundance estimation (`CountSketchInstance.h:124-136`). | UNMEASURED. | **No fit.** The analysis needs hash independence, which we can't supply, and its canonical step is exact already. |
| **Kmer-db** (refresh-bio/kmer-db@dc74e04). "kmer0db" is almost certainly this: 0 sits next to - on the keyboard. The other candidate is Ralston's `kmerdb` (PyPI), a k-mer count database. | Canonical = lexicographic min: `src/kmer_extract.h:84`. Hash = fixed-constant Murmur3-style mix plus fmix64 (`src/filter.h:96-114`). | **FracMinHash-style sampling:** keep a k-mer iff its hash falls in [start, start + fraction) (`filter.h:48-51`), for Jaccard/containment over many genomes. Table hashing in `hashmap_lp.h:36-44`. | UNMEASURED. The paper reports end-to-end time against Mash only (44× faster at 2% sampling). | **Partial.** Its sampling on packed codes (k ≤ 31) is already collision-free and cheap. Our sampler matches that and adds (Sy) by algebra, which it already has by lexicographic min. The real fit is (L): a subtractable per-genome sketch for near-duplicate genomes. See idea 3. |
| **Mash** (marbl/Mash@5df1f98) | Canonical = `memcmp` of the ASCII k-mer and its RC: `src/mash/Sketch.cpp:571`. MurmurHash3 (x64_128, or x86_32 when not use64) with seed 42, recomputed from scratch per k-mer: `src/mash/hash.cpp:10-23`, call at `Sketch.cpp:574`. | Bottom-s MinHash, giving the Jaccard → Mash distance. | UNMEASURED share. A non-rolling Murmur over k bytes per k-mer is plausibly the main per-k-mer cost. Dashing and Kmer-db report faster end-to-end sketching but not a profile. | **Weak.** Bottom-s needs the hash order to be min-wise independent-ish. S is not k-wise independent, so a mixer must follow, and then only (B) and (Sy) survive. Neither changes Mash's output, because memcmp canonicalization is already exact. |
| **Dashing** (dnbaker/dashing@0635bea, bonsai@534a5f8 HEAD rather than the pinned 5dfaf38, sketch@ba6e3f8) | Three encoders (`src/sketch_and_cmp.h:185-190`): (a) bonsai's 2-bit canonical, the lexicographic min at `bonsai/include/bonsai/kmerutil.h:156-159`, then WangHash in `addh` (`sketch/include/sketch/hash.h:54-65`); (b) ntHash v1 min(f, r) (`bonsai/include/bonsai/encoder.h:432-444`); (c) Rabin–Karp/cyclic rolling. | HyperLogLog (default) and other sketches for Jaccard and containment. | UNMEASURED share. | **No fit for HLL.** HLL's error analysis assumes ideal hashing. A bijective Wang on exact canonical codes is collision-free at k ≤ 32 and cheaper than S. |
| **Dashing 2** (dnbaker/dashing2@3906ebd) | k ≤ 32: bonsai 2-bit canonical; k > 32: `RollingHasher` with canonical = **min(cyclic f, cyclic rc)** (`bonsai/include/bonsai/encoder.h:904,948`); dispatch at `src/fastxsketch.cpp:399-421`. Keys get a seeded XOR mask plus WangHash (`src/enums.h:136-140`). | SetSketch, ProbMinHash, BagMinHash and one-permutation hashing, with multiplicities. | UNMEASURED. | **Narrow fit, k > 32 only.** There the canonical key is a 64-bit min of two cyclic hashes, with no bound, so distinct long k-mers can merge silently. S is an 8-byte key with (B) and (Sy) for any k, and a mixer can follow for the sketch. Whether that changes any Dashing 2 output is UNMEASURED. |
| **Kitsune** (natapol/kitsune@b4effe3) | Runs `jellyfish count`, with `-C` (canonical) **only if `--canonical` is given; default off** (`kitsune/modules/kitsunejf.py:191-195`). | Optimal-k selection (CRE, ACF, OFC) from Jellyfish counts. | UNMEASURED. | **No fit** (it delegates to Jellyfish). Note: by default it counts strand-specifically, so it is orientation-dependent on assemblies. |
| **minimap2** (lh3/minimap2@3c28777) | Canonical = lexicographic min of the 2-bit codes; k-mers equal to their own RC are skipped (`sketch.c:108-109`). `hash64` is Wang's invertible mix masked to 2k bits (`sketch.c:28-38`), and the key is `hash << 8 \| span` (`sketch.c:112`). | (w,k)-minimizer selection (random order) and index bucketing. Because the hash is a bijection, the index can recover the k-mer. | UNMEASURED for the hash. Mapping time is dominated by chaining and alignment (Li 2018). | **No fit.** It is collision-free by bijection, and density is equal for any good hash (our cs10 test). (SS) costs sensitivity, and minimap2 uses no spaced seeds. |
| **LexicMap** (shenwei356/LexicMap@0cc0839, lexichash@f58087b) | **No canonicalization.** For each of 20,000 random masks (`lexicmap/cmd/index.go:574`), it keeps the k-mer minimizing `kmer XOR mask` on each strand separately and records the strand in the location's low bit (`lexichash/lexichash.go:589-612` forward, `:614-640` reverse). Wang `Hash64` is used only to generate the masks (`lexichash.go:330`, `util.go:25-33`). k ≤ 32, default 31 (`index.go:120,178,571`). | **Probe seeding:** XOR-minimum means longest common prefix with the mask, so seeds are variable-length prefix matches. Suffix matching is a second index. | UNMEASURED for hashing. The cost is the mask-compare loop, a 3-base prefix index of masks (`m3`, `lexichash.go:591,620`). | **No fit for the core.** LexicHash's value IS lexicographic order: XOR with a mask preserves prefix length, and any field hash destroys it. Strand symmetry is already exact, by scanning both strands at 2× compare cost. (SS) would replace variable-length prefix seeds with fixed spaced seeds, a different seeding scheme with the 0.99^29-type cost. For k > 32 the limit is the uint64 key, and S can't extend a lexicographic prefix order. The one conceivable use, an 8-byte exact key for a k > 32 seed, is unneeded, since seeds are verified by alignment anyway. |

### sourmash (added on scope extension)

| tool | hash, and RC canonicalization (source) | hashing used for | hashing a measured bottleneck? | fit of our properties |
|---|---|---|---|---|
| **sourmash** (sourmash-bio/sourmash@43f64cb, Rust core; branchwater/manysketch call the same core) | Canonical = byte-wise `min(kmer, revcomp)` of the uppercase ASCII k-mer: `src/core/src/signature.rs:485-487`. Hash = first 64 bits of MurmurHash3_x64_128, default seed 42 (Mash-compatible): `src/core/src/lib.rs:69-70`. Non-rolling: k bytes are rehashed per k-mer. | **FracMinHash:** keep the hash iff it is ≤ max_hash = 2⁶⁴/scaled (`src/core/src/sketch/minhash.rs:24-30`, test at `:356`); default scaled = 1000. Then containment, ANI, `compare` and `gather` (min-set-cover) run on the kept hash sets. Hashes are stored in public databases (GTDB and GenBank sketches, and the SRA metagenome index behind branchwater). | **Partly.** Issue #860 (2020) posted a profile image (not read here); the fix it led to was DNA validation, not the hash. khmer PR #1792 reported a 40-50% speed-up from a rolling hash, quoted by ctb in #1659. No sourmash profile of the Murmur share was found, so UNMEASURED for sourmash itself. | **(a) Sampler:** a real speed fit, but **incompatible**. FracMinHash's kept set is defined by Murmur(k-mer) ≤ max_hash, and Murmur isn't monotone in any packed-code gate, so no pre-filter can reproduce it. Our sampler (about 9.4 cycles/k-mer at 1/1000) gives a *different* kept set, which would mean re-sketching every public database. That is only acceptable for private or new collections. **(b) Uniformity:** not provably met. At a fixed x, the map λ ↦ (F+λ)(F+λ+d) is 2-to-1, so S's marginal isn't uniform on F_q, and pairs aren't independent. The keep rate is about 1/s and **must be checked empirically**. Murmur is unproven too: sourmash's estimators assume an ideal hash. With a mixer after S it is equally heuristic, not better. **(c) Power sums:** a fit that **keeps compatibility**. KP-10 needs only field elements, so power sums of sourmash's *own* Murmur hashes (top 61 bits) give a linear, subtractable sketch of an existing FracMinHash signature. Our hash is then optional: orientation is already handled by the min() canonical. See idea 3. `gather` itself is no fit: it needs many-reference overlap queries, and power sums answer only "decode A △ B when it is small". |

### Pangenome tools (added on scope extension)

These rows come from a delegated source read: one agent, file:line checked by it, with test programs in
`tmp/survey/` (`cyc_test.cpp`, `ggcat_nt.c`, `mx3inv.py`). Submodules were followed at their pinned commits under
`bench/data/kmer_tool_survey/_deps/`. **No tool or paper profiles the hash itself, so column 4 is UNMEASURED for
every row.** No tool here uses a sum/XOR-of-strands canonical, so none has the KmerStream defect. Bifrost's embedded
KmerStream copy hashes the sorted pair {h, ht} instead, and avoids it: XOR of its own h and ht collides 10,000/10,000,
its `hash()` 0/10,000.

| tool | hash, and RC canonicalization (source) | hashing used for | fit of our properties |
|---|---|---|---|
| **Panagram** (panagram@58a3233, KMC fork skovaka/KMC@16c3988) | No hash of its own. Per-genome `kmc` gives a one-hot `set_counts` database (`workflow/Snakefile:101-108`), summed into ≤ 32-genome bitvector databases (`index.py:407-426`) and queried through py_kmc_api (`index.py:932-936`). Canonical = lexicographic min (`kmc_api/kmc_file.cpp:1228,1237`). Signatures as in KMC (`mmer.h:81-90`). | None in the main path: lexicographic and sort-based. | **Partial, idea 1.** The canonical step: no fit. A per-genome power-sum sketch beside the one-hot databases gives exact differences between near-identical assemblies or haplotypes without rescanning the bitmaps. |
| **PanKmer** (gitlab salk-tm/pankmer@b94dc99) | Canonical = max(fwd, rc) under A=3, C=2, G=1, T=0, which is the lexicographic min in ACGT order (`rust/src/decompose_kmers.rs:42-69`). mix64 (`helpers.rs:100-108`) on every canonical k-mer for downsampling (cutoff at `get_kmers.rs:255`). FxHashMap tables (`lib.rs:55`). k odd, 8 < k < 32 (`pankmer.py:699-704`). | Downsampling and hash tables. Memory blocks split by value range. | **Sampling: no fit**, since mix64 on a u64 it already holds is cheaper than S. **Index diffs: plausible**, idea 1. |
| **kmtricks** (kmtricks@9bccc77) | Canonical = lexicographic min (`include/kmtricks/gatb/sorting_count.hpp:426-431`). Bloom bit = XXH64(canonical) mod w + w·p (`:347-357`). GATB frequency-or-lex canonical minimizers (`task.hpp:278-283`, `gatb-core-stripped/.../Model.hpp:294,911,1051`), mapped to partitions by `XXH64(m) % nb_partitions` (`repartition.hpp:46-53`). | Partitioning and the k-mer → Bloom-bit map (HowDeSBT-style). | **Good architectural fit, idea 1.** Per-partition power sums add, so they could be emitted at merge time for exact sample-vs-sample or run-vs-run differences, and they survive hash mode, which loses k-mer identity. (Dead code only: `kmer_hash.hpp:50-66,239-241` has a word-count-as-bytes length and a swap-symmetric XOR fold, never built.) |
| **Bifrost** (bifrost@eaf7506) | k-mer hash: wyhash of packed `rep()` (`src/Kmer.hpp:120-122`), where `rep()` = lexicographic min (`Kmer.cpp:109-114`). Minimizer hash: RepHash with canonical `wyhash(sorted{h, ht})` (`RepHash.hpp:92-99`; 128-bit state for k > 64 at `:146-318`). | Minimizers (`CompactedDBG.tcc:744-752`), blocked Bloom filter (`:3918-3921`), k-mer table (`KmerHashTable.hpp:261`), color-set slots (`DataStorage.tcc:380`, `ColorSet.cpp:1028-1077`), embedded KmerStream (`KmerStream.hpp:126-165`). | **Canonical: no fit** (it is exact and strand-safe already). **Per-color sketch:** only across separate graphs or versions. Within one graph, walking the graph is exact. |
| **GGCAT** (ggcat@806405e) | Minimizers: an ntHash variant, per-base (2c+1)·0x397f178c6ae330f9, canonical = min(fw, rc) << 1 (`crates/hashes/src/nthash_base.rs:49-61`, `cn_nthash.rs:29-31,52-56,96-99`). k-mer key: packed code (k ≤ 64) or RabinKarp128 mod 2¹²⁸ (k > 64), both canonical by min (`api/src/utils.rs:17-65`, `cn_seqhash_base.rs:104-106`, `cn_rkhash_base.rs:148-149`). Color-set dedup: SipHash-1-3, 128-bit, with the hash as the key and **no verification** (`colors/src/colors_memmap_writer.rs:65-69`). | Minimizer bucketing, merge hash maps, color-set dedup. | **No fit.** The SIMD minimizers cost less than 24.7 cycles, and the key is the exact packed code. Its ntHash canonical had 0 collisions, exhaustively, for m = 10-13. Unverified color-set dedup by a 128-bit hash is sound in practice, and our bound wouldn't change it. |
| **Cuttlefish 2** (v2.2.0@30eef21; main is now Rust "Cuttlefish 3"@c6f27f2) | BBHash MPHF (`include/Kmer_Hash_Table.hpp:28`, `src/Kmer_Hash_Table.cpp:92`) over XXH3_64 of packed bytes (`Kmer_Hasher.hpp:19-21`, `Kmer.hpp:260-264`). Canonical = min(kmer, rc) (`Kmer.hpp:571-573`). k-mers come from KMC. Cuttlefish 3: minimizer canonical = min(wyhash(fwd), wyhash(rev)) (`partition.rs:1702-1704`). | The MPHF over the state table. | **No fit.** Keys all come from the KMC set, so an MPHF needs mixing, not a collision bound. |
| **Metagraph** (metagraph@353cc00; SSHash fork ratschlab@4b9fe58) | BOSS graph built by sorting, no hash. Optional Bloom: Lemire CyclicHash, 64-bit words, canonical = min(fwd, rev) (`kmer/kmer_bloom_filter.cpp:34-51`, `common/hashers/rolling_hasher.hpp:14-19`), with Kirsch-Mitzenmacher probes (`common/bloom_filter.cpp:27-48`). Hash graphs use std::hash of the bitset (`hash.hpp:11-17`). SSHash backend: mixer_64 minimizers plus CityHash128 PTHash. | Query Bloom prefilter (`cli/query.cpp:937`), hash graphs, row dedup, SSHash. | **Concrete defect, our hash fixes it.** CyclicHash rotates mod 64 (`cyclichash.h:33,56-58,100-105`), so for **k > 64**, swapping the bases at positions i and i+64 leaves both strand hashes unchanged. **Tested 1,000/1,000 at k = 65, 81 and 127; 0/1,000 at k = 63.** The impact is only extra graph lookups, since Bloom hits are verified. A field polynomial has no rotation period. See idea 4. |
| **Themisto** (themisto@af37939) | No hash of its own: `std::hash<Node>` = node id (`include/DBG.hh:299-303`). Search goes through SBWT. Both strands are indexed, with no canonicalization (`pseudoalign.hh:147`). Colors come from GGCAT (`build_index_main.cpp:452-466`) or an external sort (`Coloring_Builder.hh:481-510`). | Only inside GGCAT. | **No fit** (rank-based index, both strands explicit). |
| **Fulgor** (fulgor@ef65419; sshash@d85762b, sketch@f7c9d3b) | SSHash: mixer_64 minimizer order (`sshash/include/hash_util.hpp:88-114`); between strands it takes the min of the raw m-mers (`dictionary.cpp:38-45`); PTHash with murmurhash2_128; canonical k-mer = lexicographic min (`dictionary.cpp:65`). HLL (WangHash, p = 10) over unitig ids per reference, for meta and differential clustering (`build_util.hpp:16,31-32,106-111`). | Minimizers, two MPHFs, HLL clustering. | **No fit.** The MPHFs are verified, and HLL clustering needs similarity at large differences, where power sums don't decode. The differential color diff is already an exact small integer-set merge. |
| **SBWT** (SBWT@b6e6830) | No k-mer hashing in the index and no canonicalization: RCs are added as inputs (`src/CLI/sbwt_build.cpp:107-121`), and KMC runs non-canonical (`run_kmc.cpp:678`). | Only the small in-memory builder's unordered_set (`Kmer.hh:232-258`). | **No fit** (as Themisto). |

**Caveat for every sketch idea.** KP-10 decodes **hash values**, not k-mers. Turning them back into k-mers takes one
rolling re-scan of the side that holds each value, which is cheap but not free. A pure-sketch pipeline can report
exact difference counts and hash sets without the sequences, but not the k-mers themselves.

## Ranked combination ideas (at most 5)

Novelty searches were all run with WebSearch on **2026-11-11**; the queries are listed under each idea.

### 1. Orientation-invariant, subtractable per-genome sketches for pangenome k-mer matrices (Panagram/PanKmer/kmtricks)

- **Claim.** Store one KP-10 power-sum sketch per assembly beside the presence/absence matrix. Then:
  (i) the exact k-mer difference between any two genomes, or between one genome and the sum of a clade, decodes in
  bytes proportional to the difference (about 11.3 B/diff), without the matrix;
  (ii) a re-assembly or draft of the same genome gives a near-zero difference whatever its contig orientation;
  (iii) adding or removing a genome from an aggregate sketch is one addition or subtraction.
  PanKmer and Panagram rebuild or query the full |K|×|G| table for these. kmtricks could emit per-partition sums
  at merge time, because sums add across partitions.
  The output is the differing **hash values** and counts. Turning them into k-mers takes one rolling re-scan of the
  owning assembly.
- **Refutation experiment.** Take the cannabis panel assemblies already in KMC form. Over all pairs (or a 20-genome
  subset), compare decode size, time and exactness against the exact differences from KMC `kmc_tools` subtract. It is
  refuted if, for typical panel pairs, D/|A| is large enough that the sketch is no smaller than the KMC databases. That
  already happened for KOMPa vs KOMPb, at 42%, in KP-9. So the claim must be limited to near-duplicate pairs:
  re-assemblies, haplotype-resolved versions, or drafts vs finals.
- **Cost.** A C power-sum builder (not yet written; numpy took about 1 h per pair) is about 1-2 days. Then about 0.5
  day to run it on the existing panel.
- **Novelty (2026-11-11).** Queries: "power sum sketch set reconciliation k-mer sets genomes characteristic
  polynomial", "pangenome k-mer presence absence subtractable sketch exact set difference between assemblies IBLT",
  "colored de Bruijn graph color set difference sketch …". Found: Shibuya et al. (IBLT, already compared), Houen's
  "Simple Set Sketching" (power sums, general theory), MEMO and PanKmer (no sketches). **No genome or pangenome use of
  power-sum or characteristic-polynomial sketches was found.** This is not proof of absence. Houen 2022
  (arXiv:2211.03683) must be cited for the power-sum theory.

### 2. KmerStream's canonical hash is degenerate for odd k; fix it with min() or S

- **Claim.** KmerStream's `h ^ ht` (`RepHash.hpp:92-94`) maps every odd-k k-mer and its middle-complement twin to one
  value under every seed (10,000/10,000 at k = 31), so F0 and f1 are biased low wherever both twins occur. In reads
  that means sequencing errors that complement the middle base: about 1/(3k) of the error k-mers of a substitution,
  about 1% at k = 31, merge with true k-mers. It is also a published, real-world instance of the "obvious symmetric
  functions fail" result in our paper.
- **Refutation experiment.** Run KmerStream on one cannabis read set at k = 31 against KMC exact F0 and f1, then patch
  `hash()` to min(h, ht) or to S. It is refuted, as far as *mattering* goes, if the F0/f1 change is below KmerStream's
  own sampling error. The algebraic defect stays true either way.
- **Cost.** Under 0.5 day. Because the bias is probably small, the payoff is mainly for the paper: a cited example in
  the strand-symmetry section, and an upstream issue/PR if the user wants one.
- **Novelty (2026-11-11).** Query "KmerStream RepHash reverse complement XOR hash collision odd k middle base", plus
  KmerStream's GitHub issues #1-#16: nothing reports it. Last upstream commit 2016-06-14, so the repo is dormant.

### 3. FracMinHash × KP-10: a sampled, subtractable sketch for near-duplicate genomes (sourmash, Kmer-db)

- **Claim.** Take the kept hashes of a FracMinHash sketch at rate 1/s: **sourmash's existing Murmur values, so public
  databases stay compatible**, or our sampler's S values for new collections. Store their power sums. Two genomes'
  sampled symmetric difference then decodes exactly at about 11.3·D/s bytes, and sketches aggregate linearly: a
  clade's sum minus one member, or a metagenome minus decoded references. The exact sampled containment of
  near-identical genomes then comes from a sketch that scales with D/s, not |A|/s.
- **Where it can lose.** If both signatures are at hand, a sorted, delta-coded list of the differing hashes is
  smaller than 11.3 B per difference. The power sums win only when one side is absent (communication or
  reconciliation) or when linear aggregation is the point.
- **Refutation experiment.** Use the cannabis pairs from KP-10 (D = 7,931 and 3.2 M) at s = 100 and 1,000. Compare
  sketch bytes with a plain FracMinHash sketch (|A|/s hashes × 8 B), and check exactness. It is refuted if the
  per-bucket floor (2^b × (m+c) entries) makes the sketch no smaller than |A|/s × 8 B at panel-typical D. At
  |A| = 4.4×10⁸ and s = 1,000, a plain sketch is about 3.5 MB, against about 11.3 × D/1,000 B for the power sums.
- **Cost.** About 1 day once the C builder from idea 1 exists.
- **Novelty (2026-11-11).** Query "FracMinHash sketch subtract genomes invertible sketch recover differing k-mers
  near-identical assemblies". Found only Shibuya et al. (IBLT on MinHash-type sets) and the FracMinHash papers (no
  invertibility). **No prior FracMinHash power-sum combination was found.**

### 4. Long-k canonical keys without rotation-period collisions (Metagraph's Bloom for k > 64, Dashing 2 for k > 32)

- **Claim.** Cyclic-polynomial hashes over w-bit words have period w. Metagraph's Bloom prefilter (Lemire CyclicHash,
  64-bit, canonical = min) therefore collides deterministically for k > 64 whenever the bases at positions i and i+64
  are swapped: **1,000/1,000 at k = 65, 81 and 127** (agent test, `tmp/survey/cyc_test.cpp`). Dashing 2's k > 32 path
  uses the same family (bonsai `RollingHasher`, min of cyclic f and rc); whether it has the same period at its word
  size is UNCHECKED. S has no rotation period, is strand-symmetric for any k, and has a per-pair bound. Swapping it in
  removes this collision class at 0.46× batched-ntHash speed.
- **Refutation experiment.** Run a Metagraph k = 81 Bloom-prefiltered query on a repeat-rich genome (cs10) with and
  without S. Count Bloom false positives that reach the graph, and query time. It is refuted as *useful* if the
  swapped-pair collisions are negligible next to the filter's ordinary FPR. That is likely, because Bloom hits are
  verified, so the cost is only extra lookups. For Dashing 2, at k = 51, count distinct keys against an exact 128-bit
  canonical set; expect 0, or the 64-bit birthday level of about 10⁻².
- **Cost.** 1 day for both.
- **Novelty.** UNCHECKED as a search. The periodicity of cyclic polynomial hashing is textbook (Lemire & Kaser
  discuss its limited independence); what we add is the k > w consequence in these specific tools.

### 5. Feed ntCard-style estimators (ntCard, KMC's estimator, KmerGenie) from the kept-only sampler

- **Claim.** ntCard needs hash values only for the about 1/2¹⁰ of k-mers in its sampled sub-streams. A packed-code gate
  plus S for kept k-mers would run at about 1.2× batched ntHash, about 2.5× the library roll ntCard uses, with
  estimator accuracy unchanged.
- **Refutation experiment.** Run ntCard against the patched version on cs10 and one read set. It is refuted if the
  histogram differs beyond ntCard's own seed-to-seed variance, or if throughput doesn't improve end to end because
  parsing dominates.
- **Caveat that limits it.** The gate is a heuristic mixer, so its uniformity, not S, decides the estimate. (B) adds
  nothing to an estimator, so this is a speed engineering idea that borrows our sampler, not our algebra. Ranked last.
- **Novelty.** UNCHECKED (engineering).

## No fit, with the reason

- **Exact counters (Jellyfish, KMC 3, meryl, DSK; Kitsune via Jellyfish).** Keys are stored or sorted exactly. Either
  there is no hash (KMC, meryl), or a collision costs only a probe (DSK), or the hash is a bijection that also enables
  quotienting (Jellyfish). Lexicographic canonicalization on packed codes is exact and cheaper than S.
- **HLL / MinHash / count-sketch estimators (Dashing, Mash, KmerLight; ntCard for accuracy).** Their guarantees assume
  hash independence. S is not k-wise independent, and (B) concerns collisions, which these estimators don't need to
  avoid at 64 bits. At best S equals a mixer, at more than twice the cost.
- **minimap2 seeding.** The minimizer hash is a bijection, so it has no collisions, and minimizer density matched in
  our tests. Symmetric spaced seeds cost sensitivity (0.99^29 against 0.99^20) and minimap2 uses none.
- **LexicMap / LexicHash.** The algorithm depends on lexicographic (prefix) order, which XOR with a mask preserves and a
  field hash destroys. Strand symmetry is already exact, by scanning both strands. Spaced seeds would replace its
  variable-length prefix seeds rather than improve them. The k ≤ 32 limit comes from the 64-bit key and prefix order,
  so S can't lift it.
- **Replacing sourmash's hash with S or our sampler.** Every stored public signature is Murmur, seed 42, and its kept
  set can't be reproduced by a faster gate. Any speed gain costs re-sketching the databases, and no provable
  uniformity is gained (see the sourmash row). If wanted, run the speed idea only on private collections.
- **Pangenome indexes on MPHFs or ranks (Cuttlefish 2/3, Fulgor/SSHash, Themisto, SBWT, GGCAT's k-mer keys).**
  Keys are verified against stored sets, or there is no hash at all, so a collision bound buys nothing. GGCAT's SIMD
  ntHash minimizers are cheaper than S. Within one graph or index (Bifrost, Metagraph, Fulgor), columns can be diffed
  exactly, so sketches help only across separate builds.
- **O(1) substitution and concatenation.** These are shared with ntHash, and no surveyed tool uses them.

## Want-list (papers the library lacks)

Already in the library: KMC 3 (`kokot_2017_counting`), Mash (`ondov_2016_mash`), minimap2 (`li_2018_minimap2`), ntHash
(`mohamadi_2016_nthash`), ntHash2 (`kazemi_2022_nthash2`), meryl via Merqury (`rhie_2020_merqury`), PanKmer
(`aylward_2023_pankmer`), the k-mer review (`jenike_2025_kmer-approaches`) and knonical (`marcais_2024_knonical`).

| paper | DOI | open access |
|---|---|---|
| Marçais & Kingsford 2011, Jellyfish | 10.1093/bioinformatics/btr011 | PMC3051319 |
| Deorowicz et al. 2015, KMC 2 (signatures) | 10.1093/bioinformatics/btv022 | arXiv:1407.1507 |
| Rizk, Lavenier & Chikhi 2013, DSK | 10.1093/bioinformatics/btt020 | likely PMC, check |
| Mohamadi, Khan & Birol 2017, ntCard | 10.1093/bioinformatics/btw832 | likely PMC, check |
| Chikhi & Medvedev 2014, KmerGenie | 10.1093/bioinformatics/btt310 | likely PMC, check |
| Melsted & Halldórsson 2014, KmerStream | 10.1093/bioinformatics/btu713 | bioRxiv 10.1101/003962 |
| Sivadasan, Srinivasan & Goyal 2016, Kmerlight | arXiv:1609.05626 | arXiv |
| Deorowicz et al. 2019, Kmer-db | 10.1093/bioinformatics/bty610 | bioRxiv 10.1101/263590 |
| Baker & Langmead 2019, Dashing | 10.1186/s13059-019-1875-0 | Genome Biology (OA) |
| Baker & Langmead 2023, Dashing 2 | 10.1101/gr.277655.123 | PMC10538361 |
| Pornputtapong et al. 2020, KITSUNE | 10.3389/fbioe.2020.556413 | Frontiers (OA) |
| Shen, Lees & Iqbal, LexicMap | Nat. Biotechnol. (2025; DOI to confirm) | bioRxiv 10.1101/2024.08.30.610459 |
| Greenberg et al., LexicHash | UNCHECKED, locate | — |
| Ertl 2021, SetSketch (Dashing 2's sketch) | 10.14778/3476249.3476276 | arXiv:2101.00314 |
| Pătraşcu & Thorup 2012, simple tabulation (for the independence caveat) | 10.1145/2220357.2220361 | arXiv:1011.5200 |
| Dodis et al. 2008, fuzzy extractors (PinSketch) | 10.1137/060651380 | cs/0602007 |
| Houen 2022, Simple Set Sketching (power sums) | arXiv:2211.03683 | arXiv |
| Brown & Irber 2016, sourmash | 10.21105/joss.00027 | JOSS (OA) |
| Irber et al. 2022, FracMinHash and gather | 10.1101/2022.01.11.475838 | bioRxiv |
| Hera, Pierce-Ward & Koslicki 2023, FracMinHash confidence intervals and ANI (for the uniformity assumption) | 10.1101/gr.277651.123 | Genome Research (OA) |

Pangenome tools. Panagram (`jenike_2026_panagram`) and PanKmer are already in the library. DOIs below were checked
against CrossRef on 2026-11-11 unless marked.

| paper | DOI | open access |
|---|---|---|
| Lemane et al. 2022, kmtricks | 10.1093/bioadv/vbac029 | OA |
| Holley & Melsted 2020, Bifrost | 10.1186/s13059-020-02135-8 | OA |
| Cracco & Tomescu 2023, GGCAT | 10.1101/gr.277615.122 | OA |
| Khan et al. 2022, Cuttlefish 2 | 10.1186/s13059-022-02743-6 | OA |
| Karasikov et al., MetaGraph | bioRxiv 10.1101/2020.10.01.322164 (Nature 2025 DOI unverified) | bioRxiv |
| Alanko et al. 2023, Themisto | 10.1093/bioinformatics/btad233 | OA (ISMB) |
| Fan et al. 2024, Fulgor | 10.1186/s13015-024-00251-9 | OA |
| Alanko, Puglisi & Vuohtoniemi 2023, SBWT | 10.1137/1.9781611977714.20 | arXiv:2206.01688 |
| Pibiri 2022, SSHash | 10.1093/bioinformatics/btac245 | OA |
| Pibiri & Trani 2021, PTHash | 10.1145/3404835.3462849 | arXiv |
| Limasset et al. 2017, BBHash | 10.4230/LIPIcs.SEA.2017.25 (CrossRef lookup failed; verify) | LIPIcs (OA) |

DOIs were checked against CrossRef on 2026-11-11, except LexicMap (journal DOI not confirmed) and LexicHash (not
located). PMC10538361 (Dashing 2) appeared in a search result; the other PMC and arXiv identifiers are from memory
(NCBI and arXiv APIs were unreachable from the sandbox) and need checking at filing.
