# Literature for IDEAS.md section 6

**Status: ACQUIRED, 2026-11-10. 24 papers found: 18 open-access PDFs downloaded, and the 6 below fetched by the owner the same day. All 24 are in `bench/data/papers/kmer_poly_oa/`, each checked by first-page text. None was already in the library.**

Drop fetched PDFs in `bench/data/papers/kmer_poly_oa/`, named as in the key column. This is acquisition only; nothing is
filed in the library yet.

## Bearing on the novelty question (KP-9 / KP-10)

- **Reconciliation of genomic k-mer sets is published, with IBLTs.** Shibuya, Belazzougui and Kucherov (WABI 2022)
  apply invertible Bloom lookup tables, with syncmer sampling, to k-mer sets from similar genomes. They estimate Jaccard
  similarity, and recover the differing k-mers, in space proportional to the difference. They tested it on SARS-CoV-2
  and *S. pneumoniae*. Only two papers cite it (Hyper-k-mers 2024, and an IBLT 2-core paper at SEA 2026).
- **Characteristic polynomials and power sums for k-mer sets: not found.** No OpenAlex or web search found the
  characteristic-polynomial (CPI, Minsky et al.) or power-sum / BCH-syndrome (PinSketch, Minisketch) form applied to
  k-mers or sequencing data. Both forms are well developed outside genomics: Minisketch is deployed in Bitcoin
  (Erlay), and Boral and Mitzenmacher compare CPI with IBLTs.
- **What a KP-9/KP-10 paper would have to show** is what the algebraic form gives over Shibuya's IBLTs:
  - no decode-failure probability
  - linear, mergeable sketches (KP-10)
  - strand symmetry built into the hash (KP-7)
  - the partitioned scaling of Minsky and Trachtenberg (2002)

  Our prototype's measured haplotype difference also matters: Shibuya's method, like ours, pays off only when the sets
  are near-identical.
- **Strand handling (KP-7):** Marçais, Elder and Kingsford (2024) show that canonical k-mers create "sketching deserts"
  in context-free sketches. That argues for a strand-symmetric hash rather than min(fwd, rc).

## Fetched by the owner (paywalled, or the download failed), received 2026-11-10

| priority | key | paper | DOI | KP | why |
|---|---|---|---|---|---|
| 1 | minsky_2003_reconciliation | Minsky, Trachtenberg & Zippel 2003, Set reconciliation with nearly optimal communication complexity, IEEE Trans. Inf. Theory | 10.1109/tit.2003.815784 | KP-8, KP-9 | the CPI algorithm KP-9 implements. The ISIT 2001 version (10.1109/isit.2001.936095) would do if it is easier to get |
| 1 | kazemi_2022_nthash2 | Kazemi et al. 2022, ntHash2: recursive spaced seed hashing for nucleotide sequences, Bioinformatics | 10.1093/bioinformatics/btac564 | KP-6, KP-7, KP-14 | the rolling canonical hash to beat, with spaced seeds. It is open access, but every download returned the supplement or failed |
| 2 | anastassiou_2001_gsp | Anastassiou 2001, Genomic signal processing, IEEE Signal Processing Magazine | 10.1109/79.939833 | 0, KP-4 | the source usually cited for the complex encoding (the 2000 Bioinformatics paper is downloaded) |
| 2 | chen_2026_hashreview | Chen et al. 2026, Hash functions in nucleotide sequence analysis, Genome Research | 10.1101/gr.281453.125 | KP-6 | the current survey of nucleotide hashing. It is open access, but the download was blocked |
| 3 | voss_1992_fractal | Voss 1992, Evolution of long-range fractal correlations and 1/f noise in DNA base sequences, Phys. Rev. Lett. | 10.1103/physrevlett.68.3805 | 0, KP-1 | the four indicator sequences that KP-1 needs for exact Hamming distance |
| 3 | das_2021_encoding | Das, Das, Mohapatra & Nanda 2021, DNA numerical encoding schemes for exon prediction: a recent history, Nucleosides Nucleotides Nucleic Acids | 10.1080/15257770.2021.1966797 | 0 | recent review of encodings, to check whether our exact map (A=1, C=i, G=−i, T=−1) has a name |

## Open access, downloaded to `bench/data/papers/kmer_poly_oa/`

| key | paper | DOI | KP | note |
|---|---|---|---|---|
| shibuya_2022_reconciliation | Shibuya, Belazzougui & Kucherov 2022, Efficient reconciliation of genomic datasets of high similarity, WABI | 10.4230/lipics.wabi.2022.14 | KP-9, KP-10 | **the closest prior work** |
| minsky_2002_practical | Minsky & Trachtenberg 2002, Practical set reconciliation, Boston Univ. tech. report (as "Scalable set reconciliation", Allerton 2002) | none | KP-9 | the partitioned protocol. The text layer is unreadable, so it needs OCR before filing |
| boral_2014_multiparty | Boral & Mitzenmacher 2014, Multi-party set reconciliation using characteristic polynomials | 10.48550/arxiv.1410.2645 | KP-8, KP-9 | CPI against IBLT |
| mitzenmacher_2017_multiparty | Mitzenmacher & Pagh 2017, Simple multi-party set reconciliation, Distributed Computing | 10.1007/s00446-017-0316-0 | KP-9, KP-10 | arXiv version |
| dodis_2008_fuzzy | Dodis, Ostrovsky, Reyzin & Smith 2008, Fuzzy extractors, SIAM J. Comput. | 10.1137/060651380 | KP-10 | PinSketch (sec. 6.3), the power-sum BCH sketch. arXiv version |
| naumenko_2019_erlay | Naumenko, Maxwell, Wuille et al. 2019, Erlay / bandwidth-efficient transaction relay in Bitcoin, ACM CCS | 10.1145/3319535.3354237 | KP-10 | Minisketch power sums in production. arXiv version |
| goodrich_2011_iblt | Goodrich & Mitzenmacher 2011, Invertible Bloom lookup tables, Allerton | 10.1109/allerton.2011.6120248 | KP-10 | arXiv version |
| eppstein_2011_difference | Eppstein, Goodrich, Uyeda & Varghese 2011, What's the difference? Efficient set reconciliation without prior context, SIGCOMM | 10.1145/2018436.2018462 | KP-9, KP-10 | author copy (ics.uci.edu) |
| mohamadi_2016_nthash | Mohamadi et al. 2016, ntHash: recursive nucleotide hashing, Bioinformatics | 10.1093/bioinformatics/btw397 | KP-6, KP-7 | |
| benbassat_2014_incremental | Ben-Bassat & Chor 2014, String graph construction using incremental hashing, Bioinformatics | 10.1093/bioinformatics/btu578 | KP-6 | rolling polynomial (Karp-Rabin) hashing in assembly |
| marcais_2024_knonical | Marçais, Elder & Kingsford 2024, k-nonical space: sketching with reverse complements, Bioinformatics | 10.1093/bioinformatics/btae629 | KP-7 | **bioRxiv preprint** (10.1101/2024.01.25.577301); the journal PDF did not download |
| martayan_2024_cbl | Martayan et al. 2024, Conway–Bromage–Lyndon: an exact, dynamic representation of k-mer sets, Bioinformatics | 10.1093/bioinformatics/btae217 | KP-8 | set-operation baseline |
| sladky_2024_fmsi | Sladký, Veselý & Břinda 2024, Towards efficient k-mer set operations via function-assigned masked superstrings | 10.1101/2024.03.06.583483 | KP-8, KP-11 | **preprint only** (bioRxiv) |
| anastassiou_2000_frequency | Anastassiou 2000, Frequency-domain analysis of biomolecular sequences, Bioinformatics | 10.1093/bioinformatics/16.12.1073 | 0, KP-4 | |
| mendizabal_2017_numerical | Mendizabal-Ruiz et al. 2017, On DNA numerical representations for genomic similarity computation, PLoS ONE | 10.1371/journal.pone.0173288 | 0, KP-1 | compares encodings, the complex one included |
| wang_2009_mapping | Wang & Schonfeld 2009, Mapping equivalence for symbolic sequences: theory and applications, IEEE Trans. Signal Processing | 10.1109/tsp.2009.2026544 | 0 | when two numeric maps are equivalent. arXiv version |
| faloutsos_1994_subsequence | Faloutsos, Ranganathan & Manolopoulos 1994, Fast subsequence matching in time-series databases, SIGMOD | 10.1145/191843.191925 | KP-4 | DFT lower bounds (GEMINI). Author copy (cmu.edu) |
| katoh_2002_mafft | Katoh et al. 2002, MAFFT: rapid multiple sequence alignment based on fast Fourier transform, NAR | 10.1093/nar/gkf436 | KP-13 | the FFT step |

## Already in the library

None of the 24 DOIs is in `library/src/biblio/papers.json` (913 DOIs checked on 2026-11-10).

## Searches that found nothing relevant (evidence for the novelty question)

OpenAlex, 2026-11-10:

- "set reconciliation genome sequencing": only assembly and phylogenetic "reconciliation", a different sense of the word
- "characteristic polynomial k-mer set": nothing algebraic (matchtigs, gapped k-mers)
- "power sums sketch k-mers": nothing algebraic (MinHash-family sketches only)
- "Gaussian integers DNA sequence representation": nothing
- "polynomial rolling hash reverse complement DNA": ntHash, Ben-Bassat & Chor, and the 2026 hash review only

Web search, same day:

- "set reconciliation k-mer sets characteristic polynomial OR power sums genomes sequencing": only Shibuya et al.
  (IBLT) and generic reconciliation papers
- no source using the exact map A=1, C=i, G=−i, T=−1 with complement as negation. Reviews list "complex" encodings
  generically, so check Das et al. 2021 and Anastassiou 2001 for the exact assignment
