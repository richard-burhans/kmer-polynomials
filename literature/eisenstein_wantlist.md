# Literature for the Eisenstein encoding: novelty check (2026-11-10)

**Status (updated 2026-11-11):** all six owner-fetched papers received (bench/data/papers/eisenstein_oa/). The owner supplied Das & Das 2017 (I2CT, Kaiser window) in place of Das et al. 2018 (Genomics); the latter is no longer needed, see the resolution below.

**Original status: OPEN. Acquisition only, not yet filed in the library.** The encoding under check: A = 1, C = ω, G = −ω,
T = −1 (ω a primitive cube root of unity). Complement is negation. In ℂ, transitions have squared distance 1,
A–C/G–T transversions 3, and complements 4. The same values serve as polynomial-hash coefficients over F_q with ω ∈ F_q
(hash/kp7e.*, results/kp7_eisenstein_2026-11-10.md). OA PDFs are in `bench/data/papers/eisenstein_oa/`.
None of the papers below was already in the library (checked by DOI against src/biblio/papers.json and refs/).

## Assessment

1. **The map itself is a special case of a published family.** Das, Nanda & Das's "trigonometric" mapping (2017,
   as restated by das_2021_dna-encoding) is A = −cosθ + j sinθ, C = −cosθ − j sinθ, G = cosθ + j sinθ,
   T = cosθ − j sinθ, so T = −A and G = −C. Our fourth-root encoding was θ = π/4. **Our Eisenstein encoding is θ = π/3,
   up to a rotation** (multiplying by ω̄ sends A → 1, C → ω, G → −ω, T → −1). Das 2021 never fixes θ.
   - Gulocak et al. 2022 (SAUCIS, downloaded) give the trigonometric mapping with a **θ = 60° worked example**
     (0.5 + 0.866i). Their sign table, though, makes A and T conjugates rather than negatives, and their example vector
     is internally inconsistent.
   - So **θ = π/3 within this family has appeared in print**, though not clearly with complement = negation.
   - **Open:** which θ the original 2017 paper and Das et al. 2018 (Genomics) used. Both are on the fetch list.
2. **Choosing θ = π/3 BECAUSE it makes transitions the cheapest substitution: not found.** No paper found analyses
   the trigonometric family, or any complex mapping, in terms of transition vs transversion distances. Our criterion is
   that transitions separate from A–C/G–T transversions iff Re(A·C̄) ≠ 0, and fourth roots are exactly the orthogonal
   case. This looks unpublished, but the trigonometric papers' own θ discussions must be read before claiming it.
3. **Using the map as polynomial-hash coefficients over a prime field with ω ∈ F_q (q ≡ 1 mod 3), with a
   strand-symmetric (F+λ)(R+λ): not found.** No hashing paper uses roots-of-unity or Eisenstein coefficients.
4. **Related but DIFFERENT algebra:** DNA error-correcting codes over the "Eisenstein field" Z₂[ω] = F₄ map
   {0, 1, ω, 1+ω} → {A, C, G, T}, with complement as x ↦ x + 1, in characteristic 2. That has no negation, no complex
   distances and no hashing, so it is not the same construction. It is the long-standing GF(4)/ring DNA-code literature.

**Bottom line:** the encoding should be credited as the θ = π/3 member of Das et al.'s trigonometric mapping. What
looks new is (a) the transition/transversion criterion for choosing θ, and (b) its use as F_q polynomial-hash
coefficients with an exact strand-symmetric product. Both are pending the fetch list below.

## Please fetch (paywalled, or the download was blocked), most relevant first

| priority | paper | DOI | why |
|---|---|---|---|
| 1 | Das, Nanda & Das 2017, A novel DNA mapping scheme for improved exon prediction using digital filters, MAMI 2017 (IEEE) | 10.1109/mami.2017.8307889 | the original trigonometric mapping: which θ, which sign convention, any rationale for θ |
| 1 | Das, Nanda & Das 2018, An integrated approach for identification of exon locations using recursive Gauss Newton tuned adaptive Kaiser window, Genomics | 10.1016/j.ygeno.2018.10.008 | uses the trigonometric mapping, probably with a fixed θ. OA at ScienceDirect, but the download was blocked |
| 2 | 2-Simplex mapping for identifying the protein coding regions in DNA, TENCON 2007 (IEEE) | 10.1109/tencon.2007.4429086 | a "2-simplex" (triangle) mapping may use 120° geometry; check whether it is Eisenstein-like |
| 2 | Roy & Barman 2016, Polyphase filtering with variable mapping rule in protein coding region prediction, Microsyst. Technol. | 10.1007/s00542-016-2884-5 | complex "twiddle factor" mappings varying by codon position; check the angles used |
| 3 | Cristea 2002, Conversion of nucleotides sequences into genomic signals, J. Cell. Mol. Med. | 10.1111/j.1582-4934.2002.tb00196.x | the classic complex (quadrantal) representation, which discusses purine/pyrimidine structure. OA (PMC6740102), download blocked |
| 3 | Yu, Li & Zeng 2018, Survey on encoding schemes for genomic data representation and feature learning, Big Data Mining and Analytics | 10.26599/bdma.2018.9020018 | survey of mappings, to check for any 120°/Eisenstein entry. OA, download blocked |

## Open-access PDFs downloaded (bench/data/papers/eisenstein_oa/)

| file | paper | DOI |
|---|---|---|
| saucis_2022_mapping-effect.pdf | Gulocak et al. 2022, The effect of numerical mapping techniques on performance in genomic research, SAUCIS | 10.35377/saucis...1191850 |
| frsip_2025_static-encoding.pdf | 2025, A computational approach for prediction of exons using static encoding methods, digital filter and wavelet…, Front. Signal Process. | 10.3389/frsip.2025.1679555 |
| yu_2016_gene-identification-review.pdf | Yu et al. 2016, A comprehensive review of emerging computational methods for gene identification, JIPS (image-only PDF, no text layer: OCR before filing) | 10.3745/jips.04.0023 |
| eisenstein_2026_dna-bch.pdf | DNA sequences from BCH codes over Eisenstein field and their decoding algorithm, AIMS Math. 2026 | 10.3934/math.20261378 |
| eisenstein_2023_bch.pdf | Eisenstein field BCH codes construction and decoding, AIMS Math. 2023 | 10.3934/math.20231508 |

## Queries that found nothing relevant (OpenAlex, 2026-11-10)

- "hexagonal DNA numerical representation", "cube root of unity DNA mapping", "sixth roots of unity DNA"
- "DNA numerical representation transitions transversions complex", "DNA numerical mapping purine pyrimidine distance transition transversion representation"
- "rolling hash roots of unity", "polynomial hash complex coefficients DNA k-mer", "Eisenstein integer hash k-mer"
- "k-mer hashing finite field reverse complement symmetric", "trigonometric mapping DNA theta 60 degrees"

The Eisenstein-integer hits were all image cryptography or GF(4)/ring DNA codes (item 4 above).


## Resolution (2026-11-11): the exact map IS published

Das, Nanda & Das 2017 (MAMI, das_2017_trigonometric-mami.pdf, eq. 2) define A = −cos θ + j sin θ, C = −cos θ − j sin θ,
G = cos θ + j sin θ, T = cos θ − j sin θ, and state "the θ value is taken as pi/3". This has complement = negation and
exactly our squared distances: 1 for A–G and C–T (transitions), 3 for A–C and G–T, and 4 for complements. It is our
Eisenstein map up to a rotation. Their stated reasons are that purines and pyrimidines share imaginary parts and that
complements sum to zero. They do NOT mention transitions or distances.

What remains unclaimed, as far as these searches go:
1. the explicit criterion (transitions separate from A–C/G–T exactly when Re(A·C̄) ≠ 0, so no fourth-root map can) and
   transition-cheapest as the REASON for θ = π/3;
2. the use as polynomial-hash coefficients over F_q with ω ∈ F_q, and the strand-symmetric product (F+λ)(R+λ).

The other papers: Grandhi & Kumar 2007's 2-simplex map puts three bases on a triangle and the fourth at 0 (different).
Roy & Barman 2017 use a codon-level map with angle π/8. Cristea 2002 and Yu et al. 2018 are background.

## Filed in the library (2026-11-11)

All 11 papers are filed on `methods-tooling` (library commits 82664ac and 4ec1833). The slugs are
das_2017_trigonometric-mapping, das_2017_kaiser-window, cristea_2002_genomic-signals, grandhi_2007_2simplex-mapping,
roy_2017_polyphase-mapping, yu_2018_encoding-survey, gulocak_2022_mapping-effect, benarjee_2025_static-encoding,
yu_2016_gene-identification (OCR re-render), sajjad_2023_eisenstein-bch and sajjad_2026_eisenstein-dna-bch. Second
reads are owed.

Filing found two more things:
- `cristea_2002_genomic-signals` equation (5) is the fourth-root map itself (up to a global sign), and its equation
  (4) is the θ = π/4 trigonometric point. Das, Nanda & Das 2017's stated rationale restates Cristea's.
- `gulocak_2022_mapping-effect`'s registered DOI is 10.35377/saucis...1191850 (three dots). Its trigonometric row
  swaps A and G relative to Das, and its θ = 60° example vector contradicts that row at every C and T.
