# Gene copies by k-mers: the cannabinoid synthase pilot (2026-11-11)

**Verdict: on this family a k-mer scan does not earn its place against lastz and miniprot.** It finds every copy at
80% DNA identity or more, the whole synthase clade, with the right orientation and no false loci. But lastz finds that
too and also the 35 deeper relatives, and miniprot adds frameshift and stop calls. Both cost less CPU than our numpy
scan (about 24–33 CPU-s per genome, against 84–260 CPU-s per seed). The k-mer scan is worth keeping only as a
cheap presence screen when it runs compiled, and only for copies above about 80%.

Code: `ideas/gene_copies/gene_copies.py` (subcommands `verify-refs stage kmer lastz miniprot blast truth score census
spike self-test`). The references and their manifest are in `ideas/gene_copies/refs/`. Small outputs are in
`results/gene_copies_synthase_2026-11-11/`. Raw outputs are in `bench/data/gene_copies/<genome>/`.

## The premise: no exons to separate

THCAS and CBDAS are intact single-exon genes. AB057805.1's CDS is one interval (27..1664), and Taura 2007 got
CBDAS-type genes from genomic DNA as "a 1635-bp intronless open reading frame" (library summary, lines 275–277). In
cs10's RefSeq annotation all four synthase-named models (LOC115697880, LOC115696909, LOC115696884, LOC115697762)
have one exon and one CDS segment. So do 30 of the 35 oxidocyclase/BBE-like models.

So "exons separately" and the intron-retained versus processed-pseudogene split do not apply. In their place:

- **Full or partial along the reference.** A copy is full if its alignment or k-mer chain spans at least 90% of the
  reference. Partial copies are reported with their reference interval.
- **Intact or pseudogene.** The k-mer scan, lastz and the truth set use an ORF test: the locus is extended to both
  reference ends plus 300 bp, and it is intact if its longest ATG..stop ORF is 95–115% of the reference protein.
  miniprot uses its own calls: intact means full, with 0 frameshifts, 0 in-frame stops and no "intron" (in an
  intronless gene, an intron call means an insertion).
- **Split copies.** Segments that are collinear along the reference and less than 15 kb apart are joined into one
  copy, so a copy interrupted by an LTR insertion is counted once.
- **Domain segments, not done.** The ~1.6-kb genes are scored as one unit. Splitting them into FAD-binding and BBE
  domain segments would mean placing domain boundaries on the protein, which was not done here. Partial copies already
  report which part of the reference they cover.

## References (all verified against NCBI, `verify-refs`)

| gene | accession | CDS | protein | note |
|---|---|---|---|---|
| THCAS | AB057805.1 (mRNA) | 27..1664, 1,638 bp with TAA | BAC41356.1 | Sirikantaramas 2004 |
| CBDAS | AB292682.1 (mRNA) | 1..1635 with TAA | BAF65033.1 | Taura 2007 |
| CBCAS | LY658671.1 (patent KR 1020190025485-A/5) | 1..1635, no stop | its translation, 545 aa | organism field "synthetic construct"; same reference as the hub track |

`verify-refs` re-fetched each record. Each CDS matches byte for byte, each translation matches its NCBI protein, and
the SHA-256 values are in `refs/manifest.json`. The CDS file is the same as cannabis-genome's
`data/reference/cannabinoid_synthases.fa`.

## Truth set, cs10: 48 loci

The truth set is the union, merged by overlap, of four sources:

- Grassa 2021 Table S9: 13 copies, from `grassa_s9_refseq_overlap.tsv`.
- RefSeq's 35 protein-coding oxidocyclase/BBE-like genes, from `cs10_tandem_genes.tsv`.
- blastn (dc-megablast) loci.
- tblastn loci.

Each source needed a reference cover of at least 150 bp. blastn and tblastn are used only to build the truth; they
are not among the methods scored.

Each locus's identity is the blastn identity (summed over its HSPs) to its closest reference, the one with the most
identical aligned bases. Status comes from the ORF test.

| bin (DNA identity) | loci | what they are |
|---|---|---|
| ≥ 95% | 7 | the 26 Mb CBCAS array: the complete CBCAS at 25.82 Mb (99.76% to CBCAS; RefSeq calls it "inactive THCAS"), 5 CBCAS pseudogenes, and CBDAS at 30.98 Mb (100%) |
| 90–95% | 6 | the 29 Mb CBDAS-like array (5 copies, CBDAS-like 1 the only intact ORF) and the 26 Mb THCAS-like pseudogene (Grassa's label; 91.5% to CBCAS here) |
| 80–90% | 0 | none in cs10 |
| < 80% | 35 | BBE-like relatives: 6 at 70–80%, 19 under 70%, 10 with no DNA hit (tblastn or RefSeq only) |

The synthase clade (≥ 90%) is exactly Grassa's 13 copies: 3 intact (CBCAS 25.82 Mb, CBDAS-like 1, CBDAS), 8 full
pseudogenes and 2 partial copies. Grassa says "all but two homologs were pseudogenes". The third intact ORF here is
CBDAS-like 1, which Grassa also calls a full ORF. miniprot's frameshift and stop calls agree with the ORF test on all
12 clade loci it reports separately. The truth set also holds 4 loci from neither annotation source: 3 tblastn-only
fragments (chr1 29.6 Mb, chr3 32.4 Mb, chr7 30.87 Mb) and 1 pseudogene record (chr7 49.33 Mb).

## Recall and precision on cs10

All methods searched the whole cs10 (876 Mb, 221 sequences). A k-mer locus is called when its hits on one diagonal
band (±150 bp, hits less than 1 kb apart) cover at least 150 reference bases. This threshold was fixed before
scoring.

The k-mer seeds are canonical (strand-symmetric) contiguous k-mers at k = 11, 15, 21 and 31. There are two more
variants:

- **k15d1:** 15-mers with every 1-mismatch neighbour of each reference 15-mer added. This is the radius-1 Hamming ball
  B_1 of KP-11, enumerated directly rather than by transform.
- **sp31w24:** the palindromic spaced seed 1110111011101110111011101110111 from ideas/algebra and ideas/eisenstein
  (KP-14). Its packed key min(F, R) is the twin of spaced_e.h's (F+λ)(R+λ), and the self-test checks strand symmetry.

Recall per identity bin (found / truth):

| method | loci | precision | ≥95 | 90–95 | 80–90 | <80 (70–80 / <70 / protein-only) | orientation | intact/pseudo agree* | CPU-s |
|---|---|---|---|---|---|---|---|---|---|
| k11 | 17 | 1.00 | 7/7 | 6/6 | – | 4/35 (3/6, 1/19, 0/10) | 17/17 | 13/13 | 112 |
| k15 | 15 | 1.00 | 7/7 | 6/6 | – | 2/35 (2/6, 0/19, 0/10) | 15/15 | 13/13 | 145 |
| k21 | 14 | 1.00 | 7/7 | 6/6 | – | 1/35 | 14/14 | 13/13 | 182 |
| k31 | 13 | 1.00 | 7/7 | 6/6 | – | 0/35 | 13/13 | 13/13 (full 9/13) | 229 |
| k15d1 | 19 | 1.00 | 7/7 | 6/6 | – | 6/35 (3/6, 3/19, 0/10) | 19/19 | 13/13 | 168 |
| sp31w24 | 14 | 1.00 | 7/7 | 6/6 | – | 1/35 | 14/14 | 13/13 | 194 |
| **lastz** | 48 | 1.00 | 7/7 | 6/6 | – | **35/35** | 48/48 | 13/13 | **26** |
| **miniprot** | 42 | 1.00 | 7/7 | 5/6† | – | 30/35 (6/6, 19/19, 5/10) | 42/42 | 12/12 | **27** |

\* The ORF test is the same classifier for the k-mer methods, lastz and the truth, so their agreement is partly
built in. miniprot's agreement on 12 of 12 is the independent check.

† miniprot spliced the partial copy at 29.545 Mb onto its neighbour at 29.576 Mb, 31 kb away, as one gene with an
"intron". One-to-one matching credits it with one copy.

Caveats on these numbers:

- **Precision of 1.00 is partly by construction.** The truth set includes every blastn and tblastn locus, so a false
  locus would have to evade both of them. No method called a locus outside the truth set.
- **CPU and wall time.** CPU-s is user plus system time over all workers. The k-mer runs used 6 processes and took
  30–50 s of wall time. lastz ran single-threaded on the whole genome, so its CPU-s is the honest figure; its wall
  time was inflated by 5 concurrent runs. miniprot ran with 16 threads on 150-Mb parts, because the whole-genome index
  does not fit in this sandbox's 15 GB. blastn's makeblastdb took 516 s; the searches took seconds.

## Spike-in: the 80–90% bin cs10 leaves empty

`spike` planted 63 copies into cs10 chr1 0–20 Mb, which holds no family member: 3 references × 7 identities (97, 92,
87, 82, 77, 72, 67%) × 3 replicates, on random strands. The mutations are 85% substitutions (61% of them transitions)
and 15% indels of 1–3 bp.

| method | ≥95 | 90–95 | 80–90 | 70–80 | <70 |
|---|---|---|---|---|---|
| k11 | 9/9 | 9/9 | 18/18 | 17/18 | 2/9 |
| k15 | 9/9 | 9/9 | 18/18 | 3/18 | 0/9 |
| k21 | 9/9 | 9/9 | 11/18 | 0/18 | 0/9 |
| k31 | 9/9 | 9/9 | 2/18 | 0/18 | 0/9 |
| k15d1 | 9/9 | 9/9 | 18/18 | 18/18 | 6/9 |
| sp31w24 | 9/9 | 9/9 | 12/18 | 0/18 | 0/9 |
| lastz | 9/9 | 9/9 | 18/18 | 18/18 | 9/9 |
| miniprot | 9/9 | 9/9 | 18/18 | 9/18 | 0/9 |

Precision was 1.00 for every method, and every orientation was right. Two points stand out:

- **Short or mismatch-tolerant seeds are required below 90%.** At 80–90%, k ≥ 21 and the weight-24 spaced seed lose
  copies (k21 11/18, sp31w24 12/18, k31 2/18). k11 and k15d1 keep all of them.
- **The spike understates miniprot.** Its substitutions are random, so they change amino acids far more often than
  real paralogs, which are under purifying selection. On real BBE-like genes at 70–80% DNA identity, miniprot found
  6/6.

## Panel: EH23a, EH23b and Finola

All methods ran on each assembly. The synthase clade is lastz loci at ≥ 80% identity. Every k-mer seed except k31
found every clade locus or all but one. k11 and k15d1 found all of them. k15, k21 and sp31w24 each missed one
228-bp, ~90% fragment in EH23b and one in Finola: below the 150-bp cover threshold. The intact counts are the same across methods to within 1.

| genome | clade copies (lastz) | ≥98% full intact (lastz + ORF) | ≥98% pseudogene | hub track "Full" (Lynch rule) | outside clade (lastz / k15d1) |
|---|---|---|---|---|---|
| cs10 | 13 (chr7) | CBDAS 1, CBCAS 1 | CBCAS 4 (1 partial) | CBDAS 1, CBCAS 1 | 35 / 6 |
| EH23a | 14 (chr7) | **THCAS 1**, CBCAS 3 | 0 | THCAS 1, CBCAS 3 | 33 / 7 |
| EH23b | 15 (chr7) | **CBDAS 1**, CBCAS 4 | CBCAS 3 (1 partial) | CBDAS 1, CBCAS 6 | 34 / 6 |
| Finola | 17 (7 contigs) | **CBDAS 1**, CBCAS 2 | CBCAS 5 | not on the hub | 43 / 9 |

What the panel shows:

- **Chemotype.** Each genome carries one functional THCAS or CBDAS, never both (EH23a THC, EH23b CBD, Finola CBD with
  no THCAS-identical copy). This agrees with Laverty 2019 for Finola and Lynch 2025's "maximum of one".
- **Copy-number differences** are in CBCAS and the pseudogene clades:
  - Near-identical CBCAS copies: 3 in EH23a, 7 in EH23b, 7 in Finola, 5 in cs10.
  - Copies at 80–98% to their closest reference, mostly the THCAS-like/CBDAS-like pseudogene clades: 7–10 per
    genome.
- **The hub's Lynch "Full" rule counts pseudogenes.** That rule is ≥ 500 bp, < 10 mismatches and no gaps. EH23b's
  6 Full CBCAS calls include a 1,166-bp partial copy (28.887 Mb) and a full-length copy at 99.63% that carries 2
  in-frame stop codons (28.933 Mb; ORF 303 aa, miniprot st=2). A third near-identical copy (28.749 Mb, 99.82%) has a
  frameshift, so the gapless rule drops it. With the ORF test, EH23b has 4 intact CBCAS, not 6. EH23a's 3 match
  exactly. This is a finding for cannabis-genome's synthase track; nothing there was changed.

## Does a k-mer method earn its place?

**Recall per identity bin.**

- At **≥ 95% and 90–95%**, every method has 100% recall, except miniprot's one splice merge.
- At **80–90%** (spike only), k11, k15, k15d1, lastz and miniprot have 100%. k21 has 61%, the spaced seed 67% and
  k31 11%.
- At **70–80%**, lastz and k15d1 have 100%. k11 is close (94% spike, 3/6 real), miniprot gets 50% on the spike but
  6/6 on real genes, and k15 has 17%.
- **Under 70%**, lastz has 100%, miniprot 19/19 on real genes (0/9 on the random spike), k15d1 6/9 on the spike
  but 3/19 on real genes, and the other seeds ≤ 2/9 and ≤ 1/19.
- **Protein-only relatives** (no DNA hit): lastz 10/10, miniprot 5/10, every k-mer seed 0/10.

**Speed.**

- Whole-genome lastz: about 24–33 CPU-s per assembly.
- miniprot: about 24–30 CPU-s.
- The numpy k-mer scan: 84–260 CPU-s per seed. A compiled roller (KP-7e) would cut that, but that was not tested here.
- Even at zero cost, the scan only matches lastz's recall above 80%, and alone it cannot call frameshifts.

**Where it could still help.**

- As a presence screen across hundreds of assemblies: k15d1 recovers the whole clade with correct orientation.
- In alignment-free classification: the ORF test on k-mer loci agreed with miniprot on every clade copy.

For copy discovery, lastz (DNA) plus miniprot (protein) is the better pair.

## Not done

- No domain-segment split (FAD and BBE boundaries).
- No second-family test.
- The KP-11 transform itself was not used: the d = 1 ball was enumerated directly, which is cheap for one family.
- Finola's clade sits on 7 contigs, so its tandem order is not resolved.
