---
parent: How It Works
title: Gene Definition Exceptions
permalink: methods/Gene-Definition-Exceptions/
nav_order: 4
---
# Gene Definition Exceptions

Genotype determination is based on [CPIC gene definition tables](https://www.pharmgkb.org/page/pgxGeneRef), with 
modifications for the following genes:

* [CYP3A4](#cyp3a4)
* [DPYD](#dpyd)
* [G6PD](#g6pd)
* [NAT2](#nat2)
* [SLCO1B1](#slco1b1)
* [RYR1](#ryr1)

---    


## CYP3A4

PharmGKB annotates PGx-based drug dosing guidelines published by the [Royal Dutch Association for the Advancement of
Pharmacy - Pharmacogenetics Working Group (DPWG)](https://www.pharmgkb.org/page/dpwg). PharmGKB curates allele function assignments and phenotype
mappings from the DPWG to provide genotype-specific DPWG guideline recommendations. Where possible, PharmGKB maps DPWG
terms to CPIC terms, as outlined on [PharmGKB](https://www.pharmgkb.org/page/dpwgMapping).

CYP3A4 is currently not part of a CPIC guideline. Since the DPWG CYP3A4 documentation includes limit variant notations
for the included alleles (only `*16`, `*20`, and `*22` have variant positions specified, document from March 2022)
PharmCAT relies on [PharmVar CYP3A4 allele definitions](https://www.pharmvar.org/gene/CYP3A4). The CYP3A4 `*16`, `*20` and `*22` definitions are the same
in the DPWG CYP3A4 gene document and PharmVar.


### DPYD

The CPIC DPYD allele definition file includes variants that have been assigned normal function (activity value 1),
decreased function (activity value 0.5), or no function (activity value 0). According to the
CPIC guideline for fluoropyrimidines and DPYD [(PMID:29152729)](https://pubmed.ncbi.nlm.nih.gov/29152729/),
the DPYD phenotype is assigned using a gene activity score, which is calculated as the sum of the activity scores of the
two DPYD variants with the lowest variant activity score. 

> If two different decreased/no function variants are present, they are presumed to be on different gene copies.
> Irrespective of the presence of decreased/no function variants, patients may carry multiple normal function variants.
> Common normal function variants may be located on the same gene copy as other normal function variants or
> decreased/no function variants.

More details are available on the [PharmGKB's CPIC DPYD reference page](https://www.pharmgkb.org/page/dpydRefMaterials).
The _DPYD Allele Functionality Table_ has function assignment information and the _DPYD Diplotype-Phenotype Table_
includes example translations considering one or two variants. 

### Calling DPYD named alleles

Note: the `combinations` research flag is ignored when calling DPYD.

{: .info}
> Effectively phased data is unphased data that is homozygous at all positions or is heterozygous at a single a
> position. Since we can effectively predict the alleles on each chromosome in this situation, we can treat this data as
> if it were phased data.


#### Phased or effectively phased data

If phased or effectively phased data is provided in the VCF file, the `Named Allele Matcher` produces an output that
lists all detected variants per strand. For example: `[c.498G>A + c.2582A>G]/[c.2846A>T + c.2933A>G]`.
If no variants are found on a strand, the `Named Allele Matcher` returns `Reference` for that allele.

#### Unphased data

If unphased data (that cannot be considered effectively phased) is provided in the VCF file, and the data are not
homozygous at all positions, the `Named Allele Matcher` will not attempt to call a diplotype. Instead, it produces a
list of all detected DPYD variants in the sample. It will, however, check if variants can be called on both strands.
If so, it will call the variant twice.  For example: `c.1627A>G (*5)`, `c.1905+1G>A (*2A)`, `c.1905+1G>A (*2A)`.
If the sample doesn't contain variants at the positions from the allele definition file, and/or if those positions are
omitted from the VCF file, the `Named Allele Matcher` returns `Reference`.

### Calling DPYD allele functionality and phenotype

The report lists the respective allele functionality for each variant and for `Reference`. If a diplotype was called
from phased/all-homozygous data, the lowest function variants on each strand will be used to determine the gene activity
score and DPYD phenotype. Otherwise, the two lowest function variants found are used to determine the gene activity
score and DPYD phenotype. The phenotype and gene activity score are utilized to retrieve the corresponding drug
recommendations.

### Calling the HapB3 allele

DPYD variants included in the CPIC and DPWG guidelines are single variants except for
`c.1129-5923C>G, c.1236G>A (HapB3)`. HapB3 consists of an exonic SNP at rs56038477 and an intronic SNP at rs75017182.
Both variants have been found in high linkage, however this is not the case for every sample. HapB3's decreased function
is thought to be due to the rs75017182 (intron) SNP, while the exon SNP is a synonymous substitution
(`p.Glu412=`). 

The `Named Allele Matcher` will rely on the intronic SNP (rs75017182) to call HapB3 if it's available, and only use the
exonic SNP (rs56038477) when it is not. If the VCF file input is conflicting between the two variants, or if the exonic
SNP is missing, the `Named Allele Matcher` will use the input for the intronic SNP to determine if HapB3 is present and
include a warning. If the intronic SNP is missing in the VCF file, the `Named Allele Matcher` will use the input for the
exonic SNP to determine if HapB3 is present and include a warning that the intronic SNP is missing and should be
genotyped to confirm the presence or absence of HapB3 and decreased function.


### DPWG recommendation

As of October 2022, DPWG recommendations are available for 4 DPYD variations:

* `c.1129-5923C>G, c.1236G>A (HapB3)`
* `2846A>T`
* `c.1905+1G>A (*2A)`
* `1679T>G (*13)`

When inferring gene activity score and phenotype from the two variants with the lowest activity value (unphased data)
or the lowest per strand (phased data), PharmCAT uses the variants that are included in both CPIC and DPWG if more than
one variant with the same activity value is found.

For example, if a sample has been called with a diplotype of `[c.1905+1G>A (*2A) + c.2933A>G]/c.498G>A`, which is
composed of `c.1905+1G>A (*2A)` (no function), `c.2933A>G` (no function, unknown to DPWG), `c.498G>A` (normal function),
the inferred diplotype used to look up the DPYD phenotype and recommendation will be `c.1905+1G>A (*2A)/c.498G>A` rather
than `c.2933A>G/c.498G>A`.

Furthermore, to increase the likelihood of a match with DPWG, PharmCAT treats any variant that is unknown to DPWG but
has a normal function in CPIC as `Reference` (i.e. normal function).

For example, in the above sample, the inferred genotype `c.1905+1G>A (*2A)/c.498G>A` will be further translated to
`c.1905+1G>A (*2A)/Reference` and used to query DPWG data. Since `c.1905+1G>A (*2A)` is a no function variant included
in the DPWG data, DPWG guidance for `c.1905+1G>A (*2A)/Reference` will be included in the report.


## G6PD

As of June 2023, PharmCAT does not call the following two _G6PD_ alleles:
- Mediterranean Haplotype (different from the Mediterranean, Dallas, Panama, Sassari, Cagliari, Birmingham)
- Villeurbanne

`Mediterranean Haplotype` is defined by rs5030868 (shared with `Mediterranean, Dallas, Panama, Sassari, Cagliari, Birmingham`
allele) and rs2230037 G>A (exclusive to `Mediterranean Haplotype`). PharmCAT excluded rs2230037 as the reference
genotype in the Mediterranean Haplotype (G) does not align with the GRCh38 reference base (A). As a result,
Mediterranean Haplotype is currently not called in PharmCAT.

Instead, individuals with the `Mediterranean Haplotype` will be recognized as having
`Mediterranean, Dallas, Panama, Sassari, Cagliari, Birmingham`. Both alleles have the same function (deficient);
therefore, the phenotype assignment and subsequent drug prescribing recommendations are not affected.

It remains challenging to properly represent SNPs and INDELs that share the same genetic positions in a VCF file.
This affects NC_000023.11:g.154532991_154532993delGGT and rs5030869 in the _G6PD_ allele definitions.
Left-aligned, normalized NC_000023.11:g.154532991_154532993delGGT is located at the same genomic position as rs5030869.
This poses nontrivial challenges to properly format variant representations for both SNP and INDEL in a VCF file and for
PharmCAT to accurately recognize any genetic variants of this kind at any allele-defining position.
Due to this challenge, NC_000023.11:g.154532991_154532993delGGT, which defines `Villeurbanne`,
is currently excluded from the _G6PD_ allele definitions in PharmCAT.


## NAT2

[NAT2 star alleles](https://www.pharmvar.org/gene/NAT2) are largely characterized by five key variants:
* rs1801279 (`191G>A`)
* rs1801280 (`341T>C`)
* rs1799930 (`590G>A`)
* rs1208 (`803G>A`)
* rs1799931 (`857G>A)`

These are found in different combinations and with other variants. Due to the number of star alleles with these five
variants, multiple diplotype assignments will often be possible for unphased data.
In such situations, PharmCAT typically defaults to selecting the diplotype containing the star allele with the most
variants (see [Scoring](/methods/NamedAlleleMatcher-101/#scoring) for details).

However, review of historical star allele frequencies and our recent analysis of the UK Biobank 200K phased data
(manuscript in preparation) confirmed that NAT2 star alleles most variants are not always the most common across
biogeographic populations.   
We have used this frequency data to change the default PharmCAT behavior for unphased NAT2 data.
Below is a list of possible diplotype results for the same data, with the default call shown in bold.

{% include_relative unphasedPriorities-NAT2.md %}

The list of potential calls assumes all variants of interest are available.
If there are missing variants, the `Named Allele Matcher` will still try to pick the priority diplotype as long as all
diplotypes in the list are present in the call.

For example, given the following variants (with all other positions being the reference allele):

| rs1801279 | rs1801280 | rs1799930 | rs1208 | rs1799931 |
|:----------|:----------|:----------|:-------|:----------|
| G/A       | T/C       | G/A       | G/A    | G/G       |

The `Named Allele Matcher` would normally call: `*6/*29` or `*30/*46` or `*5/*15`.

However, with the NAT2 prioritization scheme, `*5/*15` will be called instead.  
A note will be added to the result whenever this occurs.


## RYR1

In December 2023, [CPIC](https://www.pharmgkb.org/guidelineAnnotation/PA166303941) added additional variants to align
with the ClinGen malignant hyperthermia susceptibility (MHS) variant curation expert panel (VCEP), which developed and
published recommendations for RYR1 pathogenicity classifications in malignant hyperthermia susceptibility.
Pathogenic and likely pathogenic variants are assigned a CPIC function of Malignant Hyperthermia associated,
variants of uncertain significance (VUS) are mapped to uncertain function, and benign and likely benign variants are
assigned normal function.

RYR1 phenotypes are determined based on the function combinations of two RYR1 variants. In case >2 variants are found,
the variants with Malignant Hyperthermia associated function should be used first. If 2 Malignant Hyperthermia
associated variants are found, those are assumed on different chromosomes.


### Calling RYR1 named alleles

Note: the `combinations` research flag is ignored when calling RYR1.

{: .info}
> Effectively phased data is unphased data that is homozygous at all positions or is heterozygous at a single a
> position. Since we can effectively predict the alleles on each chromosome in this situation, we can treat this data as
> if it were phased data.

#### Phased or effectively phased data

If phased or effectively phased data is provided in the VCF file, the `Named Allele Matcher` produces an output that
lists all detected variants per strand. For example: `[c.97A>G + c.103T>C]/[c.418G>A + c.14918C>T]`.
If no variants are found on an allele, the `Named Allele Matcher` returns `Reference` for that allele.

#### Unphased data

If unphased data (that cannot be considered effectively phased) is provided in the VCF file, and the data are not
homozygous at all positions, the `Named Allele Matcher` will not attempt to call a diplotype. Instead, it produces a
list of all detected RYR1 variants in the sample. It will, however, check if variants can be called on both strands.
If so, it will call the variant twice.  For example: `c.97A>G`, `c.103T>C`, `c.103T>C`.
If the sample doesn't contain variants at the positions from the allele definition file, and/or if those positions are
omitted from the VCF file, the `Named Allele Matcher` returns `Reference`.

### Calling RYR1 allele functionality and phenotype

The report lists the respective allele functionality for each variant and for `Reference`. If a diplotype was called
from phased/all-homozygous data, PharmCAT determines if variants assigned the function
"_Malignant Hyperthermia associated_" are present, and those variants are used to assign the phenotype. Similarly, for
unphased data, variants assigned "_Malignant Hyperthermia associated_" function are used to determine the
"_Malignant Hyperthermia susceptible_" phenotype. The phenotype is used to retrieve the corresponding drug
recommendations.


## SLCO1B1

CPIC provides recommendations based on the SLCO1B1 star allele genotype. The CPIC guideline for statins and SLCO1B1,
ABCG2, and CYP2C9 [(PMID:35152405)](https://pubmed.ncbi.nlm.nih.gov/35152405/) includes the following excerpt:

> The most common and well-studied variant in SLCO1B1 is c.521T>C (rs4149056), and can be genotyped alone
> (e.g., PCR-based single SNV assay) or multiplexed on a variety of array-based platforms.
> All SLCO1B1 genetic tests should interrogate c.521T>C; however, while other less common variants in this gene may have
> limited evidence to guide action, they may also be important.

PharmCAT attempts to determine the star allele genotype for SLCO1B1, but in cases where no call can be determined, it
provides the CPIC recommendation based on the rs4149056 variant genotype.
