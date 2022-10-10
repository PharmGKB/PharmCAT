---
parent: Methods
title: Gene Definition Exceptions
permalink: methods/Gene-Definition-Exceptions/
nav_order: 6
---
# Gene Definition Exceptions

The genotype determination is based on [CPIC gene definition
tables](https://www.pharmgkb.org/page/pgxGeneRef), with
the following modifications:
    
## SLCO1B1

CPIC provides recommendations based on the SLCO1B1 star allele genotype. The CPIC guideline for statins and SLCO1B1, ABCG2, and CYP2C9 [(PMID:35152405)](https://pubmed.ncbi.nlm.nih.gov/35152405/) includes the following excerpt:

> The most common and well-studied variant in SLCO1B1 is c.521T>C (rs4149056), and can be genotyped alone (e.g., PCR-based single SNV assay) or multiplexed on a variety of array-based platforms. All SLCO1B1 genetic tests should interrogate c.521T>C; however, while other less common variants in this gene may have limited evidence to guide action, they may also be important. 

PharmCAT attempts to determine the star allele genotype for SLCO1B1, in case no call can be determined it provides the CPIC recommendation based on the rs4149056 variant genotype.


## DPYD

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
The_DPYD Allele Functionality Table_ has function assignment information and the _DPYD Diplotype-Phenotype Table_
includes example translations considering one or two variants. 

### Calling DPYD named alleles

Note: the combination research flag is ignored when calling DPYD.

#### Phased data

If phased data is provided in the VCF file, or if the data are homozygous at all positions, the `Named Allele Matcher`
produces an output that lists all detected variants per allele. For example:
`[c.498G>A + c.2582A>G]/[c.2846A>T + c.2933A>G]` . If no variants are found on an allele, the `Named Allele Matcher`
returns `Reference` for that allele.

#### Unphased data

If unphased data is provided in the VCF file, and the data are not homozygous at all positions, the
`Named Allele Matcher` will not attempt to call a diplotype. Instead, it produces a list of all detected DPYD variants
in the sample. It will, however, check if variants can be called on both strands. If so, it will call the variant twice.
For example: `c.1627A>G (*5)`, `c.1905+1G>A (*2A)`, `c.1905+1G>A (*2A)`. If the sample doesn’t contain variants at the
positions from the allele definition file, and/or if those positions are omitted from the vcf file, the 
`Named Allele Matcher` returns `Reference`.

### Calling DPYD allele functionality and phenotype

The report lists the respective allele functionality for each variant and for `Reference`. If a diplotype was called
from phased/all-homozygous data, the lowest function variants on each strand will be used to determine the gene activity
score and DPYD phenotype. Otherwise, the two lowest function variants found are used to determine the gene activity
score and DPYD phenotype. The phenotype and gene activity score are utilized to retrieve the corresponding drug
recommendations.

### DPWG recommendation

As of October 2022, DPWG recommendations are available for 4 DPYD variation:

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
