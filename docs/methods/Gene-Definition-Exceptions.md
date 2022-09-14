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

The CPIC DPYD allele definition file includes variants that have been assigned normal (activity value 1), decreased (activity value 0.5), and no function (activity value 0) (see DPYD functionality file). The CPIC Guideline for Fluoropyrimidines and DPYD [(PMID:29152729)](https://pubmed.ncbi.nlm.nih.gov/29152729/) states:

> The DPD phenotype is assigned using a gene activity score, which is calculated as the sum of the activity scores of the two DPYD variants with the lowest variant activity score.
> 
> If two different decreased/no function variants are present, they are presumed to be on different gene copies.
>
> Irrespective of the presence of decreased/no function variants, patients may carry multiple normal function variants. Common normal function variants may be located on the same gene copy as other normal function variants or decreased/no function variants.

The CPIC DPYD Genotype to Phenotype file includes example translations considering one or two variants. 

All variants in the CPIC DPYD allele definition file are considered for the genotype assignment by the PharmCAT Named Allele Matcher. This potentially results in the detection of more than two DPYD variants. To be able to determine a DPYD phenotype the variants are reported and translated as followed:

### Phased data

If phased data is provided in the vcf file, the Named Allele Matcher produces an output that lists all detected variant per allele, e.g., c.498G>A + c.2582A>G/c.2846A>T + c.2933A>G (sample includes more than one DPYD variant per allele) or c.1627A>G (\*5)/c.1905+1G>A (\*2A) (sample one variant per allele).

The report lists the respective allele functionality for each variant or the reference. If the sample only includes none or one variant per allele, the genotype can be directly used to detemined the gene activity score and DPYD phenotype. In case the sample includes more than one DPYD variant per allele, the genotype to determine the gene activity score and DPYD phenotype is inferred based on the variant with the lowest activity value (variant activity scores) per allele. The phenotype and gene activity score are utelized to retrieve the corresponding drug recommendations. 

### Unphased data

If unphased data is provided in the vcf file, the Named Allele Matcher produces a list of all detected DPYD variants in the sample. 

The report includes the variants and the respective allele functionality. If the sample only includes none, one or two variants (variants assumed on different gene copies - see above), the genotype can be directly used to detemined the gene activity score and DPYD phenotype. In case the sample includes more than two DPYD variants, the genotype to determine the gene activity score and DPYD phenotype is inferred based on the two variants with the lowest activity value. The phenotype and gene activity score are utelized to retrieve the corresponding drug recommendations. 

### DPWG recommendation 


### Mark's original

Note: the combination research flag is ignored when calling DPYD.

When calling alleles, the `Named Allele Matcher` will only attempt to call a diplotype if the data is effectively phased (i.e. actually phased or unphased but homozygous at all positions).  If there is a simple match, it will return that.  If not, it will attempt to find combinations.  Unlike normal combination scoring, it prioritizes combinations and will return the diplotype with the most combinations.

If no diplotype call is possible, or the data is effectively unphased, the `Named Allele Matcher` will only call possible haplotypes based on key positions only (i.e. will not assume reference on positions with no alleles), and will _not_ attempt to call a diplotype.
