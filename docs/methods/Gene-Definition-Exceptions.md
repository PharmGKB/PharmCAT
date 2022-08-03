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

CPIC provides recommendations based on the SLCO1B1 star allele genotype. The CPIC guideline for statins and SLCO1B1, ABCG2, and CYP2C9 [(PMID:35152405)](https://pubmed.ncbi.nlm.nih.gov/35152405/) includes the following excerpt: "The most common and well-studied variant in SLCO1B1 is c.521T>C (rs4149056), and can be genotyped alone (e.g., PCR-based single SNV assay) or multiplexed on a variety of array-based platforms. All SLCO1B1 genetic tests should interrogate c.521T>C; however, while other less common variants in this gene may have limited evidence to guide action, they may also be important”. 

PharmCAT attempts to determine the star allele genotype for SLCO1B1, in case no call can be determined it provides the CPIC recommendation based on the rs4149056 variant genotype.

## DPYD

The CPIC DPYD allele definition file includes variants that have been assigned normal (activity value 1), decreased (activity value 0.5), and no function (activity value 0) (see DPYD functionality file). The CPIC Guideline for Fluoropyrimidines and DPYD [(PMID:29152729)](https://pubmed.ncbi.nlm.nih.gov/29152729/) states:

“The DPD phenotype is assigned using a gene activity score, which is calculated as the sum of the activity scores of the two DPYD variants with the lowest variant activity score.”

“If two different decreased/no function variants are present, they are presumed to be on different gene copies.”

“Irrespective of the presence of decreased/no function variants, patients may carry multiple normal function variants. Common normal function variants may be located on the same gene copy as other normal function variants or decreased/no function variants.”

The CPIC DPYD Genotype to Phenotype file includes example translations considering one or two variants. 

All variants in the CPIC DPYD allele definition file are considered for the genotype assignment by the PharmCAT Named Allele Matcher. This potentially results in the detection of more than two DPYD variants. To be able to determine a DPYD phenotype the variants are reported and translated as followed:

__Phased data:__

If phased data is provided in the vcf file, the Named Allele Matcher produces an output that lists all detected variant per allele. 

 - If the sample only includes none or one variant per allele the genotype is given as e.g., c.1627A>G (\*5)/c.1905+1G>A (\*2A) with a corresponding functionality (One no function allele and one normal function allele) and the phenotype assignment (Intermediate Metabolizer) in the json output and HTML/PDF report summary table. The genotype and phenotype are also included in the drug sections of the report.
  
 - In case the sample includes more than one DPYD variant per allele e.g., c.498G>A + c.2582A>G/c.2846A>T + c.2933A>G, the genotype is listed as such in the json output and HTML/PDF report summary table. The report drug section will include the genotype with the two lowest activity values (variant activity scores), which is used to determine the phenotype and gene activity score and the corresponding drug recommendations. This information is also available in the json format.

__Unphased data:__

If unphased data is provided in the vcf file, the Named Allele Matcher produces a list of all detected DPYD variants in the sample. The Phenotyper assigns function to each variant which is provided in the json output and HTML/PDF report summary table. The report drug section will include the genotype with the two lowest activity values (variant activity scores), which is used to determine the phenotype and gene activity score and the corresponding drug recommendations. For normal function variants, reference is used in the calculated genotype to present normal function alleles. This information is also available in the json format.
