---
title: Disclaimers & Other Information
permalink: Disclaimers/
nav_order: 7
---
# Disclaimers and Other Information

__Liability: PharmCAT assumes no responsibility for any injury to person or damage to persons or property arising out of,
or related to any use of PharmCAT, or for any errors or omissions. The user recognizes that PharmCAT is a research tool
and that they are using PharmCAT at their own risk.__

## A. Allele and Genotype Determination

1. PharmCAT uses gene allele definitions included in the CPIC database, with exceptions as noted
   in [Gene Definition Exceptions document](https://pharmcat.org/methods/Gene-Definition-Exceptions/). For allele
   definitions and the positions used in PharmCAT, see
   the [gene definition tables](https://github.com/PharmGKB/PharmCAT/releases).

2. Genes with DPWG recommendations that are not included in CPIC are discussed in Section C.   

3. PharmCAT results are dependent on the supplied VCF calls for the queried positions (for technical information about
   PharmCAT input formatting and requirements, please go to [pharmcat.org](https://pharmcat.org/)). PharmCAT does not
   assume any reference calls for positions missing from the submitted VCF file; all missing queried positions are not
   considered in the allele determination process. See the
   [*pharmcat_positions* file](https://github.com/PharmGKB/PharmCAT/releases) for which positions are
   queried in the VCF file. Missing positions might alter the assigned genotype and subsequent phenotype prediction and
   CPIC recommendation. If the supplied VCF file has missing positions, those positions will be noted in
   Section III: Allele Matching Details for each gene of this report. For the most reliable allele determination,
   reference calls as well as variant calls in the VCF file for every queried position must be provided by the user.
   If an allele that includes a missing position is defined by an additional position(s) for which calls are provided,
   the allele will be considered in the matching process based on the available information. This might lead to the output
   of multiple possible genotypes that received the same highest matching score. In addition, alternate calls with a lower
   score is also possible. For instructions on getting PharmCAT to output all possible matching genotypes, consult
   [the documentation](https://pharmcat.org/using/Running-PharmCAT/#named-allele-matcher).

4. For cytochrome P450 genes, TPMT, NUDT15, UGT1A1, and SLCO1B1, the \*1 allele is defined by the absence of variation
   specified in the gene definition tables. This allele cannot be identified by variants; rather, \*1 is assigned by
   default when no variation for the queried positions is reported in the submitted VCF file. The same is true for all
   other genes with multiple variant positions in the definition table (CACNA1S, CFTR, DPYD, RYR1): the reference
   sequence is the default result when variants are not reported in the VCF file. It is always possible un-interrogated
   variation can occur which could potentially affect allele function, but because it is undetected, the assignment
   would be defaulted to a \*1 (or reference) allele and normal function.

5. For all genes, variation reported in the VCF file but NOT included in the gene definition table will not be
   considered during allele assignment. There is a possibility that any such variation results in a reduced or no
   activity allele which could lead to inaccurate phenotype and CPIC recommendation, similar to the situation in point
   3, above.

6. Nucleotide base calls are displayed on the positive chromosomal strand regardless of the gene strand; further
   information is provided under Gene-specific warnings in Section III: Allele Matching Details.
   
7. Structural variation star alleles that cannot be detected using VCF file data:
   - CYP2B7-CYP2B6 hybrids: CYP2B6\*29, CYP2B6\*30
   - Partial and whole gene deletions: CYP2C19\*36, CYP2C19\*37, CYP4F2\*16, SLCO1B1\*48, SLCO1B1\*49  

8. PharmCAT matches variants to genotypes assuming unphased data (unless phased data is provided in the VCF file and
   noted as such, see [pharmcat.org](https://pharmcat.org/) for details). The assumption is that defined alleles exist
   in trans configuration, i.e. on opposite chromosomes, with exceptions noted in Section III: Allele Matching Details
   under "Gene-specific warnings." However, in cases where an allele is defined by a combination of two or more
   variants, where each variant alone also defines an allele, the match is based on the longer allele. For example,
   TPMT\*3B is defined by one SNP, \*3C is defined by another SNP, and \*3A is defined by the combination of those two
   SNPs. In the case of unphased data that is heterozygous for both SNPs, the \*1/\*3A genotype is returned though the
   possibility of \*3B/\*3C cannot be ruled out.
   
   Below cases are summarized for which two calls with different scores are possible when provided unphased data and
   heterozygous calls for the variants that define the two alleles. The genotype with the higher score (longer allele)
   will be used to determine allele functionality, phenotype, and recommendation but the genotype with the lower score
   cannot be ruled out.

Table 1: Cases for which there is an overlap in the allele definitions.

| Gene    | Genotype (Higher Score)    | Phenotype          | Genotype (Lower Score)               | Phenotype                   |
|---------|----------------------------|--------------------|--------------------------------------|-----------------------------|
| UGT1A1  | \*1/\*80+\*28              | Intermediate       | \*28/\*80                            | Indeterminate               |
| UGT1A1  | \*1/\*80+\*37              | Intermediate       | \*37/\*80                            | Indeterminate               |
| TPMT    | \*1/\*3A                   | Intermediate       | \*3B/\*3C                            | Poor                        |
| NUDT15  | \*1/\*2                    | Intermediate       | \*3/\*6                              | Possible Intermediate       |
| CYP2C9  | \*1/\*71                   | N/A                | \*10/\*22                            | Indeterminate               |
| CYP2B6  | \*1/\*36                   | Intermediate       | \*6/\*22                             | Intermediate                |
| CYP2B6  | \*1/\*34                   | Intermediate       | \*33/\*36                            | Indeterminate               |
| CYP2B6  | \*1/\*6                    | Intermediate       | \*4/\*9                              | Intermediate                |
| CYP2B6  | \*1/\*7                    | Intermediate       | \*5/\*6                              | Intermediate                |
| CYP2B6  | \*1/\*13                   | Intermediate       | \*6/\*8                              | Intermediate                |
| SLCO1B1 | \*1/\*46                   | Decreased function | \*15/\*45                            | Possible Decreased Function |
| SLCO1B1 | \*1/\*20                   | Normal Function    | \*19/\*37                            | Indeterminate               |
| SLCO1B1 | \*1/\*12                   | Indeterminate      | \*2/\*10                             | Indeterminate               |
| SLCO1B1 | \*1/\*13                   | Indeterminate      | \*3/\*11                             | Indeterminate               |
| SLCO1B1 | \*1/\*14                   | Normal Function    | \*4/\*37                             | Indeterminate               |
| SLCO1B1 | \*1/\*15                   | Decreased function | \*5/\*37                             | Decreased function          |
| SLCO1B1 | \*1/\*25                   | Indeterminate      | \*4/\*28                             | Indeterminate               |
| SLCO1B1 | \*1/\*31                   | Decreased function | \*9/\*37                             | Decreased Function          |
| SLCO1B1 | \*1/\*32                   | Indeterminate      | \*4/\*24                             | Indeterminate               |
| SLCO1B1 | \*1/\*40                   | Indeterminate      | \*5/\*19                             | Possible Decreased Function |
| SLCO1B1 | \*1/\*43                   | Indeterminate      | \*4/\*44                             | Indeterminate               |
| CYP4F2  | \*1/\*4                    | N/A                | \*2/\*3                              | N/A                         |
| CYP3A4  | \*1/\*37                   | N/A                | \*3/\*22                             | N/A                         |
| CYP3A4  | \*1/\*38                   | N/A                | \*3/\*11                             | N/A                         |
| G6PD    | A- 202A_376G/B (reference) | Variable           | A/Asahi                              | Variable                    |
| G6PD    | B (reference)/Mt Sinai     | Variable           | A/Guadalajara                        | Variable                    |
| G6PD    | B (reference)/Santa Maria  | Variable           | A/Malaga                             | Variable                    |
| G6PD    | Ananindeua/B (reference)   | Variable           | A/Viangchan, Jammu                   | Variable                    |
| G6PD    | B (reference)/Hechi        | Variable           | Asahi/Viangchan, Jammu               | Deficient                   |
| G6PD    | B (reference)/Hermoupolis  | Variable           | Cassano/Union,Maewo, Chinese-2, Kalo | Deficient                   |


Table 2: Cases for which there is an overlap in the allele definitions because the definition of the non-\*1 allele in
the genotype with the higher score allows for reference or variant at the position that defines the first allele listed
in the genotype with the lower score. Both genotypes cannot be ruled out with unphased data if the position that
overlaps between the respective alleles is heterozygous (0/1) in addition to heterozygous calls for the other variants
that define the non-\*1 allele in the genotype with the higher score.

| Gene    | Genotype (Higher Score) | Phenotype    | Genotype (Lower Score)| Phenotype     |
| ------- | ----------------------- |--------------| --------------------- |---------------|
| CYP2C19 | \*1/\*4                 | Intermediate | \*17/\*4              | Intermediate  |
| CYP2C19 | \*1/\*2                 | Intermediate | \*11/\*2              | Intermediate  |
| CYP2C19 | \*1/\*35                | Intermediate | \*15/\*35             | Intermediate  |
| CYP2B6  | \*1/\*18                | Intermediate | \*4/\*18              | Indeterminate |


## B. CPIC Allele Function, Phenotype and Recommendation

All content is sourced from the [CPIC database](https://github.com/cpicpgx/cpic-data).


## C. DPWG Allele Function, Phenotype and Recommendation

1. PharmGKB annotates PGx-based drug dosing guidelines published by the [Royal Dutch Association for the Advancement of
   Pharmacy - Pharmacogenetics Working Group (DPWG)](https://www.pharmgkb.org/page/dpwg). PharmGKB curates allele 
   function assignments and phenotype mappings from the DPWG to provide genotype specific DPWG guideline
   recommendations. Where possible, PharmGKB maps DPWG terms to CPIC terms, as outlined on
   [PharmGKB](https://www.pharmgkb.org/page/dpwgMapping).

2. CYP3A4 is currently not part of a CPIC guideline. Since the DPWG CYP3A4 documentation includes limit variant
   notations for the included alleles (only \*16, \*20, and \*22 are specified) PharmCAT relies on
   [PharmVar CYP3A4 allele definitions](https://www.pharmvar.org/gene/CYP3A4). However, the CYP3A4\*16, \*20 and \*22
   definitions are the same in both sources.

3. The CPIC UGT1A1 allele definition file includes \*6, \*27, \*28, \*36, \*37, and \*80. Since the DPWG UGT1A1 document
   does not include allele definitions besides for the UGT1A1 TA box promoter polymorphism, PharmCAT only includes the
   UGT1A1 positions from the CPIC UGT1A1 allele definition file. Other UGT1A1 alleles can be supplied as
   [outside calls](https://pharmcat.org/using/Outside-Call-Format/) but not be determined from the VCF file by the
   Named Allele Matcher.

4. IMPORTANT: As of March 2022, gene information documents from the DPWG are no longer publicly available from the KNMP
   website. PharmGKB is currently providing PDF versions of these documents to users. These files were downloaded in
   February 2022. As such, it cannot be guaranteed that these documents match the mappings DPWG may use internally as
   there have been no publicly accessible updates since February 2022.


## D. PharmCAT Exceptions to the CPIC Guideline Gene List

1. PharmCAT does not determine CYP2D6, MT-RNR1, HLA-A, or HLA-B genotypes from the VCF file, but genotypes for CYP2D6,
   MT-RNR1, HLA-A, or HLA-B can be provided to PharmCAT from an outside source and the corresponding phenotype 
   prescribing recommendations will be included in the generated report. For the required format of the outside calls
   refer to the PharmCAT [outside calls documentation](https://pharmcat.org/using/Outside-Call-Format/).

2. CPIC has assigned function to the following CYP2D6 CNV alleles: \*1x2, \*1x≥3, \*2x2, \*2x≥3, \*3x2, \*4x2, \*4x≥3,
   \*6x2, \*9x2, \*10x2, \*17x2, \*29x2, \*35x2, \*36x2, \*41x2, \*41x3, \*43x2, \*45x2, *146x2. These alleles are part
   of the CPIC diplotype to phenotype translation and can be connected to recommendations. Other CNV notations from
   outside calls need to be mapped accordingly. 


## E. CPIC Guideline Disclaimers and Caveats

1. A version of the following quoted disclaimer is part of each CPIC guideline and applies to the CPIC recommendations
   as used in PharmCAT. For the full description of potential benefits and risks, additional considerations (general and
   specific to gene-drug pairs), limitations, information about respective gene nomenclature systems, potential
   drug-drug interactions and clinical factors to consider, please see individual CPIC guidelines at
   ([cpicpgx.org](https://cpicpgx.org)).
     
   1. "CPIC guidelines reflect expert consensus based on clinical evidence and peer-reviewed literature available at
       the time they are written and are intended only to assist clinicians in decision-making and to identify questions
       for further research. New evidence may have emerged since the time a guideline was submitted for publication.
       Guidelines are limited in scope and are not applicable to interventions or diseases not specifically identified.
       Guidelines do not account for all individual variations among patients and cannot be considered inclusive of all
       proper methods of care or exclusive of other treatments. It remains the responsibility of the health-care
       provider to determine the best course of treatment for a patient. Adherence to any guidelines is voluntary, with
       the ultimate determination regarding its application to be made solely by the clinician and the patient. CPIC
       assumes no responsibility for any injury to persons or damage to persons or property arising out of or related to
       any use of CPIC's guidelines, or for any errors or omissions." (PMID:
       [27997040](https://www.ncbi.nlm.nih.gov/pubmed/27997040))

   2. "Caveats: appropriate use and/or potential misuse of genetic tests. The application of genotype-based dosing is
       most appropriate when initiating therapy with a tricyclic. Obtaining a pharmacogenetic test after months of drug
       therapy may be less helpful in some instances, as the drug dose may have already been adjusted based on plasma
       concentrations, response, or side effects. Similar to all diagnostic tests, genetic tests are one of several
       pieces of clinical information that should be considered before initiating drug therapy." (PMID:
       [27997040](https://www.ncbi.nlm.nih.gov/pubmed/27997040))

2. CPIC guidelines reflect the alleles/genotypes known and considered by the guideline authors for inclusion by the time
   of publication, however they may be updated online at
   [cpicpgx.org](https://cpicpgx.org) and in the CPIC database in between publications. Additional alleles and/or more
   extensive allele definitions might exist by the representative gene nomenclatures for various genes.

3. CPIC is a registered service mark of the U.S. Department of Health & Human Services (HHS).


## F. PharmGKB Disclaimers and Caveats

PharmGKB is a registered service mark of the U.S. Department of Health & Human Services (HHS).
