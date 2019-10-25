---
title: Gene Definition Exceptions
permalink: methods/gene-definition-exceptions/
---

# Gene Definition Exceptions

The genotype determination is based on [CPIC gene definition
tables](https://www.pharmgkb.org/page/pgxGeneRef), with
the following modifications:

1.  CYP3A5
    1.  The variant 12952T>C (rs55965422) is the defining variant of the \*5 allele but can also exist on the \*3 background (see [Pharmacogene Variation Consortium](https://www.pharmvar.org/gene/CYP3A5) CYP3A5\*3G). To account for this possibility an 'R' at position 12952T>C (rs55965422) is included in the \*3 definition. The inclusion of the ambiguous nucleotide 'R' in the \*3 definition results in a \*1/\*3 genotype call in case rs55965422 A\>G is present together with rs776746 T\>C since the match is based on the longer allele.
2. CYP2C19
    1.  rs3758581 (80161A>G, I331V) is not included in the PharmCAT allele definition table given the presence of these variant in almost all star alleles, including sub alleles of \*1.
    2.  \*36 and \*37 are not included in the table since these alleles are a full or a partical deletion of CYP2C19. 
4. SLCO1B1
    1.  CPIC provides recommendations based on the rs4149056 genotype or the SLCO1B1 star allele genotype as reporting       options in Table 1 of the SLCO1B1/simvastatin guideline \[PMID: 24918167\]. While PharmCAT attempts to determine the star allele genotype for SLCO1B1, in case no call can be determined it provides the CPIC recommendation based on the rs4149056 variant genotype.
5.  UGT1A1
    1.  Phased data: The CPIC recommendations are provided using the subject's phased UGT1A1 genotype.
    2.  Unphased data: All variant alleles found are listed, based on the vcf input and the UGT1A1 definition table. The CPIC recommendations are provided based on an [exception logic](calling/UGT1A1).
    3.  The allele definition table does not include a row for only \*80 since the presence of the rs887829 variant without the detection of \*28 or \*37 is assigned uncertain function and there are not enough clinical data to predict metabolizer status with certainty. In the UGT1A1 gene section rs887829 is noted.
6. CFTR
    1.  In the CPIC CFTR allele definition table the F508del variation is represented as F508del(CTT) rs113993960 and F508del(TCT) rs199826652. 
7. TPMT
    1.  \*1S variant can be observed in combination with other TPMT variants on the same strand which confounds star allele calling. To reduce TPMT 'no calls' in PharmCAT the TPMT genotype is determined without considering the \*1S allele. However, the variant (rs2842934) is provided in the TPMT gene section.
8. DPYD
    1.  HGVS variant names are used to connect recommendation to DPYD variation in the 2017 CPIC guideline update for DPYD and fluoropyrimidines.
    2.  The PharmCAT DPYD gene definition table includes decreased and no function variants with strong and moderate evidence, see [Allele Functionality Table](https://www.pharmgkb.org/page/dpydRefMaterials). If none of these variants are detected in the VCF, recommendations are provided based on two functional alleles.
    3.  c.1129-5923C>G is likely the HapB3 haplotype causal variant and is currently listed as single variant repesenting this haplotype. Proxy SNPs are c.1236G>A (E412E) is currently not included in the PharmCAT DPYD allele definition table.
    4.  The 2017 CPIC DPYD guideline update includes normal function variants with strong, moderate or weak evidence. These variants are not included in the PharmCAT allele definition table. However, these variants are listed in the DPYD gene section of the PharmCAT report.
    5.  The CPIC DPYD Allele Functionality Table includes further no and decreased function variants with in-vitro data only evidence. According to the guideline, to date, there are no studies linking the variants listed in the "in vitro data only and/or limited clinical/ex vivo data" category as decreased or no function variants directly to toxicity related to fluoropyrimidines and therefore, these variants are not specifically included in the CPIC recommendations (DPYD Diplotype to Phenotype table). These variants are not included in the PharmCAT DPYD gene definition table, however, a list of these variants is provided in the DPYD gene section under "Other Positions of Interest".
