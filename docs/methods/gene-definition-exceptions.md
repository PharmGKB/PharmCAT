---
title: Gene Definition Exceptions
permalink: methods/gene-definition-exceptions/
---

# Gene Definition Exceptions

The genotype determination is based on [CPIC gene definition
tables](https://www.pharmgkb.org/page/pgxGeneRef), with
the following modifications:

1. CYP2C19
    1.  rs3758581 (80161A>G, I331V) is not included in the allele definition used in PharmCAT given the presence of this variant in almost all star alleles, including  \*1. Therefore, there is no distinction between \*1 and \*38 (which is the reference allele without any variants) and PharmCAT reports \*1.
    
2. SLCO1B1
    1.  CPIC provides recommendations based on the rs4149056 genotype or the SLCO1B1 star allele genotype as reporting options in Table 1 of the SLCO1B1/simvastatin guideline \[PMID: 24918167\]. While PharmCAT attempts to determine the star allele genotype for SLCO1B1, in case no call can be determined it provides the CPIC recommendation based on the rs4149056 variant genotype.
