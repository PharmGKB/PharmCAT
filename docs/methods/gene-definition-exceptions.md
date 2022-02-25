---
title: Gene Definition Exceptions
permalink: methods/gene-definition-exceptions/
---

# Gene Definition Exceptions

The genotype determination is based on [CPIC gene definition
tables](https://www.pharmgkb.org/page/pgxGeneRef), with
the following modifications:
    
## SLCO1B1

CPIC provides recommendations based on the SLCO1B1 star allele genotype. The CPIC guideline for statins and SLCO1B1, ABCG2, and CYP2C9 [(PMID:35152405)](https://pubmed.ncbi.nlm.nih.gov/35152405/) includes the following excerpt: "The most common and well-studied variant in SLCO1B1 is c.521T>C (rs4149056), and can be genotyped alone (e.g., PCR-based single SNV assay) or multiplexed on a variety of array-based platforms. All SLCO1B1 genetic tests should interrogate c.521T>C; however, while other less common variants in this gene may have limited evidence to guide action, they may also be important”. 

PharmCAT attempts to determine the star allele genotype for SLCO1B1, in case no call can be determined it provides the CPIC recommendation based on the rs4149056 variant genotype.
