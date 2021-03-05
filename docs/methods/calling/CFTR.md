---
title: Calling CFTR
permalink: methods/calling/CFTR/
---

# Calling and Reporting on CFTR in PharmCAT

### If caller returns reference/reference

Display \'No CPIC variants found\' in the diplotype or genotype field
and N/A for allele functionality and phenotype.

In the recommendation section:  
\'No CPIC variants are found. CPIC does not provide recommendations for
this genotype. CPIC guidelines reflect the alleles/genotypes known and
considered by the guideline authors for inclusion by the time of
publication. CPIC guidelines are periodically updated. Read the [full
guideline](https://www.pharmgkb.org/guideline/PA166114461).\'

### If caller returns variant/variant with variant being any name we have in definition file 

Display caller output in diplotype or genotype field and N/A for allele
functionality and phenotype.

In recommendation section:  
For combinations including F508del(CTT), F508del(TCT), G551D, G1244E,
G1349D, G178R, G551S, S1251N, S1255P, S549N, S549R(A\>C), S549R(T\>G),
R117H use recommendation as in guideline tool with \'other\' being Reference.

### If caller returns variant/reference with variant being any name we have in definition file 

Display name from caller output (heterozygous), e.g. R553X
(heterozygous), in diplotype or genotype field and N/A for allele
functionality and phenotype.

Diplotype or genotype = \*variant\* (heterozygous)  
Allele Functionality = N/A  
Phenotype = N/A

In recommendation section:

For a heterozygous variant of any of F508del(CTT), F508del(TCT), G551D,
G1244E, G1349D, G178R, G551S, S1251N, S1255P, S549N, S549R(A\>C),
S549R(T\>G), R117H use recommendation as in guideline tool in
combination with \'other\'.
