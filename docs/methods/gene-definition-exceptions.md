---
title: Gene Definition Exceptions
permalink: methods/gene-definition-exceptions/
---

# Gene Definition Exceptions

The genotype determination is based on [CPIC gene definition
tables](https://www.pharmgkb.org/page/pgxGeneRef), with
the following modifications:

1.  CYP3A5
    1.  Includes the a 'Y' at position rs28383479 to reflect the
        statement 'Existence of the CYP3A5\*3 polymorphism 6986A\>G on
        the same allele cannot be excluded' from the [Pharmacogene
        Variation Consortium](https://www.pharmvar.org/gene/CYP3A5).
        The inclusion of the ambiguous nucleotide 'Y' in the \*3
        definition results in a \*1/\*3 genotype call in case
        rs28383479 C\>T is present together with rs776746 T\>C since
        the match is based on the longer allele.
2. CYP2C9
    1.  rs1057911 (G475=) has been removed from the CYP2C9\*3 and 18
        definition. In the supporting article cited for the CYP2C9\*18
        allele definition (PMID: 15371982,
        https://www.pharmvar.org/gene/CYP2C9) this synonymous variant
        is described in linkage with rs1057910, the defining variant
        for \*3. However, the rs1057911 (G475=) variant seems to be
        present in further alleles beside \*3 and \*18.
3. CYP2C19
    1.  rs17885098 (99C>T) and rs3758581 (80161A>G, I331V) are not included in the PharmCAT allele definition table given the presence of these variants in almost all star alleles, including sub alleles of \*1.
    2.  rs4917623 (87106T\>C) is an intronic variant in \*18 and \*19. Based on Table 3 of the reference (PMID: 16141610) for both alleles in the [Pharmacogene Variation Consortium](https://www.pharmvar.org/gene/CYP2C19), the variant exists in other star alleles besides \*18 and \*19. Therefore, this intronic variant is not included in the PharmCAT allele definition table.
    3.  Sequence information in the CYP2C19 5 prime region was not consistently included in the star allele definitions (see archived CYP2C19 version on the Pharmacogene Variation Consortium website). Five prime variations are included in the PharmCAT allele definition table if a) currently linked to functional relevance (-806C>T, \*17) or b) the only identifying variant of the star allele (-1041G>A, \*27).
    4.  The variant 12802G>A (R150H) is the defining variant of the \*11 allele but can also exist on the \*2 background (see PharmVar CYP2C19\*2.010). To account for this possibility an 'R' at position 12802G>A (R150H) is included in the \*2 definition.
4. SLCO1B1
    1.  CPIC provides recommendations based on the rs4149056 genotype or
        the SLCO1B1 star allele genotype as reporting options in Table 1 of the
        SLCO1B1/simvastatin guideline \[PMID: 24918167\]. While
        PharmCAT attempts to determine the star allele genotype for
        SLCO1B1, in case no call can be determined it provides the
        CPIC recommendation based on the rs4149056 variant genotype.
5.  UGT1A1
    1.  Phased data: The CPIC recommendations are provided using the
        subject's phased UGT1A1 genotype.
    2.  Unphased data: All variant alleles found are listed, based on
        the vcf input and the UGT1A1 definition table. The CPIC
        recommendations are provided based on an [exception logic](calling/UGT1A1).
    3.  The allele definition table does not include a row for only \*80
        since the presence of the rs887829 variant without the
        detection of \*28 or \*37 is assigned uncertain function and
        there are not enough clinical data to predict metabolizer
        status with certainty. In the UGT1A1 gene section rs887829 is
        noted.
6. CFTR
    1.  In the CPIC CFTR allele definition table the F508del variation
        is represented as F508del(CTT) rs113993960 and F508del(TCT)
        rs199826652. In the PharmCAT CFTR definition file only
        F508del(CTT) rs113993960 is used.
7. TPMT
    1.  \*1S variant can be observed in combination with other TPMT
        variants on the same strand which confounds star allele
        calling. To reduce TPMT 'no calls' in PharmCAT the TPMT genotype
        is determined without considering the \*1S allele. However,
        the variant (rs2842934) is provided in the TPMT gene section.
8. DPYD
    1.  HGVS variant names are used to connect recommendation to DPYD variation in the 2017 CPIC guideline update for DPYD and fluoropyrimidines.
    2.  The PharmCAT DPYD gene definition table includes decreased and no function variants with strong and moderate evidence, see [Allele Functionality Table](https://www.pharmgkb.org/page/dpydRefMaterials). If none of these variants are detected in the VCF, recommendations are provided based on two functional alleles.
    3.  The 2017 CPIC DPYD guideline update includes normal function variants with strong, moderate or weak evidence. These variants are not included in the PharmCAT allele definition table. However, these variants are listed in the DPYD gene section of the PharmCAT report.
    4.  The CPIC DPYD Allele Functionality Table includes further no and
    decreased function variants with in-vitro data only evidence.
    According to the guideline, to date, there are no studies linking
    the variants listed in the "in vitro data only and/or limited
    clinical/ex vivo data" category as decreased or no function
    variants directly to toxicity related to fluoropyrimidines and
    therefore, these variants are not specifically included in the
    CPIC recommendations. These variants are not included in the
    PharmCAT DPYD gene definition table, however, a list of these
    variants is provided in the DPYD gene section under "Other
    Positions of Interest".
