---
title: Changelog
permalink: changelog/
nav_order: 6
---

# Change Log

This page highlights _major updates_ to PharmCAT (starting with v2.10.0).

If you are looking for details on releases, please visit the
[PharmCAT GitHub Releases page](https://github.com/PharmGKB/PharmCAT/releases).  If you'd like to learn more about how
we version PharmCAT, check out our page on [versioning and releases](/versioning).  Finally, you can subscribe to
PharmCAT release notifications by following
[these instructions](/versioning#subscribing-to-release-notifications). This method will require a GitHub account.


## v2.13.0

Added support for FDA drug guidance.


## v2.11.0

DPWG removed the recommendation for _F5_ and hormonal contraceptives for systemic use.
The guideline annotation was retired on PharmGKB, and resulted in the removal of _F5_ from PharmCAT.


## v2.10.0

#### _DPYD_

_DPYD_ HapB3 is defined by two variants, `c.1129-5923C>G` (an intronic variant) and `c.1236G>A` (an exonic variant). PharmVar and [CPIC](https://cpicpgx.org/guidelines/guideline-for-fluoropyrimidines-and-dpyd/) now include and list the variant `c.1129-5923C>G` separately as this is likely the causal variant leading to a decreased function. CPIC still retains the _DPYD_ HapB3 haplotype definition for cases where only the exonic variant `c.1236G>A` is tested (e.g. whole-exome sequencing). For more information, see CPIC Guideline for fluoropyrimidines and _DPYD_.

This change prompted the following PharmCAT updates:

- If `c.1236G>A` (the exonic defining variant of _DPYD_ HapB3) is missing, PharmCAT will use the intronic variant `c.1129-5923C>G`.
- If `c.1129-5923C>G` and `c.1236G>A` do not agree, PharmCAT will use `c.1129-5923C>G` and also report the presence of `c.1236G>A` in the input.
- _DPYD_ HapB3 will be reported if both of its defining variants `c.1236G>A` and `c.1129-5923C>G` are present and “in sync” with each other in the input VCF file.

#### _RYR1_

In December 2023, [CPIC](https://cpicpgx.org/guidelines/cpic-guideline-for-ryr1-and-cacna1s/) added additional 291 variants that were included in the ClinGen variant curation expert panel (VCEP) recommendations for RYR1 pathogenicity. PharmCAT accommodates these changes by updating the Named Allele Matcher module. The Name Allele Matcher now will report all _RYR1_ variants found in an individual. 

CPIC _RYR1_ phenotypes are determined based on the function combinations of two _RYR1_ variants. PharmCAT will prioritize the variants with Malignant Hyperthermia-associated function for individuals with more than two variants. If the input VCF file is not phased, PharmCAT assumes an individual to have two or more Malignant Hyperthermia-associated variants on different chromosomes. PharmCAT prioritizes variants with Malignant Hyperthermia-associated function when determining the RYR1 phenotype for an individual.

