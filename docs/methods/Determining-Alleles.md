---
title: Determining Alleles
permalink: methods/determining-alleles/
parent: Methods
---
# Determining Alleles

PharmCAT has an algorithm to match allele definitions to variant call data. This algorithm in the PharmCAT codebase is called [NamedAlleleMatcher](/methods/NamedAlleleMatcher-101). The `NamedAlleleMatcher` was written specifically for use in PharmCAT but can be run independantly.

PharmCAT can also incorporate output from other diplotype caller software (e.g. [Astrolabe](https://www.childrensmercy.org/childrens-mercy-research-institute/research-areas/genomic-medicine-center/data-and-software-resources/), [Stargazer](https://stargazer.gs.washington.edu/stargazerweb/index.html), etc...). These tools are run separately from PharmCAT but the results can be fed into PharmCAT to add unmatched genes or supercede diplotypes matched by PharmCAT's `NamedAlleleMatcher`.

Currently, CYP2D6 diplotypes must be supplied by external caller software and the `NamedAlleleMatcher` is the default for all other gene diplotype matching.
