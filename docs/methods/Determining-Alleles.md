---
title: Determining Alleles
permalink: methods/determining-alleles
parent: Methods
nav_order: 3
---
# Determining Alleles

PharmCAT has an algorithm to match allele definitions to variant call data. The module in the PharmCAT responsible for this is the [Named Allele Matcher](/methods/NamedAlleleMatcher-101). The `Named Allele Matcher` was written specifically for use in PharmCAT but can be run independently.

PharmCAT can also incorporate output from other diplotype calling software (e.g. [Astrolabe](https://www.childrensmercy.org/childrens-mercy-research-institute/research-areas/genomic-medicine-center/data-and-software-resources/), [Stargazer](https://stargazer.gs.washington.edu/stargazerweb/index.html), [StellarPGx](https://github.com/SBIMB/StellarPGx) etc.). These tools are run separately from PharmCAT but the results can be fed into PharmCAT to add unmatched genes or supercede diplotypes matched by PharmCAT's `Named Allele Matcher`. See [Outside Call Format](/specifications/Outside-Call-Format) for details.

For details on calling CYP2D6, plesae see [Calling CYP2D6](/using/Calling-CYP2D6).
