---
parent: Methods
title: Determining Alleles
permalink: methods/Determining-Alleles/
nav_order: 3
---
# Determining Alleles

PharmCAT has an algorithm to match allele definitions to variant call data. The module in the PharmCAT responsible for this is the [Named Allele Matcher](NamedAlleleMatcher-101). The `Named Allele Matcher` was written specifically for use in PharmCAT but can be run independently.

PharmCAT can also incorporate output from other diplotype calling software (e.g. [Stargazer](https://stargazer.gs.washington.edu/stargazerweb/index.html) and [StellarPGx](https://github.com/SBIMB/StellarPGx)). These tools are run separately from PharmCAT but the results can be fed into PharmCAT to add unmatched genes or supercede diplotypes matched by PharmCAT's `Named Allele Matcher`. See [Outside Call Format](/using/Outside-Call-Format) for details.

For details on calling CYP2D6, plesae see [Calling CYP2D6](/using/Calling-CYP2D6).
