---
title: Genes & Drugs
permalink: Genes-Drugs/
nav_order: 5
---

# PharmCAT Genes & Drugs

## Genes

The following tables list genes used by PharmCAT to find drug recommendation, along with the sources that use the gene.

### Genes PharmCAT will attempt to match

The `Named Allele Matcher` will search the given sample file for locations associated with these genes and attempt to match them to known allele definitions.

%s

### Genes handled by outside callers

These genes will not get allele matches from PharmCAT<sup>*</sup>. However, you can use an outside caller like
[Stargazer](https://stargazer.gs.washington.edu/stargazerweb/index.html) or
[StellarPGx](https://github.com/SBIMB/StellarPGx) to get diplotype calls and then supply that to PharmCAT for use in
matching recommendation data.  See [Outside Call Format](/using/Outside-Call-Format) for details.

%s

<sup>*</sup> Except for CYP2D6 if the requisite [research mode](/using/Running-PharmCAT#research-only-options) is enabled.


## Drugs

The following table lists drugs for which PharmCAT has recommendations for, along with their sources. 

%s
