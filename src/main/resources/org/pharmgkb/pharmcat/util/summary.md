---
title: Genes & Drugs
permalink: Genes-Drugs/
nav_order: 4
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
[StellarPGx](https://github.com/SBIMB/StellarPGx) to get CYP2D6 diplotype calls, or
[Optitype](https://github.com/FRED-2/OptiType) to get HLA calls, and then supply that to PharmCAT for use in
matching recommendation data.

See [Outside Call Format](/using/Outside-Call-Format) for formatting details and
[Calling CYP2D6(/using/Calling-CYP2D6) or [Calling HLA](/using/Calling-HLA) for how to obtain CYP2D6 or
HLA calls, respectively, using sequencing data.

%s

<sup>*</sup> Except for CYP2D6 if the requisite [research mode](/using/Running-PharmCAT#research-only-options) is enabled.


## Drugs

The following table lists drugs for which PharmCAT has recommendations for, along with their sources. 

%s
