---
title: Genes & Drugs
permalink: Genes-Drugs/
nav_order: 5
---

# PharmCAT Genes & Drugs

## Genes

The following genes are used by PharmCAT to find drug recommendation data.

### Genes PharmCAT will attempt to match

The `Named Allele Matcher` will search the given sample file for locations associated with these genes and attempt to match them to known allele definitions. Each of these genes will appear in the "Genotype Summary" section of the final output report.

%s

### Genes handled by outside callers

These genes will not get allele matches from PharmCAT. However, you can use an outside caller like [Astrolabe](https://www.childrensmercy.org/childrens-mercy-research-institute/research-areas/genomic-medicine-center/data-and-software-resources/) or [Stargazer](https://stargazer.gs.washington.edu/stargazerweb/index.html) to get diplotype calls and then supply that to PharmCAT for use in matching recommendation data.

%s


## Drugs

The following drugs are included in the PharmCAT report "Prescribing Recommendations" section. The sources for the drug prescribing information are in brackets after the drug name.

%s
