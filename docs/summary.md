---
title: PharmCAT Data Summary
permalink: summary/
---

# PharmCAT Data Summary

## Genes

The following genes are used by PharmCAT to find drug recommendation data.

### Genes PharmCAT will attempt to match

The `NamedAlleleMatcher` will search the given sample file for locations associated with these genes and attempt to match them to known allele definitions. Each of these genes will appear in the "Genotype Summary" section of the final output report.

- CACNA1S
- CFTR
- CYP2B6
- CYP2C19
- CYP2C9
- CYP3A5
- CYP4F2
- DPYD
- IFNL3
- NUDT15
- RYR1
- SLCO1B1
- TPMT
- UGT1A1
- VKORC1

### Genes handled by outside callers

These genes will not get allele matches from PharmCAT. However, you can use an outside caller like [Astrolabe](https://www.childrensmercy.org/childrens-mercy-research-institute/research-areas/genomic-medicine-center/data-and-software-resources/) or [Stargazer](https://stargazer.gs.washington.edu/stargazerweb/index.html) to get diplotype calls and then supply that to PharmCAT for use in matching recommendation data.

- CYP2D6


## Drugs

The following drugs have been read from CPIC and will have "recommendation" sections in the final output report.

- amitriptyline
- atazanavir
- atomoxetine
- azathioprine
- capecitabine
- celecoxib
- citalopram
- clomipramine
- clopidogrel
- codeine
- desflurane
- desipramine
- dexlansoprazole
- doxepin
- efavirenz
- enflurane
- escitalopram
- fluorouracil
- flurbiprofen
- fluvoxamine
- halothane
- hydrocodone
- ibuprofen
- imipramine
- isoflurane
- ivacaftor
- lansoprazole
- lornoxicam
- meloxicam
- mercaptopurine
- methoxyflurane
- nortriptyline
- omeprazole
- ondansetron
- pantoprazole
- paroxetine
- peginterferon alfa-2a
- peginterferon alfa-2b
- piroxicam
- sertraline
- sevoflurane
- simvastatin
- succinylcholine
- tacrolimus
- tamoxifen
- tenoxicam
- thioguanine
- tramadol
- trimipramine
- tropisetron
- voriconazole
- warfarin
