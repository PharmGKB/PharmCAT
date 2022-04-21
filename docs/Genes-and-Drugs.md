---
title: Genes & Drugs
permalink: genes-drugs/
nav_order: 5
---

# PharmCAT Genes & Drugs

## Genes

The following genes are used by PharmCAT to find drug recommendation data.

### Genes PharmCAT will attempt to match

The `NamedAlleleMatcher` will search the given sample file for locations associated with these genes and attempt to match them to known allele definitions. Each of these genes will appear in the "Genotype Summary" section of the final output report.

- ABCG2 (2 alleles)
- CACNA1S (3 alleles)
- CFTR (40 alleles)
- CYP2B6 (35 alleles)
- CYP2C19 (34 alleles)
- CYP2C9 (85 alleles)
- CYP3A5 (6 alleles)
- CYP4F2 (3 alleles)
- DPYD (83 alleles)
- IFNL3 (2 alleles)
- NUDT15 (20 alleles)
- RYR1 (49 alleles)
- SLCO1B1 (42 alleles)
- TPMT (46 alleles)
- UGT1A1 (9 alleles)
- VKORC1 (2 alleles)

### Genes handled by outside callers

These genes will not get allele matches from PharmCAT. However, you can use an outside caller like [Astrolabe](https://www.childrensmercy.org/childrens-mercy-research-institute/research-areas/genomic-medicine-center/data-and-software-resources/) or [Stargazer](https://stargazer.gs.washington.edu/stargazerweb/index.html) to get diplotype calls and then supply that to PharmCAT for use in matching recommendation data.

- CYP2D6 (149 alleles)
- G6PD (187 alleles)
- MT-RNR1 (25 alleles)


## Drugs

The following drugs have been read from CPIC and will have "recommendation" sections in the final output report.

- abacavir
- allopurinol
- amikacin
- amitriptyline
- atazanavir
- atomoxetine
- atorvastatin
- azathioprine
- capecitabine
- carbamazepine
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
- fluvastatin
- fluvoxamine
- fosphenytoin
- gentamicin
- halothane
- hydrocodone
- ibuprofen
- imipramine
- isoflurane
- ivacaftor
- kanamycin
- lansoprazole
- lornoxicam
- lovastatin
- meloxicam
- mercaptopurine
- methoxyflurane
- nortriptyline
- omeprazole
- ondansetron
- oxcarbazepine
- pantoprazole
- paromomycin
- paroxetine
- peginterferon alfa-2a
- peginterferon alfa-2b
- phenytoin
- piroxicam
- pitavastatin
- plazomicin
- pravastatin
- rasburicase
- rosuvastatin
- sertraline
- sevoflurane
- simvastatin
- streptomycin
- succinylcholine
- tacrolimus
- tamoxifen
- tenoxicam
- thioguanine
- tobramycin
- tramadol
- trimipramine
- tropisetron
- voriconazole
- warfarin
