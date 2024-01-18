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

| Gene | CPIC | PharmGKB-DPWG |
| :--- | :---: | :---: |
| [ABCG2](/Phenotypes-List#abcg2) | :heavy_check_mark: | :heavy_check_mark: |
| [CACNA1S](/Phenotypes-List#cacna1s) | :heavy_check_mark: |  |
| [CFTR](/Phenotypes-List#cftr) | :heavy_check_mark: |  |
| [CYP2B6](/Phenotypes-List#cyp2b6) | :heavy_check_mark: | :heavy_check_mark: |
| [CYP2C19](/Phenotypes-List#cyp2c19) | :heavy_check_mark: | :heavy_check_mark: |
| [CYP2C9](/Phenotypes-List#cyp2c9) | :heavy_check_mark: | :heavy_check_mark: |
| [CYP3A4](/Phenotypes-List#cyp3a4) |  | :heavy_check_mark: |
| [CYP3A5](/Phenotypes-List#cyp3a5) | :heavy_check_mark: | :heavy_check_mark: |
| [CYP4F2](/Phenotypes-List#cyp4f2) | :heavy_check_mark: |  |
| [DPYD](/Phenotypes-List#dpyd) | :heavy_check_mark: | :heavy_check_mark: |
| [G6PD](/Phenotypes-List#g6pd) | :heavy_check_mark: |  |
| [IFNL3](/Phenotypes-List#ifnl3) | :heavy_check_mark: |  |
| [NUDT15](/Phenotypes-List#nudt15) | :heavy_check_mark: | :heavy_check_mark: |
| [RYR1](/Phenotypes-List#ryr1) | :heavy_check_mark: |  |
| [SLCO1B1](/Phenotypes-List#slco1b1) | :heavy_check_mark: | :heavy_check_mark: |
| [TPMT](/Phenotypes-List#tpmt) | :heavy_check_mark: | :heavy_check_mark: |
| [UGT1A1](/Phenotypes-List#ugt1a1) | :heavy_check_mark: | :heavy_check_mark: |
| [VKORC1](/Phenotypes-List#vkorc1) | :heavy_check_mark: | :heavy_check_mark: |


### Genes handled by outside callers

These genes will not get allele matches from PharmCAT<sup>*</sup>. However, you can use an outside caller like
[Stargazer](https://stargazer.gs.washington.edu/stargazerweb/index.html) or
[StellarPGx](https://github.com/SBIMB/StellarPGx) to get CYP2D6 diplotype calls, or
[Optitype](https://github.com/FRED-2/OptiType) to get HLA calls, and then supply that to PharmCAT for use in
matching recommendation data.

See [Outside Call Format](/using/Outside-Call-Format) for formatting details and
[Calling CYP2D6](/using/Calling-CYP2D6) or [Calling HLA](/using/Calling-HLA) for how to obtain CYP2D6 or
HLA calls, respectively, using sequencing data.

| Gene | CPIC | PharmGKB-DPWG |
| :--- | :---: | :---: |
| [CYP2D6](/Phenotypes-List#cyp2d6) | :heavy_check_mark: | :heavy_check_mark: |
| [F5](/Phenotypes-List#f5) |  | :heavy_check_mark: |
| [HLA-A](/Phenotypes-List#hla-a) | :heavy_check_mark: | :heavy_check_mark: |
| [HLA-B](/Phenotypes-List#hla-b) | :heavy_check_mark: | :heavy_check_mark: |
| [MT-RNR1](/Phenotypes-List#mt-rnr1) | :heavy_check_mark: |  |


<sup>*</sup> Except for CYP2D6 if the requisite [research mode](/using/Running-PharmCAT#research-only-options) is enabled.


## Drugs

The following table lists drugs for which PharmCAT has recommendations for, along with their sources. 

| Drug | CPIC | PharmGKB-DPWG |
| :--- | :---: | :---: |
| abacavir | :heavy_check_mark: | :heavy_check_mark: |
| acenocoumarol |  | :heavy_check_mark: |
| allopurinol | :heavy_check_mark: | :heavy_check_mark: |
| amikacin | :heavy_check_mark: |  |
| amitriptyline | :heavy_check_mark: | :heavy_check_mark: |
| aripiprazole |  | :heavy_check_mark: |
| atazanavir | :heavy_check_mark: |  |
| atomoxetine | :heavy_check_mark: | :heavy_check_mark: |
| atorvastatin | :heavy_check_mark: | :heavy_check_mark: |
| azathioprine | :heavy_check_mark: | :heavy_check_mark: |
| brexpiprazole |  | :heavy_check_mark: |
| capecitabine | :heavy_check_mark: | :heavy_check_mark: |
| carbamazepine | :heavy_check_mark: | :heavy_check_mark: |
| celecoxib | :heavy_check_mark: |  |
| citalopram | :heavy_check_mark: | :heavy_check_mark: |
| clomipramine | :heavy_check_mark: | :heavy_check_mark: |
| clopidogrel | :heavy_check_mark: | :heavy_check_mark: |
| codeine | :heavy_check_mark: | :heavy_check_mark: |
| dapsone | :heavy_check_mark: |  |
| desflurane | :heavy_check_mark: |  |
| desipramine | :heavy_check_mark: |  |
| dexlansoprazole | :heavy_check_mark: |  |
| doxepin | :heavy_check_mark: | :heavy_check_mark: |
| efavirenz | :heavy_check_mark: | :heavy_check_mark: |
| eliglustat |  | :heavy_check_mark: |
| enflurane | :heavy_check_mark: |  |
| escitalopram | :heavy_check_mark: | :heavy_check_mark: |
| flecainide |  | :heavy_check_mark: |
| flucloxacillin |  | :heavy_check_mark: |
| flucytosine |  | :heavy_check_mark: |
| fluorouracil | :heavy_check_mark: | :heavy_check_mark: |
| flurbiprofen | :heavy_check_mark: |  |
| fluvastatin | :heavy_check_mark: |  |
| fluvoxamine | :heavy_check_mark: |  |
| fosphenytoin | :heavy_check_mark: |  |
| gentamicin | :heavy_check_mark: |  |
| haloperidol |  | :heavy_check_mark: |
| halothane | :heavy_check_mark: |  |
| hormonal contraceptives for systemic use |  | :heavy_check_mark: |
| hydrocodone | :heavy_check_mark: |  |
| ibuprofen | :heavy_check_mark: |  |
| imipramine | :heavy_check_mark: | :heavy_check_mark: |
| irinotecan |  | :heavy_check_mark: |
| isoflurane | :heavy_check_mark: |  |
| ivacaftor | :heavy_check_mark: |  |
| kanamycin | :heavy_check_mark: |  |
| lamotrigine |  | :heavy_check_mark: |
| lansoprazole | :heavy_check_mark: | :heavy_check_mark: |
| lornoxicam | :heavy_check_mark: |  |
| lovastatin | :heavy_check_mark: |  |
| meloxicam | :heavy_check_mark: |  |
| mercaptopurine | :heavy_check_mark: | :heavy_check_mark: |
| methoxyflurane | :heavy_check_mark: |  |
| methylene blue | :heavy_check_mark: |  |
| metoprolol |  | :heavy_check_mark: |
| nitrofurantoin | :heavy_check_mark: |  |
| nortriptyline | :heavy_check_mark: | :heavy_check_mark: |
| omeprazole | :heavy_check_mark: | :heavy_check_mark: |
| ondansetron | :heavy_check_mark: |  |
| oxcarbazepine | :heavy_check_mark: | :heavy_check_mark: |
| pantoprazole | :heavy_check_mark: | :heavy_check_mark: |
| paromomycin | :heavy_check_mark: |  |
| paroxetine | :heavy_check_mark: | :heavy_check_mark: |
| peginterferon alfa-2a | :heavy_check_mark: |  |
| peginterferon alfa-2b | :heavy_check_mark: |  |
| pegloticase | :heavy_check_mark: |  |
| phenprocoumon |  | :heavy_check_mark: |
| phenytoin | :heavy_check_mark: | :heavy_check_mark: |
| pimozide |  | :heavy_check_mark: |
| piroxicam | :heavy_check_mark: |  |
| pitavastatin | :heavy_check_mark: |  |
| plazomicin | :heavy_check_mark: |  |
| pravastatin | :heavy_check_mark: |  |
| primaquine | :heavy_check_mark: |  |
| propafenone |  | :heavy_check_mark: |
| quetiapine |  | :heavy_check_mark: |
| rasburicase | :heavy_check_mark: |  |
| ribavirin | :heavy_check_mark: |  |
| risperidone |  | :heavy_check_mark: |
| rosuvastatin | :heavy_check_mark: |  |
| sertraline | :heavy_check_mark: | :heavy_check_mark: |
| sevoflurane | :heavy_check_mark: |  |
| simvastatin | :heavy_check_mark: | :heavy_check_mark: |
| siponimod |  | :heavy_check_mark: |
| streptomycin | :heavy_check_mark: |  |
| succinylcholine | :heavy_check_mark: |  |
| tacrolimus | :heavy_check_mark: | :heavy_check_mark: |
| tafenoquine | :heavy_check_mark: |  |
| tamoxifen | :heavy_check_mark: | :heavy_check_mark: |
| tegafur |  | :heavy_check_mark: |
| tenoxicam | :heavy_check_mark: |  |
| thioguanine | :heavy_check_mark: | :heavy_check_mark: |
| tobramycin | :heavy_check_mark: |  |
| toluidine blue | :heavy_check_mark: |  |
| tramadol | :heavy_check_mark: | :heavy_check_mark: |
| trimipramine | :heavy_check_mark: |  |
| tropisetron | :heavy_check_mark: |  |
| venlafaxine | :heavy_check_mark: | :heavy_check_mark: |
| voriconazole | :heavy_check_mark: | :heavy_check_mark: |
| vortioxetine | :heavy_check_mark: |  |
| warfarin | :heavy_check_mark: | :heavy_check_mark: |
| zuclopenthixol |  | :heavy_check_mark: |

