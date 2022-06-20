---
parent: Methods
title: PharmCAT Modules
permalink: methods/PharmCAT-Modules/
nav_order: 2
render_with_liquid: false
---
# PharmCAT Modules

The PharmCAT tool is really a pipeline comprised of multiple modules:

The `Named AlleleMatcher` module is responsible for calling diplotypes for a sample based on a VCF file.  For an overview on how this works, see [Named Allele Matcher 101](/specifications/NamedAlleleMatcher-101)

The `Phenotyper` module is reponsible for translating diplotypes and other genotype information into gene-specific phenotypes. For example, CYP2C19 `*2/*4` is a poor metabolizer.

The `Reporter` module is responsible for generating a report with genotype-specific expert-reviewed drug prescribing recommendations for clinical decision support.

![process diagram](/images/flowchart.png)
