---
title: How It Works
permalink: methods/
nav_order: 2
has_children: true
has_toc: false
---
# How PharmCAT Works

<img src="/images/pipeline-h.svg" alt="PharmCAT pipeline" class="img-pipeline" />

PharmCAT (Pharmacogenomics Clinical Annotation Tool) is a bioinformatics tool that analyzes genetic variants to predict
drug response and tailor medical treatment to an individual patient’s genetic profile.
It does this in two phases:

1. Processes VCF files from next generation sequencing (NGS) or genotyping methods and identifies pharmacogenomic (PGx)
   genotypes and infers haplotypes, typically called star alleles.
2. Uses the pharmacogene diplotypes (combination of maternal and paternal star alleles) to predict PGx phenotypes and
   reports the corresponding drug-prescribing recommendations from [CPIC guidelines](https://cpicpgx.org/guidelines/),
   [ClinPGx-annotated DPWG guidelines](https://www.clinpgx.org/page/dpwg) and
   [ClinPGx-annotated FDA-approved drug labels](https://www.clinpgx.org/page/drugLabelLegend).


### Phase 1 - Determining Alleles

The `VCF Preprocessor` module is responsible for normalizing the input VCF.
See [VCF Preprocessor](/using/VCF-Preprocessor) for details.

The `Named Allele Matcher` module is responsible for calling diplotypes for a sample based on a VCF file.  See:

* [Named Allele Matcher 101](/methods/NamedAlleleMatcher-101) for an overview
* [Named Allele Matcher 201](/methods/NamedAlleleMatcher-201) for details

Special cases and exceptions to how alleles are called are cataloged in
[Gene Definition Exceptions](/methods/Gene-Definition-Exceptions/).
For details on calling CYP2D6, please see [Calling CYP2D6](/using/Calling-CYP2D6).


### Phase 2 - Matching Recommendations

The `Phenotyper` module is responsible for translating diplotypes and other genotype information into gene-specific
phenotypes. For example, CYP2C19 `*2/*4` is a poor metabolizer.

The `Reporter` module is responsible for generating a report with genotype-specific expert-reviewed drug prescribing
recommendations for clinical decision support.

See [Matching Recommendations](/methods/Matching-Recommendations) for details on how the `Phenotyper` and `Reporter`
works. 


## Video Tutorial

Prefer video tutorials?
We have [an introduction to PharmCAT, modules, and reports](https://youtu.be/PjVdtMp8oRI?si=mRaiaU6EVEEd6dJL)
on YouTube.

<iframe width="560" height="315" src="https://www.youtube-nocookie.com/embed/PjVdtMp8oRI?si=LkWbiIN-L2A3gZeJ&amp;controls=0" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>



## Further Reading

Useful links for more information on drug recommendations. 

* [CPIC Site](https://cpicpgx.org)
* [CPIC Guideline Publications](https://cpicpgx.org/publications/)
* [ClinPGx Site](https://www.clinpgx.org)
* [List of ClinPGx clinical guideline annotations](https://www.clinpgx.org/guidelineAnnotations)
