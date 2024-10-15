# PharmCAT

[![PharmCAT CI](https://github.com/PharmGKB/PharmCAT/actions/workflows/ci-pharmcat.yml/badge.svg)](https://github.com/PharmGKB/PharmCAT/actions/workflows/ci-pharmcat.yml)
[![Preprocessor CI](https://github.com/PharmGKB/PharmCAT/actions/workflows/ci-preprocessor.yml/badge.svg)](https://github.com/PharmGKB/PharmCAT/actions/workflows/ci-preprocessor.yml)
[![codecov.io](https://codecov.io/github/PharmGKB/PharmCAT/coverage.svg?branch=main)](https://codecov.io/github/PharmGKB/PharmCAT?branch=development)


PharmCAT (Pharmacogenomics Clinical Annotation Tool) is a bioinformatics tool that analyzes genetic variants to predict
drug response and tailor medical treatment to an individual patient’s genetic profile. It does this in two phases:

1. Processes VCF files from next generation sequencing (NGS) or genotyping methods and identifies pharmacogenomic (PGx)
   genotypes and infers haplotypes, typically called star alleles.
2. Uses the pharmacogene diplotypes (combination of maternal and paternal star alleles) to predict PGx phenotypes and
   reports the corresponding drug-prescribing recommendations from [CPIC guidelines](https://cpicpgx.org/guidelines/),
   [PharmGKB-annotated DPWG guidelines](https://www.pharmgkb.org/page/dpwg) and
   [PharmGKB-annotated FDA-approved drug labels](https://www.pharmgkb.org/page/drugLabelLegend).

For more information:

* The [PharmCAT website](https://pharmcat.org) will have the latest documentation
* _Commentary:_ TE Klein, MD Ritchie.
  [PharmCAT: A Pharmacogenomics Clinical Annotation Tool](https://dx.doi.org/10.1002/cpt.928).
  Clinical Pharmacology & Therapeutics (2018) 104(1):19-22.
* _Methods paper:_ K Sangkuhl & M Whirl-Carrillo, et al.
  [Pharmacogenomics Clinical Annotation Tool (PharmCAT)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6977333).
  Clinical Pharmacology & Therapeutics (2020) 107(1):203-210.
* _Tutorial paper:_ B Li & K Sangkuhl et al.
  [How to Run the Pharmacogenomics Clinical Annotation Tool (PharmCAT)](https://ascpt.onlinelibrary.wiley.com/doi/10.1002/cpt.2790).
  Clinical Pharmacology & Therapeutics (2022).


## Status

PharmCAT is available for general use, but it is still under _active development_. 
New features, data updates, and bug fixes will be released.
Watch this repository or check the [releases](../../releases) page for new releases.

All technical requirements and documentation are available on [PharmCAT.org](https://pharmcat.org).

PharmCAT is managed at Stanford University & University of Pennsylvania (NHGRI U24HG013077).


## Contact

For technical questions or bug reports, [file an issue](https://github.com/PharmGKB/PharmCAT/issues).

For general questions about the PharmCAT project, contact [pharmcat@pharmgkb.org](mailto:pharmcat@pharmgkb.org).


## Liability

:warning: PharmCAT assumes no responsibility for any injury to person or damage to persons or property arising out of,
or related to any use of PharmCAT, or for any errors or omissions.
The user recognizes they are using PharmCAT at their own risk.
