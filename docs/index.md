---
title: Home
permalink: /
nav_order: 1
---

# PharmCAT:<br />Pharmacogenomics Clinical Annotation Tool

<span class="logoDiv">
<span class="logoDiv__logo">![PharmCAT logo](/images/pharmcat_logo.svg)</span>
<span class="logoDiv__links">
<span>[Download v{{ site.pharmcat_version }}](https://github.com/PharmGKB/PharmCAT/releases/latest){: .btn .btn-blue .umami--click--download-button }</span>
<span>[View on GitHub](https://github.com/PharmGKB/PharmCAT){: .btn }</span>
</span>
</span>

PharmCAT (Pharmacogenomics Clinical Annotation Tool) is a bioinformatics tool that analyzes genetic variants to predict
drug response and tailor medical treatment to an individual patient’s genetic profile. It does this in two phases:

1. Processes VCF files from next generation sequencing (NGS) or genotyping methods and identifies pharmacogenomic (PGx)
   genotypes and infers haplotypes, typically called star alleles.
2. Uses the pharmacogene diplotypes (combination of maternal and paternal star alleles) to predict PGx phenotypes and
   reports the corresponding drug-prescribing recommendations from [CPIC guidelines](https://cpicpgx.org/guidelines/),
   [ClinPGx-annotated DPWG guidelines](https://www.clinpgx.org/page/dpwg) and
   [ClinPGx-annotated FDA-approved drug labels](https://www.clinpgx.org/page/drugLabelLegend).

This is a very high-level example of this process:

<img src="/images/translation.svg" class="img-translation" alt="from variation to recommendation" />

For details, take a look at our documentation on [how PharmCAT works](/methods).

PharmCAT was developed in a collaboration between the PharmGKB (now [ClinPGx](https://www.clinpgx.org))
and the former [PGRN Statistical Analysis Resource (P-STAR)](https://ritchielab.org/pgrn-star/), with input from other
groups. The work was originally based on established guidelines from the
[Clinical Pharmacogenetics Implementation Consortium (CPIC)](https://cpicpgx.org). 

References:
- _Commentary:_ TE Klein, MD Ritchie. [PharmCAT: A Pharmacogenomics Clinical Annotation Tool](https://dx.doi.org/10.1002/cpt.928). Clinical Pharmacology & Therapeutics (2018) 104(1):19-22.
- _Methods paper:_ K Sangkuhl & M Whirl-Carrillo, et al. [Pharmacogenomics Clinical Annotation Tool (PharmCAT)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6977333). Clinical Pharmacology & Therapeutics (2020) 107(1):203-210.
- _Tutorial paper:_ B Li & K Sangkuhl et al. [How to Run the Pharmacogenomics Clinical Annotation Tool (PharmCAT)](https://ascpt.onlinelibrary.wiley.com/doi/10.1002/cpt.2790). Clinical Pharmacology & Therapeutics (2022).

PharmCAT is under active development.

Note: [PharmGKB is now ClinPGx](https://blog.clinpgx.org/pharmgkb-is-now-clinpgx/).
While the name and the website have changed, the code repositories and distribution channels have not.
It will remain this way for the foreseeable future, and you will still see references to PharmGKB in the documentation
while this is the case.

[Subscribe to PharmCAT Updates](https://pharmgkb.us10.list-manage.com/subscribe?u=c46dea014a68524407fdbffa1&id=d0d1ec73ab){: .btn .btn-blue}


## Documentation

See [Genes & Drugs](/Genes-Drugs) for a list of all genes and drugs supported by PharmCAT.
We have detailed documentation in [Gene Definition Exceptions](/methods/Gene-Definition-Exceptions) for genes that
require special handling.

We also have a whole section explaining [how PharmCAT works](/methods), including how it
[matches sample data to allele definitions](/methods/NamedAlleleMatcher-101)
and [matches genotypes to drug recommendations](/methods/Matching-Recommendations).

Finally, [learn how to run PharmCAT](/using) and the different components that make up the PharmCAT pipeline.
Please make sure to also read and understand PharmCAT's [VCF requirements](/using/VCF-Requirements).

Tutorial videos are available on the [PharmGKB YouTube channel](https://www.youtube.com/channel/UCnYHYK_5HD1Lt2N_B4FsTYQ).
These videos provide step-by-step demonstrations on running PharmCAT:

1. [An introduction to PharmCAT, modules, and reports](https://youtu.be/PjVdtMp8oRI?si=mRaiaU6EVEEd6dJL)
2. [How to run PharmCAT](https://youtu.be/d1IZPLOrPOE?si=LREY8RI-wz-5PoqN) - a hands-on example that walks through the setup and running PharmCAT


### Examples

Interested in seeing the kinds of reports that PharmCAT produces?

We have a collection of [sample reports and the data files that generated them](examples). 


## Contact

[Ask a Question](mailto:pharmcat@clinpgx.org){: .btn .btn-blue }
[Submit a Bug / Feature Request](https://github.com/PharmGKB/PharmCAT/issues/new){: .btn }


## Contributing

We welcome contributions from anyone interested in helping to improve PharmCAT.

If you notice a bug in PharmCAT, have an idea for a new feature, or just have a question about how PharmCAT works,
please [submit an issue on GitHub](https://github.com/PharmGKB/PharmCAT/issues).

If you want to make a code contribution to the project, please
[check out the code repo](https://github.com/PharmGKB/PharmCAT) and
[read the developer wiki](https://github.com/PharmGKB/PharmCAT/wiki) for information about our development process.


## License

PharmCAT is licensed under the [Mozilla Public License 2.0 (MPL-2.0)](https://github.com/PharmGKB/PharmCAT/blob/main/LICENSE).
