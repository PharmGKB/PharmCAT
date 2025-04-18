## Overview

This directory contains the tests and test data for the PharmCAT VCF Preprocessor.

## Explanation of files

### 1. Test files for bcftools and the preprocessing scripts

`raw.vcf` is the primary VCF file.  From this file, the following files are generated by running `tests/prep.py`:

1. `raw.vcf.bgz`
2. `raw-p1.vcf.bgz`
3. `raw-p2.vcf.bgz`

The following files are the expected results of processing the `raw.vcf`:

1. `raw.preprocessed.vcf`
2. `raw.Sample_1.preprocessed.vcf`
3. `raw.Sample_2.preprocessed.vcf`

The remaining test files are used to test specific corner cases. 


### 2. Test files for the VCF preprocessing script - Performance

I tested the performance of the VCF preprocessing script, including run time, multi-sample VCF processing, etc.

The data was the 1000 Genomes Project sequences of Coriell samples with corresponding Genetic Testing Reference Materials Coordination Program (GeT-RM) sample characterization. This dataset was generated by Adam Lavertu, Ryan and Mark for the paper, [_Pharmacogenomics Clinical Annotation Tool (PharmCAT)_](https://doi.org/10.1002/cpt.1568); and is hosted on [the Stanford Digital Repository](https://purl.stanford.edu/rd572fp2219).

```
# Follow the script 01-03 downloaded from the Stanford Digital Repository. You will need to modify the codes, such as file paths.

# run VCF preprocessing
pharmcat_vcf_preprocessor \
-vcf PharmCAT_calling_pipeline-master/data/1kg_data/GeT-RM_sample_data/PGx.chrAllPGx.GRCh38.genotypes.20170504.vcf.gz \
-refFna GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
-refVcf pharmcat_positions_0.8.0_updated_06222021.vcf.gz
```



