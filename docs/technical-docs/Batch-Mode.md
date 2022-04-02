---
title: Batch Processing Multiple Samples
permalink: technical-docs/batch-mode/
parent: Using PharmCAT
nav_order: 5
nav_exclude: true
---

# Batch Annotating Multiple Samples

As of March 2022, PharmCAT only takes a single-sample VCF. If a multi-sample VCF is provided, only the first sample will be annotated. However, noteworthy, PharmCAT can generate a PGx report in the matter of seconds for a preprocessed VCF. The fast runtime of the PharmCAT enables the users to batch-annotate VCFs, with help from as few as one simple additional script, from the scale of a few dozen samples to a biobank-scale of cohort in an efficient manner.

## Example using the Stanford Sherlock HPC

This part of the tutorial uses the VCFs on multiple GeT-RM samples [1] generated from the 30x high-coverage WGS sequencing performed by the New York Genome Center (NYGC) [2]. The job was run on the Stanford Sherlock, a High-Performance Computing (HPC) cluster which uses Slurm as an open-source resource manager and job scheduler. The scripts to run the PharmCAT on the Stanford Sherlock is written in shell programming language, which can be easily adapted to another HPC environment supported by a different job scheduler.

### Preprocessing the data

We presume that you have familiarized yourself with using the PharmCAT VCF preprocessor. If not, please check out [Preprocessing VCF Files for PharmCAT](). 

A separate PharmCAT tutorial GitHub repository provides exemplary commands working with real genetic data sets. The following cases provide quick exemplary codes for different types of input VCFs. For full commands and a thorough walk-through, please check out [the PharmCAT tutorial GitHub repo](https://github.com/PharmGKB/PharmCAT-tutorial), specifically _**src/02_VCF_preprocessing.sh**_.

### Preprocessing case 1 - single-sample VCFs

If your genetic data is already stored in single-sample VCFs, you are one step closer to running the PharmCAT. We recommend the users to still run the PharmCAT VCF preprocessor to ensure the appropriate variant representation format, which can be achieved by the following command; and the complete commands can be found in [PharmCAT-tutorial/src/02_VCF_preprocessing.sh](https://github.com/PharmGKB/PharmCAT-tutorial/blob/main/src/).

```commandline
# a file listing all your single-sample VCFs, say single_sample_vcf_list.txt
# the following command shows the content of the single_sample_vcf_list.txt
$ cat single_sample_vcf_list.txt
PharmCAT_tutorial_get-rm_wgs_30x_grch38.NA18526.vcf.gz
PharmCAT_tutorial_get-rm_wgs_30x_grch38.NA18565.vcf.gz
PharmCAT_tutorial_get-rm_wgs_30x_grch38.NA18861.vcf.gz

# run the preprocessor iteratively across single-sample VCFs
for SINGLE_VCF in $(cat single_sample_vcf_list.txt)
do 
    # run the PharmCAT VCF preprocessor
    python3 PharmCAT_VCF_Preprocess.py --input_vcf "$SINGLE_VCF"
done
```

### Preprocessing case 2 - multi-sample VCF

Population- or biobank-scale VCFs most likely come in multi-sample format. The PharmCAT VCF preprocessor is designed to help the users with this case and produce multiple single-sample VCFs that PharmCAT requires. The simplest command to preprocess a multi-sample VCF is as following but the complete example can be found in [PharmCAT-tutorial/src/02_VCF_preprocessing.sh](https://github.com/PharmGKB/PharmCAT-tutorial/blob/main/src/02_VCF_preprocessing.sh).

```commandline
# run the PharmCAT VCF preprocessor for a multi-sample VCF
python3 PharmCAT_VCF_Preprocess.py --input_vcf <multi_sample.vcf.gz>
```

### Preprocessing case 3 - multi-sample VCF divided by chromosome or into consecutive genetic blocks

As sometimes seen with large-scale genetic studies, the genetic data may be divided into multiple by-chromosome VCFs or VCFs with consecutive genetic blocks. The PharmCAT VCF preprocessor can manage this type of genetic data sets by taking a list of VCFs as the input. The complete example can be found in [PharmCAT-tutorial/src/02_VCF_preprocessing.sh](https://github.com/PharmGKB/PharmCAT-tutorial/blob/main/src/02_VCF_preprocessing.sh)

```commandline
# run the PharmCAT VCF preprocessor for multiple VCFs with non-overlapping genetic regions of the same cohort
python3 PharmCAT_VCF_Preprocess.py --input_list <input_vcf_list.txt>
```

##  Running PharmCAT

Assuming the users have run the PharmCAT VCF preprocessor to generate multiple single-sample VCFs which is named such as “pharmcat_ready_vcf.<sample_id>.vcf”. Use the following command to batch annotate multiple VCFs using PharmCAT. A full example can be found at [PharmCAT-tutorial/src/03_PharmCAT.sh](https://github.com/PharmGKB/PharmCAT-tutorial/blob/main/src/03_PharmCAT.sh) in the PharmCAT tutorial GitHub repository.

```commandline
bcftools query -l multi_sample.vcf.gz > sample_list.txt
# run PharmCAT for all samples
for SINGLE_SAMPLE in $(cat sample_list.txt)
do
    # always use the latest PharmCAT
    java  -jar  <path_to_the_latest_pharmcat_jar>  \
    -vcf  pharmcat_ready_vcf."$SINGLE_SAMPLE".vcf  \
    -o  <output_dir>  -f   pharmcat."$SINGLE_SAMPLE"
done
```

The output is a set of PGx reports in HTML format named as “pharmcat.<sample_id>.html”.

## Batch outside calls

To incorporate outside PGx calls with PharmCAT for multiple samples, the users have to prepare the outside PGx calls of each sample in separate files and supply the individual file as outside calls to the PharmCAT.

## Running individual components

PharmCAT users can run the individual components of the PharmCAT on multiple samples in a similar manner to how a user runs the whole PharmCAT. This behavior is desirable by the users who are interested in understanding population PGx and specifically obtaining PGx frequencies (named alleles, diplotypes, or metabolizer phenotypes) in their cohort. To achieve such statistics, one can run the following commands substituted with specific content. More details and a tutorial can be found at [PharmCAT-tutorial/src/03_PharmCAT.sh](https://github.com/PharmGKB/PharmCAT-tutorial/blob/main/src/03_PharmCAT.sh)

### Named Allele Matcher
```commandline
for SINGLE_SAMPLE in "$(cat <sample_list.txt>)"
do
    java  -cp  <path_to_the_latest_pharmcat_jar>  org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher  \
    -vcf  pharmcat_ready_vcf."$SINGLE_SAMPLE".vcf   \
    -json  pharmcat_named_allele_matcher."$SINGLE_SAMPLE".json
done
```

### Phenotyper - use the Named Allele Matcher JSON data
```commandline
for SINGLE_SAMPLE in "$(cat <sample_list.txt>)"
do
    java -cp  <path_to_the_latest_pharmcat_jar>  org.pharmgkb.pharmcat.phenotype.Phenotyper  \
    -c pharmcat_named_allele_matcher."$SINGLE_SAMPLE".json  \
    -f pharmcat_phenotyper."$SINGLE_SAMPLE".json
done
```

### Phenotyper - use VCF directly
```commandline
for SINGLE_SAMPLE in "$(cat <sample_list.txt>)"
do
    java -cp  <path_to_the_latest_pharmcat_jar>   org.pharmgkb.pharmcat.phenotype.Phenotyper   \
    -vcf  pharmcat_ready_vcf."$SINGLE_SAMPLE".vcf  \
    -f  pharmcat_phenotyper."$SINGLE_SAMPLE".json
done
```

# Reporter
```commandline
for SINGLE_SAMPLE in "$(cat <sample_list.txt>)"
do
    java  -cp  <path_to_the_latest_pharmcat_jar>  org.pharmgkb.pharmcat.reporter.Reporter  \
    -p pharmcat_phenotyper."$SINGLE_SAMPLE".json  \
    -o pharmcat_reporter."$SINGLE_SAMPLE".html  \
    -j pharmcat_reporter."$SINGLE_SAMPLE".json  \
    -t  'Report for '"$SINGLE_SAMPLE"
done
```

These commands yield Named Allele Matcher, Phenotyper, and Reporter results for each individual separately in the format of JSON files. We encourage the users to explore and perform data analysis using the rich content in these JSON files, which can be easily achieved using auxiliary JSON libraries and data analysis packages in R or python.

## Extracting PharmCAT JSON content into TSV

We also provide accessory R scripts that organize and extract the content from the Named Allele Matcher or Phenotyper JSON outputs into tab-separated values (TSV) files. The accessory result-organizing R scripts can be found in the [PharmCAT-tutorial/src/](https://github.com/PharmGKB/PharmCAT-tutorial/tree/main/src/). The commands are as follows:

### Extracting the PharmCAT Named Allele Matcher JSON data into a tabular file
```commandline
SCRIPT_PATH=src/organize_pharmcat_named_allele_matcher_results.R
MATCHER_DIR=results/pharmcat_named_allele_matcher/
MATCHER_PATTERN=pharmcat_named_allele_matcher*json
PROJECT_DIR="$PWD"
Rscript  “$SCRIPT_PATH” \
--input-dir "$MATCHER_DIR" \
--input-file-pattern "$MATCHER_PATTERN" \
--output-dir "$PROJECT_DIR"
```

### Extracting the PharmCAT Phenotyper JSON data into a tabular file

```commandline
SCRIPT_PATH=src/organize_pharmcat_phenotyper_results.R
PHENOTYPER_DIR=results/pharmcat_phenotyper/
PHENOTYPER_PATTERN=pharmcat_phenotyper*json
PROJECT_DIR="$PWD"
Rscript  "$SCRIPT_PATH" \
--input-dir "$PHENOTYPER_DIR" \
--input-file-pattern "$PHENOTYPER_PATTERN" \
--output-dir "$PROJECT_DIR"
```

You can find another example, applied to Penn Medicine Biobank data, about how to convert JSON data into a CSV file using python on the PharmCAT wiki.

## Reference
1. Pratt, V. M. et al. Characterization of 137 Genomic DNA Reference Materials for 28 Pharmacogenetic Genes: A GeT-RM Collaborative Project. J Mol Diagn 18, 109–123 (2016).
2. Byrska-Bishop, M. et al. High coverage whole genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios. 2021.02.06.430068 (2021) doi:10.1101/2021.02.06.430068.
