---
parent: Using PharmCAT
title: Calling HLA
permalink: using/Calling-HLA/
nav_order: 9
---

# Calling HLA-A and B Alleles
PharmCAT supports HLA allele data with the following option, which is discussed in more detail below:
* Using a different tool to determine HLA alleles, which can be imported into PharmCAT to provide a translation of the phenotypes and recommendations.

## HLA Complications

### Problems with calling HLA from VCF

The HLA region is highly polymorphic and complex region filled with genes, retrotransposons, transposons, regulatory elements, and pseudogenes that make mapping using with [short sequence reads problematic](https://doi.org/10.1038/jhg.2008.5). The current human reference genomes are linear references (hg19 and hg38); therefore, [a single allele is designated as a reference allele](https://www.ebi.ac.uk/ipd/imgt/hla/help/genomics.html) to serve as part of the human genome. While efforts have been made to improve mapping the reads by including ALT contigs to the reference, many reads relating to HLA remain in the unmapped files. Therefore, calling HLA from a VCF is difficult given that most of the variation necessary to type an HLA allele is missing from the files. Unless you have targeted sequence data for HLA, your other options are either HLA imputation from a VCF file or calling from alignment files/raw sequence reads containing both mapped and unmapped reads.

### HLA imputation

We tested imputation accuracy using 258 in-house samples for which we had targeted amplicon HLA sequencing and low-pass whole genomes (~2.54X). The low-pass whole genomes were phased using GLIMPSE and phase-corrected using RFMix. The HLA region was imputed using the Michigan Imputation Server’s Four-digit Multi-ethnic HLA v1 (2021) panel, which uses the SNP2HLA program as part of the [HLA-TAPAS](https://github.com/immunogenomics/HLA-TAPAS) software to call HLA. The reference panel used contains [36,586 HLA haplotypes belonging to a multi-ethnic cohort](https://imputationserver.readthedocs.io/en/latest/reference-panels/#four-digit-multi-ethnic-hla-v1-2021). For our testing purposes, we looked at 4 important PGx HLA alleles (HLA-A*31:01, HLA-B*15:02, HLA-B*57:01, HLA-B*58:01). All 4 alleles passed QC as described in [Luo et al. 2021](https://doi.org/10.1038/s41588-021-00935-7). A summary of our results are as follows:

![Accuracy of imputed HLAs](/images/hla_imputation_accuracy.png)

| HLA Allele  | Sequence | Imputed |
|-------------| -------- | ------- |
| HLA-A*31:01 | 25 | 14 |
| HLA-B*15:02 | 3 | 2 |
| HLA-B*57:01 | 13 | 13 |
| HLA-B*58:01 | 8 | 7 |

**We do NOT recommend imputing HLA from a VCF**, even though imputation is possible from VCF. Our results indicate that HLA-B allele imputation was more accurate than the single HLA-A allele tested. Nonetheless, there are many factors that contribute to inaccuracy (e.g., marker density, an appropriate reference panel, etc.). Although [the methods and reference panels](https://doi.org/10.1038/tpj.2017.7) have drastically improved over the past few years, we still have reservations based on our own internal testing of the accuracy of allele typing. Instead, we recommend that you call HLA directly from either whole-exome or whole-genome data. 

## HLA from WES and WGS alignment files

For demonstration purposes, we will use the commonly used program _[Optitype](https://github.com/FRED-2/OptiType)_, available through the _Nextflow_ analysis pipeline [HLA Typing](https://nf-co.re/hlatyping) (which is also where StellarPGx is hosted). Optitype is a freely available software package used to call HLA Class I alleles. Please note that the references are based on the IMGT/HLA Release 3.14.0, which does not affect the HLA alleles we tested for PGx. Nextflow is currently working on building a solution to bring in [newer IMGT releases](https://nf-co.re/hlatyping/2.0.0/usage#hla-references). If you are looking to call HLA Class II alleles, please look at [HLA-LA](https://github.com/DiltheyLab/HLA-LA) but note that they only type HLA based on the protein binding domains and therefore are only able to output [G-group resolution](https://hla.alleles.org/alleles/g_groups.html). 

For our test, we used the [1000 genome HLA calls](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/) that are available through [The International Genome Sample Resource](https://www.internationalgenome.org/category/hla/) and described in [Abi-Rached et al, 2018](https://doi.org/10.1371/journal.pone.0206512). These HLA calls were generated using the commercial software PolyPheMe from WGS. We picked several samples (minimum 10 per allele) across different populations that had one of our alleles of interest: _HLA-A_*31:01, _HLA-B_*15:02, _HLA-B_*57:01, _HLA-B_*58:01.

We then accessed the 30X WGS CRAM (Compressed Reference-oriented Alignment Map) files through the Globus Endpoint. Since Optitype takes either FASTQ files (Raw Sequence files with Quality Scores) or BAM files (Binary Alignment Map) as input. Therefore, we converted the CRAM files to BAM files using Samtools. CRAM files are a lot smaller than BAM files but do not contain all the information in a BAM file and therefore require the FASTA file used in the alignment process. Converting a CRAM to BAM is easy, but be prepared for a much larger disk space footprint of the BAM file:

```commandline
samtools view -b -T <refgenome.fa> -o <output_file.bam> <input_file.cram>
```

For a 1000 genome use case scenario, such as HG00140.final.cram, we can convert it as:

```commandline
samtools view -b -T GRCh38_full_analysis_set_plus_decoy_hla.fa -o HG02420.final.bam HG02420.final.cram
```

You can visit the [samtools documentation](http://www.htslib.org/doc/samtools.html) if you would like to modify the commands or look into converting the CRAM file to FASTQ files. Once converted to BAM files, you can follow the [Full Sample sheet documentation](https://nf-co.re/hlatyping/2.0.0/usage#full-samplesheet) and [the usage docs](https://nf-co.re/hlatyping/2.0.0/usage) to see how to proceed to generate a sample list and submit it for Nextflow. We are not affiliated with either Nextflow or Optitype and we suggest you contact them for [questions or for troubleshooting](https://github.com/nf-core/hlatyping). 

Regardless of whether it is a single sample per sheet or multi-sample per sample sheet, Optitype will generate individual files for each sample. It will generate a plot showing the overall coverage of the allele and a TSV (tab-separated values) file with the HLA allele calls. Although not necessary, it is always a good idea to check and ensure you had good coverage across the HLA region. For our 4 alleles tested, our results indicated a 1:1 relationship for those specific alleles. Similarly, HLA-LA also had a 1:1 relationship for those specific alleles once we truncated the G-group resolution to a 2-field, but the G-group nomenclature also includes other HLA proteins that share the same nucleotide sequence for the protein binding domain (exons 2 and 3 for HLA class I and exon 2 only for HLA class II alleles). For example:

A*33:01:01G includes alleles such as A*33:01:01:01, A*33:01:01:02, A*33:01:01:03 which share the same protein, but it also contains A*33:220, A*33:222, A*33:228, A*33:234, which have different proteins, but share the same nucleotide sequence in exons 2 and 3 as the others.

Although we did not simulate coverage, a recent study by [Thuesen et al. 2022](https://doi.org/10.3389/fimmu.2022.987655), evaluated software performances of the most common HLA calling software, including Optiype, at different simulated coverage and degraded conditions such as those with  aDNA. While performance was relatively stable across different coverages for Whole Exome Sequence Data and Optitype, we still recommend working with as high coverage as possible. 

### Working with HLA calls from WES/WGS in PharmCAT

PharmCAT supports incorporating results from your favorite HLA programs for phenotype translations through “[outside calls](https://pharmcat.org/using/Outside-Call-Format/)”.

To incorporate the outside calls, you would run PharmCAT as you normally would, and add the -po flag which signals the program to look for an external call file. Please note that PharmCAT requires that your HLA call have only two fields for phenotype translations. Therefore, if your external calls have more than two fields, you should truncate your output to only two fields. For more information on HLA nomenclature and what each field means in HLA, please visit the page for [Nomenclature for Factors of the HLA System page](https://hla.alleles.org/nomenclature/naming.html). For example:

```commandline
java -jar pharmcat.jar -vcf test.vcf -po test_sample_hla.txt
```

The format of this file (test_sample_hla.txt) should be:
```text
HLA-A	*32:01/*68:03
HLA-B	*07:02/*35:01
```
Where the first column is the gene name, followed by the allele calls.

For more information:
- see [Running PharmCAT](https://pharmcat.org/using/Running-PharmCAT#phenotyper) for details on the -po flag

- see [Outside Call Format](https://pharmcat.org/using/Outside-Call-Format) for details on the outside call file

### Formatting Optitype output for PharmCAT

We recommend calling HLA from Whole Genome Sequencing or Whole Exome Sequencing (WGS or WES) from BAM or FASTQ files using software similar to Optitype. We also recommend that you use Optitype through the [Nextflow pipeline](https://nf-co.re/hlatyping) for reproducibility and ease of use and installation. This is especially ideal if you are already using the pipeline to call CYP2D6 using StellarPGx.

While this tutorial is specific to Optitype, it can be modified for integrating HLA calls from other programs. When you run Optitype, whether a single sample or multiple sample run, the program outputs a single file TSV file per individual, which looks like:

```text
	A1	A2	B1	B2	C1	C2	Reads	Objective
0	A*32:01	A*68:03	B*07:02	B*35:01	C*07:02	C*07:02	10191.0	9915.832999999959
```

Each column represents a separate HLA allele call. To convert this into a format digestible by PharmCAT, you must reorder the TSV file. This can be done using a simple Python script. The following is a simple script to parse Optitype’s output and output a file format usable for PharmCAT. Currently, we only support HLA A and B calls in PharmCAT, therefore only those two genes are output.

```python
#!/usr/bin/env python

import pandas as pd
import sys
import os

# Prompt for the Sample Name, which will be used as the file name and input file path to the TSV file
sample_name = sys.argv[1]
input_path = sys.argv[2]

# Load TSV file as a dataframe
df = pd.read_csv(input_path, sep="\t")

# Reformat Data
a_values = df.iloc[0, 1:3].str.replace("[AB](?=\*)", "", regex=True).str.cat(sep="/")
b_values = df.iloc[0, 3:5].str.replace("[AB](?=\*)", "", regex=True).str.cat(sep="/")

# Create new dataframe with reformatted data
new_df = pd.DataFrame([["HLA-A", a_values], ["HLA-B", b_values]])

# Export the data frame into a Tab-delimited file
new_df.to_csv(f"{sample_name}_HLA.txt", sep="\t", index=False, header=False)
```

Below is a usage example of how it would be called from through the command path: 

```shell
python3 Optitype_to_PharmCAT.py <sample_ID> <input file path or input file> 
python3 Optitype_to_PharmCAT.py test_sample test_sample_results.tsv
```

The output file (test_sample_hla.txt) should look like:

```text
HLA-A	*32:01/*68:03
HLA-B	*07:02/*35:01
```

This script can be used in conjunction with bash scripting to automate the processing of multiple sample files, where you can use the sample names as variables.
We are not responsible or can guarantee how these programs will perform and should primarily be used for research purposes. We suggest you consult a clinical laboratory if you need clinical-grade HLA typing, particularly high-resolution typing based on NGS. Generally, clinical laboratories can provide an [HML file](https://www.sciencedirect.com/science/article/pii/S0198885915004346?via%3Dihub), and we will explore how to bring that data into PharmCAT in the future.