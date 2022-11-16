---
parent: Using PharmCAT
title: Running the VCF Preprocessor
permalink: using/VCF-Preprocessor/
nav_order: 4
render_with_liquid: false
---
# PharmCAT VCF Preprocessor

The PharmCAT VCF preprocessor is a script that can preprocess VCF files for PharmCAT.

This tool will:

1. Strip out PGx positions that PharmCAT does not care about.
2. Break down a multi-sample VCF to multiple single-sample VCF files.
3. Automatically download the necessary Human Reference Genome Sequence FASTA and index files from the NIH FTP site if they are not provided.
4. __Perform VCF normalization__ - a standardization process that turns VCF into a parsimonious, left-aligned variant representation format (as discussed in [Unified Representation of Genetic Variants](https://doi.org/10.1093/bioinformatics/btv112) by Tan, Abecasis, and Kang).
5. Normalize the multiallelic variant representation to PharmCAT's expectation.
6. Process a subset of samples if a sample file is provided.

The PharmCAT VCF preprocessing produces two types of **output**:

1. One or more PharmCAT-ready, single-sample VCF file(s)
2. A report of missing pharmacogenomics core allele defining positions in user's input


## How to run the PharmCAT VCF preprocessing tool

### Prerequisites

We assume that the input VCF files are prepared following the [Variant Call Format (VCF) Version >= 4.1](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

To run the tool, you need to download the following bioinformatic tools:
* [bcftools >= v1.16](http://www.htslib.org/download/)
* [htslib >= v1.16](http://www.htslib.org/download/)

We assume a working python3 installation with necessary dependencies:
* python >= 3.9.4
* pandas >= 1.4.2
* scikit-allel >= 1.3.5 (warning: testing compatibility with python >= 3.11)

To install necessary python packages, run the following code
```console
$ pip3 install -r requirements.txt
```

### Command line

To normalize and prepare a VCF file (single or multiple samples) for PharmCAT, run the following code substituted with proper arguments/inputs:

```console
$ python3 pharmcat_vcf_preprocessor.py -vcf path/to/file.vcf(.bgz)
```

**Mandatory** argument: `-vcf`.

-vcf
: Path to a single VCF file or a file containing the list of VCF file paths (one per line), sorted by chromosome position. All VCF files must have the same set of samples.  Use this when data for a sample has been split among multiple files (e.g. VCF files from large cohorts, such as UK Biobank).

  Example valid list file:
  ```
  chr1_set1.vcf
  chr1_set2.vcf
  chr2_set1.vcf
  chr2_set2.vcf
  ...
  ```
  Example invalid list file:
  ```
  chr3_set2.vcf
  chr2_set2.vcf
  chr1_set1.vcf
  chr1_set2.vcf
  ...
  ```

VCF files can have more than 1 sample and should be [bgzip](http://www.htslib.org/doc/bgzip.html) compressed. If not bgzip compressed, they will be automatically bgzipped.


**Optional** arguments:

-refVcf `<vcf_file>` <span class="altArg"><br />or --reference-pgx-vcf `<vcf_file>`</span>
: A sorted, compressed VCF of PGx core allele defining positions used by PharmCAT.  By default, the preprocessor will
look for `pharmcat_positions.vcf.bgz` under the current working directory.  You can find this VCF in the
`pharmcat_preprocessor-<release_version>.tar.gz` available from the PharmCAT GitHub releases page.

-refFna `<fna_file>` <span class="altArg"><br />or --reference-genome `<fna_file>`</span>
: The [GRCh38.p13](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/) FASTA file. The FASTA file has to be
decompressed and indexed (.fai). These mandatory files will be automatically downloaded (~0.9 GB) to the same directory
as the reference PGx VCF file (`-refVcf`) if not provided by user (see [Notes](#notes) for details).

-S `<txt_file>` <span class="altArg"><br />or --sample-file `<txt_file>`</span>
: The list of samples to be processed and prepared for PharmCAT. The file should contain one sample per line.

-bcftools `</path/to/bcftools>` <span class="altArg"><br />or --path-to-bcftools `</path/to/bcftools>`</span>
: bcftools must be installed. This argument is optional if bcftools is available in your PATH. 
If not, you can download and compile [bcftools](http://www.htslib.org/download/) and provide the path to the bcftools
program.

-bgzip `</path/to/bgzip>` <span class="altArg"><br />or --path-to-bgzip `</path/to/bgzip>`</span>
: bgzip must be installed. This argument is optional if bgzip is available in your PATH.
If not, bgzip is a part of the [htslib](http://www.htslib.org/download/). You can download and compile it and provide
the path to the bgzip program.

-o `<dir>` <span class="altArg"><br />or --output-dir `<dir>`</span>
: Directory to save preprocessed VCF to.  Default is the parent directory of the input VCF.

-bf `<name>` <span class="altArg"><br />or --base-filename `<name>`</span>
: Prefix of the output VCF files. Default is sample IDs from the input VCF(s).

-k <span class="altArg"><br />or --keep-intermediate-files</span>
: This option will help you save useful intermediate files, for example, a normalized, multiallelic VCF named `<base_input_file_name>.pgx_regions.normalized.multiallelic.vcf.bgz`, which will include all PGx regions from the first position to the last one in each chromosome as listed in the reference PGx VCF.

-0 <span class="altArg"><br />or --missing-to-ref</span>
: This option will add missing PGx positions to the output. Missing PGx positions are those whose genotypes are all missing "./." in every single sample.
  * This option will not convert "./." to "0/0" if any other sample has non-missing genotype at this position as these missing calls are likely missing for good reasons.
  * This **SHOULD ONLY BE USED** if you are sure your data is reference at the missing positions
    instead of unreadable/uncallable at those positions. Running PharmCAT with positions as missing vs reference can lead to different results.

-c <span class="altArg"><br />or --concurrent-mode</span>
: Enable concurrent mode.  This defaults to using one less than the number of CPU cores available.
Note that this is only useful if processing many files/samples.  With only a few files/samples, the overhead of
using concurrent mode is more than the benefit it may provide.

-cp `<num processes>` <span class="altArg"><br />or --max-concurrent-processes `<num processes>`</span>
: The maximum number of processes to use if concurrent mode is enabled.

### Output

All preprocessor output files will use the base filename of the input file unless otherwise specified using the `-bf`/`--base-filename` argument.  For example, if the input file is "study.vcf", then the base filename is "study".  If the input file is "biobank_files.txt" then the base filename is "biobank_files".

The preprocessor will produce one PharmCAT-ready VCF file per sample.  If there is only one sample, the output file is named `<base_filename>.preprocessed.vcf`.  If there are more than one samples, the output files are named `<base_filename>.<sample_id>.preprocessed.vcf`

If there are missing PGx positions, it will also produce a report named `<base_filename>.missing_pgx_var.vcf`.  This file only reports positions that are missing in _all_ samples.  If `-0`/`--missing-to-ref` is turned on, you can use this report to trace positions whose genotypes are missing in all samples (`./.`) in the original input but have now been added into the output VCF(s) as reference (`0/0`).


## Tutorial

### Case 1 - single-sample VCF
Imagine we have a VCF named *"test_1.vcf.bgz"* to be used in PharmCAT.
```console
$ gunzip -c test_1.vcf.bgz
$ cat test_1.vcf
<...header truncated...>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample_1
2	233760233	rs3064744	C	CAT	.	PASS	.	GT	1/0
2	233760233	rs3064744	CAT	C	.	PASS	.	GT	0/0
2	233760233	rs3064744	C	CATAT	.	PASS	.	GT	0/1
7	117548628	.	GTTTTTTTA	GTTTTTA	.	PASS	.	GT	0/1
```

Command to run the PharmCAT VCF preprocessor:
```console
$ python3 pharmcat_vcf_preprocessor.py -vcf test_1.vcf.bgz
```

VCF preprocessor will return two files in this test case.
1. one named *"test_1.preprocessed.vcf"*, which is a PharmCAT-ready VCF
2. the other named *"test_1.missing_pgx_var.vcf"* as a report of missing PGx positions.

Note that the chr7 variant is not used in PharmCAT and was removed by the PharmCAT VCF preprocessor.

```console
$ cat test_1.preprocessed.vcf
<...header truncated...>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample_1
chr2	233760233	rs3064744	CAT	C,CATATAT,CATAT	.	PASS	PX=UGT1A1	3/2

$ cat test_1.missing_pgx_var.vcf
<...header truncated...>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	PharmCAT
chr1	97078987	rs114096998	G	T	.	PASS	PX=DPYD	GT	0/0
chr1	97078993	rs148799944	C	G	.	PASS	PX=DPYD	GT	0/0
chr1	97079005	rs140114515	C	T	.	PASS	PX=DPYD	GT	0/0
<...truncated...>
```

### Case 2 - multi-sample VCF
Imagine we have a VCF named *"test_2.vcf.bgz"* that has two samples.
```console
$ gunzip -c test_2.vcf.bgz
<...header truncated...>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample_1	Sample_2
1	97740414	rs72549309	AATGA	A	.	PASS	.	GT	1/0	0/1
2	233760233	rs3064744	C	CAT	.	PASS	.	GT	1/0	0/0
2	233760233	rs3064744	CAT	C	.	PASS	.	GT	0/0	0/1
2	233760233	rs3064744	C	CATAT	.	PASS	.	GT	0/1	1/0
7	117548628	.	GTTTTTTTA	GTTTTTA	.	PASS	.	GT	0/1	1/0
10	94942212	rs1304490498	AAGAAATGGAA	A	.	PASS	.	GT	1/0	0/1
13	48037826	rs777311140	G	GCGGG	.	PASS	.	GT	1/0	0/1
19	38499645	rs121918596	GGAG	G	.	PASS	.	GT	1/0	0/1
22	42130727	.	AG	A	.	PASS	.	GT	1/0	0/1
M	1555	.	G	A	PASS	.	GT	1/0	0/1
```

Command to run the PharmCAT VCF preprocessor:
```console
$ python3 pharmcat_vcf_preprocessor.py -vcf test_2.vcf.bgz
```

VCF preprocessor will return three (3) files in this test case:
1. *"test_2.Sample_1.preprocessed.vcf"* 
2. *"test_2.Sample_2.preprocessed.vcf"*
3. *"test_2.missing_pgx_var.vcf"*

Note that the PharmCAT-ready VCFs will use the sample names from the input VCF.

```console
$ cat test_2.Sample_1.preprocessed.vcf
<...header truncated...>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample_1
chr1    97740410        rs72549309      GATGA   G       .       PASS    PX=DPYD       GT      1/0
chr2    233760233       rs3064744       CAT     C,CATATAT,CATAT .       PASS    PX=UGT1A1 GT      3/2
chr10   94942205        rs1304490498    CAATGGAAAGA     C       .       PASS    PX=CYP2C9     GT      1/0
chr13   48037825        rs777311140     C       CGCGG   .       PASS    PX=NUDT15     GT      1/0
chr19   38499644        rs121918596     TGGA    T       .       PASS    PX=RYR1       GT      1/0

$ cat test_2.Sample_2.preprocessed.vcf
<...header truncated...>
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample_2
chr1    97740410        rs72549309      GATGA   G       .       PASS    PX=DPYD       GT      0/1
chr2    233760233       rs3064744       CAT     C,CATATAT,CATAT .       PASS    PX=UGT1A1 GT      2/1
chr10   94942205        rs1304490498    CAATGGAAAGA     C       .       PASS    PX=CYP2C9     GT      0/1
chr13   48037825        rs777311140     C       CGCGG   .       PASS    PX=NUDT15     GT      0/1
chr19   38499644        rs121918596     TGGA    T       .       PASS    PX=RYR1       GT      0/1


$ cat test_2.missing_pgx_var.vcf
<...header truncated...>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	PharmCAT
chr1	97078987	rs114096998	G	T	.	PASS	PX=DPYD	GT	0/0
chr1	97078993	rs148799944	C	G	.	PASS	PX=DPYD	GT	0/0
chr1	97079005	rs140114515	C	T	.	PASS	PX=DPYD	GT	0/0
<...truncated...>
```

## Notes

PharmCAT uses [**GRCh38.p13**](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39).  It is available through the [NCBI RefSeq FTP site](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz).

PharmCAT takes this file and prepares it for use with the following commands:

```console
# curl -#fSL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz -o genomic.fna.gz
# gunzip genomic.fna.gz
# awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' genomic.fna | grep -v "^>chr\S*_" - | tr "\t" "\n" > genomic.short.fna
# bgzip -c genomic.short.fna > reference.fna.bgz
# samtools faidx reference.fna.bgz
# tar -czvf GRCh38_reference_fasta.tar reference.fna.bgz reference.fna.bgz.fai reference.fna.bgz.gzi
```

PharmCAT makes this indexed FASTA files available on [Zenodo](https://zenodo.org/record/7288118).
