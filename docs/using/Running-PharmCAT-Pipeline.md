---
parent: Using PharmCAT
title: Running the PharmCAT Pipeline
permalink: using/Running-PharmCAT-Pipeline/
nav_order: 5
render_with_liquid: false
---
# Running the PharmCAT Pipeline

For convenience, we have a single script (`pharmcat_pipeline`) that simplifies the process of running the entire PharmCAT pipeline (the [VCF Preprocessor](/using/VCF-Preprocessor) and the [core PharmCAT tool](/using/Running-PharmCAT)).

This script tries to keep things as simple as possible.  If you need advanced features, you will need to run the individual parts of the pipeline directly.


## Prerequisites

You must have both the VCF preprocessor and the core PharmCAT tool installed correctly.  All the dependencies (python3, java, bcfools, bgzip) must already be in your PATH.  If necessary, you can customize which versions of these dependencies to use using environment variables:

* `JAVA_HOME` - the directory where the version of Java you want to use is installed
* `BCFTOOLS_PATH` - the full path to the `bcftools` program
* `BGZIP_PATH` - the full path to the `bgzip` program

The `pharmcat_pipeline` script is in the PharmCAT VCF Preprocessor tar file that's available on our [releases page](https://github.com/PharmGKB/PharmCAT/releases/).


## Usage

Standard use case:

```console
# pharmcat_pipeline <vcf_file>
```


### Details 

```
usage: pharmcat_pipeline [-s <samples> | -S <txt_file>]
                         [-0]
                         [-matcher] [-ma] [-matcherHtml] [-research <type>]
                         [-phenotyper]
                         [-reporter] [-rs <sources>] [-re] [-reporterJson]
                         [-o <dir>] [-bf <name>] [-del]
                         [-cp <num processes>]
                         [-v] [-V]
                         input file or directory

options:
  -h, --help            Show this help message and exit
  -v, --verbose         Print more verbose messages
  -V, --version         Show program's version number and exit

Input arguments:
  input file or directory
                        Path to a VCF file or a file of paths to VCF files (one file per line), sorted by
                        chromosome position.
  -s <sample id>, --samples <samples>
                        A comma-separated list of sample.
  -S <txt_file>, --sample-file <txt_file>
                        A file containing a list of samples, one sample per line.

Preprocessor arguments:
  -0, --missing-to-ref  Assume genotypes at missing PGx sites are 0/0. DANGEROUS!.

Named allele matcher arguments:
  -matcher              Run named allele matcher independently.
  -ma, --matcher-all-results
                        Return all possible diplotypes, not just top hits.
  -matcherHtml, --matcher-save-html
                        Save named allele matcher results as HTML.
  -research <type>, --research-mode <type>
                        Comma-separated list of research features to enable: [cyp2d6, combinations]

Phenotyper arguments:
  -phenotyper           Run phenotyper independently.

Reporter arguments:
  -reporter             Run reporter independently.
  -rs <sources>, --reporter-sources <sources>
                        Comma-separated list of sources to limit report to: [CPIC, DPWG]
  -re, --reporter-extended
                        Output extended report.
  -reporterJson, --reporter-save-json
                        Save reporter results as JSON.

Output arguments:
  -o <dir>, --output-dir <dir>
                        Directory for outputs. Defaults to the directory of the input VCF.
  -bf <name>, --base-filename <name>
                        Prefix for output files. Defaults to the same base name as the input.
  -del, --delete-intermediate-pharmcat-files
                        Delete intermediate PharmCAT files (saved by default).

Concurrency/Memory arguments:
  -cp <num processes>, --max-concurrent-processes <num processes>
                        The maximum number of processes to use when concurrent mode is enabled.
  -cm <size>, --max-memory <size>
                        The maximum memory PharmCAT should use (e.g. "64G"). This is passed on to Java
                        using the -Xmx flag.  Alternatively, set using the JAVA_MAX_HEAP environment
                        variable.
```

#### Inputs

* A VCF file
    * can be single or multisample
    * can be compressed (with gzip/bgzip)
* A text file containing a list of VCF file paths.  Files should be listed one per line, sorted by chromosome position.  Use this when data has been split among multiple files (e.g. VCF files from large cohorts, such as UK Biobank).
* A directory.  The script will look for VCF files to process, and they will be treated individually.

If the provided VCF file contains multiple samples, you can limit which samples get processed with either the `-s` or `-S` flag.


#### Concurrent Processing

The script will automatically attempt to use concurrent processing if possible.


#### Outside Calls

If you need to provide [outside calls](/using/Outside-Call-Format), you can do so with the following file naming 
convention: `<sample_id>.outside.tsv`.  For example, if your sample is `Sample_1`, then use `Sample_1.outside.tsv`.

If you have a single sample VCF file, you can also just use the basename of the VCF file.  For example, use
`mydata.outside.vcf` if you have a single sample VCF file called `mydata.vcf`.

If you have a multisample VCF file, (e.g. `multisample.vcf`) and have outside calls for `Sample_1`, then can use
`multisample.Sample_1.outside.tsv` instead.

These files need to be in the same directory as the VCF file.


#### Naming Conventions

PharmCAT uses the following sub-extensions in filenames to indicate which part of the pipeline it comes from:

* `.preprocessed` 
* `.match`
* `.phenotype`
* `.report`

These sub-extensions will be stripped off any filename to derive a base name for the file.

For example: the basename for `mydata.preprocessed.vcf` is `mydata`. 

In addition, outside call files with a `.outside` sub-extension will have the sub-extension stripped to derive the base
name for the file.
