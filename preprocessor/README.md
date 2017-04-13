# PharmCAT Preprocessor

This Dockerfile can be used to prepare a VCF for input into PharmCAT.


## Introduction

:exclamation: :exclamation: :exclamation: __FOR TESTING ONLY__ :exclamation: :exclamation: :exclamation:

This application is __under development__ and has __not been officially released__. Watch this repo or check the [releases](../../releases) page for the official release

:warning: PhamCAT assumes no responsibility for any injury to person or damage to persons or property arising out of, or related to any use of PharmCAT, or for any errors or omissions. The user recognizes they are using PharmCAT at their own risk.


In order to run PharmCAT there are several requriements for the input VCF described in the [wiki](../../wiki). Briefly the variants must be aligned to grc38, must contain all positions (even reference or missing positions), and the variants must exactly match the CPIC defintions.


This docker script provides a template that you can adapt to your own data. It has mainly been used for in-house testing on public data sets like 1k genomes and Genome In a Bottle, where we only have access to the variant VCF files. Due to this there are several compromises made, such as presuming reference for all missing positions. Therefore this procedure is suitable for testing, but we heavily suggest additional steps if running in on inhouse data, such as using a VCF with all sites ('EMIT_ALL_SITES') flag in GATK, or a tool like [jvarkit](https://github.com/lindenb/jvarkit/wiki/FixVcfMissingGenotypes) to fill in the reference/missing positions.

The steps taken in the script are shown in the prepare.sh.


##  Running the script:

In order to run this script you will need to have docker [installed](https://docs.docker.com/engine/getstarted/step_one/).

Once you have docker installed change directories into the directory containing the Dockerfile and type the following at the command line:
```
docker build -t pharmcat .
```

You can now change directories into the workng directory or a fresh test directory. For instance if you would like some test data you can download the Genome in a Bottle grc38 vcf files:
```
    wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/platinum_genomes/2016-1.0/hg38/small_variants/NA12878/NA12878.vcf.gz
    wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/platinum_genomes/2016-1.0/hg38/small_variants/NA12878/NA12878.vcf.gz.tbi
```

And can then run it using the following command (if in the /work/inhouse. directory):
```
docker run -i -v /work/inhouse/:/input/ pharmcat /input/NA12878.vcf.gz
```

In general the script can be run using when the current working directory contains the file input file:
```
docker run -i -v <current path>:/input/ pharmcat /input/NA12878/<input_file>
```

Multi-sample and gzipped/ungzipped vcf files can be handled.  In a multi-sample file all samples will be processed.


## Output
In the above example the out put will be the following files:
   * <filename>.<sample>.final.vcf - a final VCF file ready for PharmCAT input.
   * <filename>.<sample>.report.html - an initial pharmcat run is also carried out. This the output report.
   * <filename>.<sample>.matcher.html - this is diplotype matching results for the PharmCAT run.
   * <filename>.<sample>.report.json - This is a json version of the report for easy parsing.

In addition the script log will produce a summary output.


## Contact
For questions about the PharmCAT project, contact [pharmcat@pharmgkb.org](mailto:pharmcat@pharmgkb.org)
