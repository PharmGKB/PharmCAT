---
parent: Using PharmCAT
title: PharmCAT Pipeline
permalink: using/Running-PharmCAT-Pipeline/
nav_order: 2
render_with_liquid: false
---
# Running the PharmCAT Pipeline

The PharmCAT pipeline is composed of two components: the [VCF Preprocessor](/using/VCF-Preprocessor) and the 
[PharmCAT tool](/using/Running-PharmCAT).

For convenience, we have a single script (`pharmcat_pipeline`) that simplifies the process of running the pipeline. 
This script tries to keep things as simple as possible. If you need advanced features, you may need to run the
individual components of the pipeline directly.


## Prerequisites

This assumes that you are either [using Docker](/using/PharmCAT-in-Docker) or have already
[set up PharmCAT](/using/Setup-PharmCAT).


## Usage

Standard use case:

```console
# pharmcat_pipeline <vcf_file>
```


Full list of options: 

```
usage: pharmcat_pipeline [-s <samples> | -S <txt_file>] [-sm <tsv_file> ]
                         [-0] [--absent-to-ref] [-unspecified-to-ref] [-G] 
                         [-R <bed_file>]
                         [-matcher] [-ma] [-matcherHtml] [-research <type>]
                         [-phenotyper]
                         [-reporter] [-rs <sources>] [-re]
                         [-reporterHtml] [-reporterJson] [-reporterCallsOnlyTsv]
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
                        A comma-separated list of sample IDs.
                        Only applicable if you have multiple samples and only want to work on specific ones. 
  -S <txt_file>, --sample-file <txt_file>
                        A file containing a list of samples, one sample per line.
                        Only applicable if you have multiple samples and only want to work on specific ones.
  -sm <tsv_file>, --sample-metadata <tsv_file>
                        A TSV file containing sample metadata.

Preprocessor arguments:
  -0, --missing-to-ref               
                        Assume genotypes at absent or unspecified PGx sites are "0/0". DANGEROUS!
                        Modifying the data in this way can lead to different results in PharmCAT.

                        This option is equivalent to using both --absent-to-ref and --unspecified-to-ref.
                        Please consult the documentation of those flags for details.
  --absent-to-ref             
                        Assume genotypes at absent PGx sites are "0/0".  DANGEROUS!
                        Modifying the data in this way can lead to different results in PharmCAT.

                        This option will add PGx positions that are absent from the input VCF into the output as
                        homozygous reference ("0/0"). This SHOULD ONLY BE USED if you are sure your data is
                        reference at the absent positions instead of being unreadable/uncallable at the those positions.
  --unspecified-to-ref           
                        Assume unspecified genotypes ("./.") as "0/0" when every sample is "./.". DANGEROUS!
                        Modifying the data in this way can lead to different results in PharmCAT.
    
                        This option will convert an unspecified PGx position to homozygous reference ("0/0").
                        Unspecified PGx positions are those whose genotypes are unspecified ("./.") in every single sample.
                        As such, this check really only makes sense when working with multi-sample VCF files.
                        This option will not convert "./." to "0/0" when there is a specified genotype at a PGx position
                        as these "./." calls are likely left unspecified for good reasons.
  -G, --no-gvcf-check                
                        Bypass the gVCF check for the input VCF.
  -R `<bed_file>`, --retain-specific-regions `<bed_file>`      
                        A sorted .bed file indicating regions to retain in VCF.
                        For research use only. Additional variants are not used by PharmCAT and will slow PharmCAT down.

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
                        Comma-separated list of sources to limit recommendations to: [CPIC, DPWG, FDA]
  -re, --reporter-extended
                        Write an extended report (includes all possible genes and drugs, even if no data is available)
  -reporterHtml, --reporter-save-html
                        Save reporter results as HTML.  This is the default if no format is specified.
                        If any format is specified, only the specified formats will be saved.
  -reporterJson, --reporter-save-json
                        Save reporter results as JSON.
  -reporterCallsOnlyTsv, --reporter-save-calls-only-tsv
                        Save call results only as TSV.

Output arguments:
  -o <dir>, --output-dir <dir>
                        Directory for outputs. Defaults to the directory of the input VCF.
  -bf <name>, --base-filename <name>
                        Prefix for output files. Defaults to the same base name as the input.
  -del, --delete-intermediate-pharmcat-files
                        Delete intermediate PharmCAT files. Defaults to saving all files.

Concurrency/Memory arguments:
  -cp <num processes>, --max-concurrent-processes <num processes>
                        The maximum number of processes to use when concurrent mode is enabled.
  -cm <size>, --max-memory <size>
                        The maximum memory PharmCAT should use (e.g. "64G"). This is passed on to Java
                        using the -Xmx flag.  Alternatively, set using the JAVA_MAX_HEAP environment
                        variable.
```

### Inputs

#### 1. VCF data 

You must specify exactly one of the following:

* A VCF file
    * can be single or multi-sample
    * can be compressed (with gzip/bgzip)
* A directory containing VCF files
* For multi-file VCFs, a text file containing the list of VCF files.
  Files should be listed one per line, sorted by chromosome position.

{: .info }
Some large datasets come in multi-file VCFs, where a single extremely large VCF file is split across multiple smaller
VCF files. Examples of this are data from the International Genome Sample Resource (IGSR) 
(e.g. [Sample NA12878](https://www.internationalgenome.org/data-portal/sample/NA12878)) and UK Biobank
(e.g. [Data-Field 23156](https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=23156)).
<br /><br />
To feed multi-file VCFs to `pharmcat_pipeline`, you will need to create a text file that lists all the VCF files.
Files should be listed one per line, sorted by chromosome position.


#### 2. Outside call file(s)

If you need to provide [outside calls](/using/Outside-Call-Format), you can do so with the following file naming 
convention: `<sample_id>.outside.tsv`.  For example, if your sample is `Sample_1`, then use `Sample_1.outside.tsv`.

If you have a single sample VCF file, you can also use the base name of the VCF file.
For example, use `mydata.outside.vcf` if you have a single sample VCF file called `mydata.vcf`.

If you have a multi-sample VCF file, (e.g. `multisample.vcf`) and have outside calls for `Sample_1`, you can use
`multisample.Sample_1.outside.tsv` instead.

These files need to be in the same directory as the VCF file.


#### Filtering Samples

If the provided VCF file contains multiple samples, you can limit which samples get processed with either the
`-s` (a comma separated list of sample IDs) or `-S` flag (a file containing a list of sample IDs, one per line)


#### Concurrent Processing

The script will automatically attempt to use concurrent processing if possible.




#### Naming Conventions

PharmCAT uses the following sub-extensions in filenames to indicate which part of the pipeline it comes from:

* `.preprocessed` 
* `.match`
* `.phenotype`
* `.report`

These sub-extensions will be stripped off any filename to derive a base name for the file.

For example, the base name for `mydata.preprocessed.vcf` is `mydata`. 

In addition, outside call files with a `.outside.tsv` extension will have the extension stripped to derive the base name
for the file.


### Examples

#### Single-sample VCF

You have a [VCF file that contains data for a single person](/examples/pharmcat.example.vcf) saved as `/data/test.vcf`.

Base command, assuming you are in the directory in which you've installed PharmCAT:

```console
# pharmcat_pipeline /data/test.vcf
```

To run via Docker:

```console
# docker run --rm -v /data:/pharmcat/data pgkb/pharmcat pharmcat_pipeline data/test.vcf
```

Explanation: you are mounting `/data` as `/pharmcat/data`.  When you run the `pharmcat_pipeline` command, you are 
running it from the `/pharmcat` directory, so you can use a relative path to your VCF file.

If running Docker in [interactive mode](/using/PharmCAT-in-Docker#interactive-mode):

```console
# docker run --rm -v /data:/pharmcat/data -it pgkb/pharmcat
/pharmcat > pharmcat_pipeline data/test.vcf
```

{: .info }
For the remainder of these examples, I will only be providing the base command.
Hopefully, the examples above are enough to show you how to use the base command when running with Docker.  


#### Single sample VCF + outside call file

* You have a [VCF file that contains data for a single person](/examples/pharmcat.example.vcf) saved as
`/data/test.vcf`, under the sample ID `Sample_1` 
* You have an [outside call file](/examples/pharmcat.example.outsideCall.tsv) saved as `/data/test.outside.tsv`

```console
# pharmcat_pipeline /data/test.vcf
```

Notice that the command is the same as the previous example.

Since the VCF file only has a single sample, `pharmcat_pipeline` will notice that you have an outside call
file with the same base name and automatically use it.

The outside call file could also have been saved as `/data/Sample_1.outside.tsv` and `pharmcat_pipeline` would have
noticed that you have an outside call file with the name of the sample and automatically use it.


#### Multi-sample VCF file (+ outside call files)

You have a [VCF file that contains data for multiple people](/examples/multisample.vcf) saved as
`/data/multisample.vcf`.

```console
# pharmcat_pipeline /data/multisample.vcf
```

Notice that except for the filename, the command is the same as the previous example.

`pharmcat_pipeline` will automatically handle each sample it finds in the VCF file.

If you want to provide outside calls, you _must_ use the sample ID file naming scheme.
You do not need to provide outside call files for every sample if you do not wish to.
For example, `split_vcf_list.txt` has 2 samples (Sample_1 and Sample_2).
You can just have a `/data/Sample_1.outside.tsv` and no `/data/Sample_2.outside.tsv`, or vice versa.


#### Using multi-file VCF data

Some large datasets come in multi-file VCFs, where a single extremely large VCF file is split across multiple smaller
VCF files. For example, data from the International Genome Sample Resource (IGSR) and UK Biobank come in this format:

* [Data-Field 23156](https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=23156) from UK Biobank
* [Sample NA12878](https://www.internationalgenome.org/data-portal/sample/NA12878) from IGSR

To feed multi-file VCFs to `pharmcat_pipeline`, you will need to create a text file that lists all the VCF files (one 
per line), sorted by chromosome position.

You can download the sample [split VCF list](/examples/split_vcf_list.txt) to test this out.
You will also need the actual [split VCF files](/examples/split_vcf.tar).
Untar this file into the same directory of the split VCF list file.

```console
# pharmcat_pipeline split_vcf_list.txt
```

Notice that except for the filename, the command is the same as the first example.

If you want to provide outside calls, you _must_ use the sample ID file naming scheme.
You do not need to provide outside call files for every sample if you do not wish to.
For example, `split_vcf_list.txt` has 2 samples (Sample_1 and Sample_2).
You can just have a `/data/Sample_1.outside.tsv` and no `/data/Sample_2.outside.tsv`, or vice versa.


#### Using multiple VCF files

If you have multiple independent VCF files (single or multi-sample) and would like to process them in one swoop, you can
just point `pharmcat_pipeline` at the directory the files are in.

If you have `[/data/test1.vcf](/examples/pharmcat.example.vcf)`, `[/data/test2.vcf](/examples/pharmcat.example2.vcf)`
and `[/data/test1.outside.tsv](/examples/pharmcat.example.outsideCall.tsv)`, then you can run:

```console
# pharmcat_pipeline /data
```

Notice that except for the directory, the command is the same as the first example.

If you want to provide outside calls, you can either use the VCF file base name approach
(e.g. `/data/test1.outside.tsv`) if the VCF file only has a single sample. 


You do not need to provide outside call files for every sample if you do not wish to.
For example, `split_vcf_list.txt` has 2 samples (Sample_1 and Sample_2).
You can just have a `/data/Sample_1.outside.tsv` and no `/data/Sample_2.outside.tsv`, or vice versa.

