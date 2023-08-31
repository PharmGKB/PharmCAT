---
parent: Using PharmCAT
title: Running PharmCAT
permalink: using/Running-PharmCAT/
nav_order: 1
render_with_liquid: false
mermaid: true
---
# Running PharmCAT

This document teaches you how to use the core PharmCAT tool.


## Prerequisites

### Install the Software
You can skip this if are [running PharmCAT in Docker](/using/PharmCAT-in-Docker).

1. You will need [Java 17 or newer](https://adoptium.net/index.html?variant=openjdk17&jvmVariant=hotspot).
2. Download the PharmCAT Jar file from our [releases page](https://github.com/PharmGKB/PharmCAT/releases/).


### Prepare Your VCF Files

PharmCAT takes VCF files as input.

{: .warn}
> Please make sure you have read and understand PharmCAT's [VCF requirements](/using/VCF-Requirements).
>
> If you are not preparing your VCF files yourself, we highly recommend you run your VCF file through PharmCAT's [VCF Preprocessor](/using/VCF-Preprocessor).


## Running PharmCAT

This is the basic command, which will take your VCF file and produce the PharmCAT report:

```console
# java -jar <path_to_pharmcat_jar_file> -vcf <vcf_file>
```

Where:

-jar `<path_to_jar_file>`
: The compiled PharmCAT Jar file

-vcf `<vcf_file>`
: Input VCF file (must comply with PharmCAT's [VCF requirements](/using/VCF-Requirements))

By default, the output will be saved to the same directory as the input VCF file and will use the same base file name.  For example:

```console
# java -jar pharmcat.jar -vcf /tmp/sample1.vcf
Saving named allele matcher JSON results to /tmp/sample1.match.json
Saving phenotyper JSON results to /tmp/sample1.phenotype.json
Saving reporter HTML results to /tmp/sample1.report.html
```

Notice that all intermediary files will be kept.

To control this behavior, provide:


-o `<dir>` <span class="altArg"><br />or --output-dir `<dir>`</span>
: Directory to save files to

-bf `<name>` <span class="altArg"><br />or --base-filename `<name>`</span>
: The base name (without file extensions) used for output files

-del <span class="altArg"><br />or --delete-intermediary-files</span>
: Delete intermediary output files

Example:

```console
# java -jar pharmcat.jar -vcf /tmp/input.vcf -o /tmp/results -bf sample
Saving named allele matcher JSON results to /tmp/results/sample.match.json
Saving phenotyper JSON results to /tmp/results/sample.phenotype.json
Saving reporter HTML results to /tmp/results/sample.report.html

# java -jar pharmcat.jar -vcf /tmp/input.vcf --output-dir /tmp/results -del
Saving reporter HTML results to /tmp/results/input.report.html
```

If the provided VCF file contains multiple samples, you can limit which samples get processed with either:

-S `<txt_file>` <span class="altArg"><br />or --sample-file `<txt_file>`</span>
: The list of samples to be processed and prepared for PharmCAT. The file should contain one sample per line.

-s `<samples>` <span class="altArg"><br />or --samples `<samples>`</span>
: A comma-separated list of sample IDs.


### Outside Calls

If you need to provide diplotypes directly to PharmCAT, you can do so using an ["outside calls" file](/using/Outside-Call-Format).  You might want to do this for genes that PharmCAT does not call directly, or to override PharmCAT's call.

To do so, provide:

-po `<tsv_file>` <span class="altArg"><br />or --phenotyper-outside-call-file `<tsv_file>`</span>
: Path to an outside call file (TSV)

Example:

```console
# java -jar pharmcat.jar -vcf /tmp/sample.vcf -po /tmp/outside_calls.tsv
Saving named allele matcher JSON results to /tmp/results/sample.match.json
Saving phenotyper JSON results to /tmp/results/sample.phenotype.json
Saving reporter HTML results to /tmp/results/sample.report.html
```


## Advanced Usage

Remember that the PharmCAT tool is composed of 3 modules:  the `Named Allele Matcher`, the `Phenotyper`, and the `Reporter`.

```mermaid
graph LR;
  MI{{vcf file}} --> M[Named Allele Matcher];
  M --> P[Phenotyper];
  PI{{"outside call file"}} -. optional .-> P;
  P --> R[Reporter];
  R --> RO([HTML report]);
  class M,P,R module;
  class PI optional;
```

Each module has its own arguments to customize its behavior.

#### Named Allele Matcher

-matcher
: run Named Allele Matcher

-ma <span class="altArg"><br />or --matcher-all-results</span>
: return all possible diplotypes, not just top hits

-matcherHtml <span class="altArg"><br />or --matcher-save-html</span>
: save named allele matcher results as HTML

#### Phenotyper

-phenotyper
: run Phenotyper

-pi `<json_file>` <span class="altArg"><br />or --phenotyper-input `<json_file>`</span>
: JSON results from named allele matcher

-po  `<tsv_file>` <span class="altArg"><br />or --phenotyper-outside-call-file `<tsv_file>`</span>
: path to an outside call file (TSV)

#### Reporter

-reporter
: run Reporter

-ri `<json_file>` <span class="altArg"><br />or --reporter-input `<json_file>`</span>
: JSON results from phenotyper

-rt `<title>` <span class="altArg"><br />or --reporter-title `<title>`</span>
: text to add to the report title

-rs `<CPIC or DPWG>` <span class="altArg"><br />or --reporter-sources `<CPIC or DPWG>`</span>
: comma-separated list of sources to limit recommendations to (defaults to both)

-re <span class="altArg"><br />or --reporter-extended</span>
: generate extended report (includes all possible genes and drugs, even if no data is available)

-reporterJson
: save reporter results as JSON


### Running Individual Modules

#### Just the `Named Allele Matcher`

This will call the diplotypes for the input VCF file.

Examples:

```console
# java -jar pharmcat.jar -matcher -vcf /tmp/sample.vcf
Saving named allele matcher JSON results to /tmp/sample.match.json
```

#### Just the `Phenotyper`

This will take matcher and/or call data and output function and phenotype information for them.

Examples:

```console
# java -jar pharmcat.jar -phenotyper -pi /tmp/sample1.match.json
Saving phenotyper JSON results to /tmp/results/sample1.phenotype.json

# java -jar pharmcat.jar -phenotyper -pi /tmp/sample2.match.json -po /tmp/outside_calls.tsv
Saving phenotyper JSON results to /tmp/results/sample2.phenotype.json

# java -jar pharmcat.jar -phenotyper -po /tmp/outside_calls.tsv
Saving phenotyper JSON results to /tmp/results/outside_calls.phenotype.json
```

#### Just the `Reporter`

This will take the phenotyper data and output the relevant drug annotations in a comprehensive HTML report.

Examples:

```console
# java -jar pharmcat.jar -reporter -ri /tmp/sample1.phenotype.json
Saving reporter HTML results to /tmp/results/sample1.report.html

# java -jar pharmcat.jar -reporter -ri /tmp/sample2.phenotype.json -rt "Awesome Report" -reporterJson
Saving reporter JSON results to /tmp/results/sample1.report.json
Saving reporter HTML results to /tmp/results/sample1.report.html
```


#### Combinations

Run `Named Allele Matcher` and `Phenotyper`:

```console
# java -jar pharmcat.jar -matcher -vcf /tmp/sample1.vcf -phenotyper
Saving named allele matcher JSON results to /tmp/sample1.match.json
Saving phenotyper JSON results to /tmp/results/sample1.phenotype.json

# java -jar pharmcat.jar -matcher -vcf /tmp/sample2.vcf -phenotyper -po /tmp/outside_calls.tsv
Saving named allele matcher JSON results to /tmp/sample2.match.json
Saving phenotyper JSON results to /tmp/results/sample2.phenotype.json
```

Run `Phenotyper` and `Reporter`:

```console
# java -jar pharmcat.jar -phenotyper -pi /tmp/sample1.phenotype.json -po /tmp/outside_calls.tsv -reporter
Saving phenotyper JSON results to /tmp/results/sample1.phenotype.json
Saving reporter HTML results to /tmp/results/sample1.report.html

# java -jar pharmcat.jar -phenotyper -po /tmp/outside_calls.tsv -reporter
Saving phenotyper JSON results to /tmp/results/outside_calls.phenotype.json
Saving reporter HTML results to /tmp/results/outside_calls.report.html
```


### Research Mode Options

PharmCAT has functionality meant for **research use only**. You can read about these features on the [Research Mode](/using/Research-Mode) page.


### Custom Definition Files

Advanced users can provide PharmCAT with custom allele definitions:

-def `<dir>` <span class="altArg"><br />or --definitions-dir `<dir>`</span>
: directory containing named allele definitions (JSON files)

This can be used to get PharmCAT to call diplotypes for genes that PharmCAT does not support by default.  Unless these are genes that PharmCAT supports through [outside calls](/Genes-Drugs/#genes-handled-by-outside-callers), PharmCAT will **NOT** be able to match them to any recommendations.
