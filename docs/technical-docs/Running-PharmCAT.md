---
title: Running PharmCAT
parent: Using PharmCAT
permalink: technical-docs/Running-PharmCAT
nav_order: 2
---
# Running PharmCAT

This will run the entire PharmCAT tool.

## Requirements

You will need [Java 14 or newer](https://adoptium.net/index.html?variant=openjdk17&jvmVariant=hotspot).

You can either build a fresh copy of the Jar file yourself or download a pre-compiled Jar file from our [releases page](https://github.com/PharmGKB/PharmCAT/releases/).

PharmCAT takes VCF files as input.

_Note_: PharmCAT does not need a network connection during runtime. All data needed by PharmCAT is either supplied by 
the user on the command-line or self-contained.

:warning: Please make sure you have read and undertand PharmCAT's [VCF requirements](/specifications/VCF-Requirements).

## Building (optional)

Checkout the repo and from the base repo directory run:

```commandline
# ./gradlew shadowJar
```

This will build a "fat" jar with bundled dependencies in `build/libs`. You can use this jar file in the following section.

For more information on building PharmCAT, check [Building PharmCAT](https://github.com/PharmGKB/PharmCAT/wiki/Building-PharmCAT).


## Running

From the command line:

```commandline
# java -jar <path_to_pharmcat_jar_file> -vcf <sample_file> -o <output_dir>
```

Where:

* __-jar__ `<path_to_jar_file>` = __required__, the compiled PharmCAT Jar file
* __-vcf__ `<sample_file>` = __required__, sample VCF file (:warning: Please read [VCF requirements](/specifications/VCF-Requirements))
* __-o__ `<output_dir>` = __required__, diretory path to write result files to
* __-f__ `<output_name>` = _optional_, a base filename to use for output files (e.g. `<output_name>.html`)
* __-a__ `<outside_call_file>` = _optional_, [gene call TSV file](/specifications/Outside-Call-Format) from the an outside tool (like Astrolabe)
* __-k__ = _optional_, keep the interim output files from the NamedAlleleMatcher
* __-j__ = _optional_, flag to write reporter JSON data (will be `<output_name>.report.json`)
* __-pj__ = _optional_, flag to write phenotyper JSON data (will be `<output_name>.phenotyper.json`)


### Custom Data

PharmCAT includes the raw data it relies on.  However, you can change this by using the following arguments:

* __-na__ `<definitions_dir>` = _optional_, a directory containing allele definitions to use instead of the default packaged allele definitions

The latest version of the dosing guideline annotations can be [downloaded from PharmGKB](https://www.pharmgkb.org/downloads).


## Running Individual Components

PharmCAT has multiple components that are run internally. Each of these components can be run individually if the 
output of the other components is not needed. Below are examples of how to run particular components.

### Running the NamedAlleleMatcher

The NamedAlleleMatcher will match given sample VCF data to the named allele definitions in PharmCAT. This does not do 
anything with outside call data.

From the command line:

```commandline
> java -cp <path_to_pharmcat_jar_file> org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher -vcf <vcf_file> -html <html_file>
```

Where:

* __-vcf__ `<sample_file>` = __required__, sample VCF file
* __-html__ `<html_file>` = __optional__, the HTML file to write the results to
* __-json__ `<json_file>` = __optional__, the JSON file to write the results to


#### Custom Data

PharmCAT includes the raw data it relies on.  However, you can change this by using the following arguments:

* __-d__ `<definitions_dir>` = _optional_, a directory containing allele definitions to use instead of the default packaged allele definitions


### Running the Phenotyper

The Phenotyper will take match or call data and return function and phenotype information for them.

There are two ways to run this. First option, you can run directly with VCF sample data:

```commandline
> java -cp <path_to_pharmcat_jar_file> org.pharmgkb.pharmcat.phenotype.Phenotyper -vcf <vcf_file> -f <path_to_output_json>
```

Second option, you can run with output from the NamedAlleleMatcher:

```commandline
> java -cp <path_to_pharmcat_jar_file> org.pharmgkb.pharmcat.phenotype.Phenotyper -c <call_file> -f <path_to_output_json>
```

Where:

* __-vcf__ `<sample_file>` = sample VCF file
* __-c__ `<call_file>` = JSON call data output from the NamedAlleleMatcher
* __-f__ `<path_to_output_json>` = the path to an output JSON file
* __-o__ `<path_to_outside_call>` = __optional__, [a TSV of outside caller information](/specifications/Outside-Call-Format)


### Running the Reporter

The reporter takes allele matcher and call data and finds drug annotations that are relevant to them. The annotations
are compiled into JSON and HTML reports.

From the command line:

```commandline
> java -cp <path_to_pharmcat_jar_file> org.pharmgkb.pharmcat.reporter.Reporter -c <call_file> -a <outside_call_file> -o <output_html>
```

Where:

* __-p__ `<phenotyper-file>` = __required__, JSON result file from the `Phenotyper`
* __-o__ `<output-file>` = __required__, file path to write HTML result to
* __-j__ `<output-json>` = _optional_, file path to write JSOn data to
* __-t__ `<title>` = _optional_, text to ad to the report title


#### Custom Data

PharmCAT includes the raw data it relies on.  However, you can change this by using the following arguments:

* __-g__ `<guidelines_dir>` = _optional_, directory containing JSON files of dosing guidelines instead of the default packaged guidelines

The latest version of the dosing guideline annotations can be [downloaded from PharmGKB](https://www.pharmgkb.org/downloads).
