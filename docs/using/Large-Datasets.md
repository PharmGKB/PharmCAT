---
parent: Using PharmCAT
title: Working with Large Datasets
permalink: using/Large-Datasets/
nav_order: 10
has_toc: false
---
# Working with Large Datasets

PharmCAT has successfully been run on large datasets such as the
[UK Biobank data](https://www.ukbiobank.ac.uk/), the [Penn Medicine Biobank data](https://pmbb.med.upenn.edu/), and the
[All of Us data](https://researchallofus.org/).

Documenting how to do so is challenging because there are many ways to run PharmCAT, and the PharmCAT team doesn't
always have direct access to these large datasets. 

We are actively working on improving PharmCAT for use on large data sets and would love to
[get your feedback](mailto:pharmcat@clinpgx.org) on how it worked for you and how you'd like to see it improved.
We would also appreciate any notes you may have on either how you used PharmCAT on different systems and the costs you
encountered.


## Recommendations

While there are many ways to run PharmCAT on large datasets that are dependent on how you plan to work with the data,
the simplest way is to just use the [PharmCAT Pipeline](/using/Running-PharmCAT-Pipeline).

### Output Control

The most impactful control you have over PharmCAT's performance is to limit its output to avoid unnecessary file I/O.

At the top of the list is to use the `-del` flag.
This will stop PharmCAT from producing intermediate data files that are only necessary if you're going to dive into
PharmCAT internals.

If you are only interested in diplotype calls and/or frequency analysis, use the `-reporterCallsOnlyTsv` flag.  This
will generate output most relevant to you.  See [Allele Frequency Analysis](/using/Allele-Frequency-Analysis)

If you are interested in getting drug recommendations, you will also need to use the `-reporterHtml` and/or
`-reporterJson` flags depending on the format of the output you need.


### Preprocessing VCFs

While using the PharmCAT Pipeline is easy, consider running the [VCF Preprocessor](/using/VCF-Preprocessor) and
the [PharmCAT](/using/Running-PharmCAT) separately, especially if you have large VCF files.

Again, disk IO is usually the bottleneck, and the VCF preprocessor will extract only what PharmCAT needs from your large
VCF files.

This will also save you a lot of time if you need to re-run PharmCAT with different parameters. 


#### Advanced Flags

This section is a work in progress.  We would like to provide recommendations for concurrency (max processes) and memory
based on different compute platforms (e.g. Amazon AWS, Google Cloud, Microsoft Azure).

If you have experience running PharmCAT, or would be willing to help us test this out, please reach out and
[let us know](mailto:pharmcat@clinpgx.org).


## Past Experiences

To the best of our ability, we have gathered steps that have been used to work with large datasets on various platforms.
However, please note that this is really just a snapshot in time when the analysis was performed and may not be the
best way to do so now:

* [Multi-Sample Analysis](/using/Multi-Sample-Analysis) - steps taken to run PharmCAT on some platforms
* [PharmCAT Run Time and Cost](/using/Runtime-Cost) - a collection of reported cost for running PharmCAT on various
  large datasets
