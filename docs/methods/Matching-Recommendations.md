---
parent: Methods
title: Matching Recommendations
permalink: methods/Matching-Recommendations/
nav_order: 6
---
# Matching Recommendations

PharmCAT was built, among other tasks, to take genotyping results and match them to applicable recommendations published
by expert groups like CPIC and DPWG. This document explains how this process works within the `Phenotyper` and
`Reporter` modules of the PharmCAT tool.

{:.info}
> All guideline information discussed below is authored by the source expert groups and not
by PharmCAT. PharmCAT attempts to use the information as the guideline authors intended and to transform the data as
little as possible to make it work with the application.


## Data Sources 

Guideline recommendation data has a typical format. A guideline can apply to one or more drugs. If the guideline applies
to multiple drugs, it's typically because they belong to the same drug class or are used for the same therapeutic case.
The guideline will also specify one or more genes that should be genotyped to get specific recommendation
information. This will lead to possible drug-gene combinations in the guideline where some may have recommendations
while others may have "no recommendation" as indicated by the guideline authors. The specific recommendations are usually a 
collection of text annotations of different types (e.g. implications, recommendation, strength of recommendation) that
apply to specific gene results like a particular phenotype or the presence of a particular allele. A guideline may have
more than one group of annotations for the same genotype if the guideline can apply to multiple populations that have 
different recommendations. The method to get from a genotype (e.g. a diplotype like CYP2C9 `*1/*3`) to gene results
vary by guideline, so it is important to closely read and understand the guideline before interpreting genotypes in that
context.

Currently, PharmCAT pulls guideline data from two sources. The first is CPIC. CPIC makes their guideline
recommendations and supporting data available in machine-readable formats through their REST API and via a file archive
service. CPIC publishes drug, phenotype, and allele definition data that PharmCAT consumes to do its work. 

The second source is DPWG (the Dutch Pharmacogenomics Working Group). DPWG does not publish their guideline information
in a well-structured machine-readable format (they only publish as PDFs), so PharmCAT relies on the annotated DPWG
guidelines in PharmGKB. PharmGKB curators read the source PDFs and curate the guideline data, so it is available via
PharmGKB's API and file archive service. PharmGKB also provides allele definition data for genes that are used in DPWG
guidelines but are not in CPIC.

PharmCAT pulls the aforementioned guideline data and reorganizes it to prepare for reporting. The annotation data pulled
from the sources may apply to the same drug (e.g. CPIC and DPWG both have annotations for amitriptyline) so PharmCAT
tags annotations with a "source" property that describes where it originated and then groups those annotations under a
"drug" object so the information can be grouped by drug in the final output. PharmCAT also takes the different formats
of annotations and translates them into a single format that can be displayed or compared in the same context.


## Recommendation Lookup

The `Phenotyper` handles the "gene result" assignment, and the `Reporter` handles the recommendation assignments.


{:.info}
> Each guideline describes the method it uses to match recommendation text to genotyping results.
> 
> If you have questions about why a lookup method works the way it does, it will most probably be answered by reading
the original guideline publication.

Each gene in PharmCAT will use one of the methods described below to come up with a "gene result" that can be used to
match genotyping data to guideline recommendations.


### Method 1: Phenotype lookup

Currently, most guidelines use a common "phenotype" lookup method. Guidelines have different recommendations that apply
to individual gene phenotypes. These phenotypes are based on a 
[standard vocabulary](https://cpicpgx.org/resources/term-standardization/) with values like "Normal Metabolizer" or
"Poor Metabolizer". The guideline will also have a translation table that describes what combinations of gene function
are assigned to a given phenotype. For example, a diplotype of one "normal function" allele and one "no function" allele
could be assigned to a "Intermediate Metabolizer". The guideline source also includes mapping for named alleles to 
clinical functional status. For example, `*1` maps to "normal function" and `*7` maps to "no function".

Using these mappings, PharmCAT: 

1. takes a diplotype from the `Named Allele Matcher` or from an outside call,
2. assigns function values to each alleles,
3. matches that function combination to a gene phenotype,
4. and uses that phenotype as the gene result to match a recommendation in a guideline.

Some genes may use the algorithm of a phenotype lookup but instead of using gene function assignments, will use other
phenotypic descriptions. For example, CACNA1S has a "Normal Function" assignment for the `Reference` allele but uses
"Malignant hyperthermia associated" for the other alleles (e.g. `c.520C>T`) assigned in PharmCAT. This is not a function
but gets treated like a function in this algorithm.

See [Phenotypes List](/Phenotypes-List) for a complete list of phenotypes used by PharmCAT.


### Method 2: Activity Score lookup

The activity score lookup method is similar to the phenotype lookup method but also includes the assignment of an
activity score in addition to a phenotype when a gene result is assigned.  The difference is that the activity score is
used for looking up recommendations. 

Example 1: A CYP2C9 diplotype of `*1/*4` has one normal function allele (`*1`) which is assigned an activity value of
`1.0` and one decreased function allele (`*4`) which is assigned an activity value of `0.5`. Combined, those two 
activity values give an activity score of `1.5`, which maps to an "Intermediate Metabolizer" phenotype for CYP2C9.
We then use the `1.5` activity score as the gene result to look up the appropriate recommendations for CYP2C9.

Example 2: A CYP2C9 diplotype of `*1/*3` has one normal function allele (`*1`) with an activity value of `1.0` and one
no function allele (`*3`) with activity value of `0.0`. This yields a total activity score of `1.0` for this diplotype,
which is also an "Intermediate Metabolizer" for CYP2C9.  Even though `*1/*3` has the same phenotype as `*1/*4` from
Example 1, it may have a different recommendation since it is the score (`1.0`) that is used for lookup up 
recommendations and not the phenotype.

See [Phenotypes List](/Phenotypes-List) for a complete list of activity scores used by PharmCAT.


### Method 3: Single position lookup

Some genes do not have defined named alleles and, instead, rely on a single SNP for their gene result. There is no final
"phenotype" assigned like the previous two methods. Instead, the "function" is used as the gene result to lookup
recommendations.

For example, ABCG2 uses the `rs2231142` SNP for recommendation lookups with two possible alleles that can be used for
lookups: `reference (G)` and `variant (T)`. Each allele is assigned a function, but the diplotype itself is also
assigned an overall function which is then used for lookup.


### Method 4: Allele status lookup

While function and phenotype assignments are common for genes used by PharmCAT, there are some genes where function and
phenotype are not assigned. For these genes, the mere presence of particular alleles is what the gene result for lookup
is based on.

For example, HLA-B is used in a few guidelines and each guideline can use a different allele for matching
recommendations. For example, the CPIC guideline for abacavir looks for the presence of the `*57:01` allele but the
allopurinol guideline looks for the presence of the `*58:01` allele. This means one diplotype given to PharmCAT could
result in more than one gene result for HLA-B, one for each allele that is used in recommendation lookup.


## Multi-gene Guidelines

Some guidelines will use more than one gene to look up recommendations. In that case, each gene uses its own lookup
method as described above, and then the combination of the gene results is used to match a guideline recommendation. 

For example, the CPIC phenytoin guideline uses both CYP2C9 (an activity score gene) and HLA-B (an allele status gene) to
match to the recommendation text. If a sample is CYP2C9 `*1/*1` and the HLA-B `*15:02` allele is not present then a combined
gene result of a CYP2C9 activity score of `2.0` and an HLA-B result of `*15:02 negative` will be used to match to a
recommendation.


## Populations & Multiple Recommendations

A gene result calculated per the methods above can result in more than one recommendation for a particular guideline
from a particular source. This happens because the guideline was written with multiple populations in mind and each
recommendation that matches should specify a different population.

For example, the CPIC atomoxetine guideline has recommendations for both an adult population and a pediatric population.
The populations are not always based on age. Some populations may be based on how long the person has been taking a 
particular medication or a particular clinical finding, or something else.

{:.info}
> Read the source guidelines to find how the population should be interpreted.


## Exceptions for recommendation lookup

There are genes that need to use special, individual logic to determine their gene results for recommendation lookup.

See [Gene Definition Exceptions](/methods/Gene-Definition-Exceptions) for details on these genes.

In the `Reporter` JSON output, the diplotype used to look up the recommendation is specified in the
`recommendationDiplotype` field.

The "real" diplotype is stored in the `sourceDiplotype` field.
This will be either what the `Named Allele Matcher` called, or was provided as an outside call.
This is the value that is displayed in the PharmCAT reports.
