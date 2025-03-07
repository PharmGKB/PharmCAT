---
parent: Using PharmCAT
title: Allele Frequency Analysis
permalink: using/Allele-Frequency-Analysis/
nav_order: 11
---
# Allele Frequency Analysis

As of PharmCAT 3.0, PharmCAT has much better built-in support for performing allele frequency analysis.


## Preparation

To do so, you will first need to generate [calls-only TSV reports](/using/Running-PharmCAT/#calls-only-tsv-reports)
using the `-reporterCallsOnlyTsv` flag.

You can also [add arbitrary sample metadata]((/using/Running-PharmCAT/#sample-metadata) (e.g. biogeographic data) to the
calls-only TSV reports using the `-sm`.


## Analysis

Once you've generated the calls-only TSV reports, you can generate allele frequencies with:

```console
# java -cp pharmcat.jar org.pharmgkb.pharmcat.stats.CalcAlleleFrequencies -i <path_to_report_tsv>  -o <directory>
```  

If you are using Docker or have installed the PharmCAT Pipeline, you can save a bit of typing and run:

```console
# calc_allele_freqs -i <path_to_report_tsv> -o <directory>
```  

Where:

-i `<path_to_report_tsv>`
: Input .report.tsv files.  The tool will look through all subdirectories to find all *.report.tsv files. 

-o `<directory>`
: Directory to output results to.  The tool will generate 1 Excel file per gene.
If not specified, it will use the same directory as the input. 


### Implementation Notes

If you are planning on doing your own frequency analysis, please note that the _Source Diplotype_ column can have
multiple entries!  If the `Named Allele Matcher` is unable to call a single diplotype (either because the data is not
phased or the data is incomplete), all the potential diplotypes will be concatenated together with an " OR ".

For example:
```tsv
CYP2B6	*4/[*6 + *10] OR *6/[*4 + *10]								no	9 / 9			
```

For frequency analysis, this should effectively be treated as unknown.

If this happens with genes that use the two lowest function variants to calculate phenotype/activity score,
such as [DPYD](/methods/Gene-Definition-Exceptions/#dpyd) and [RYR1](/methods/Gene-Definition-Exceptions/#ryr1),
the calls will be separated with an " AND ".

For example:
```tsv
DPYD	c.775A>G AND c.1627A>G (*5)									no		c.775A>G/c.1627A>G (*5)	Normal Metabolizer	2.0
```

For frequency analysis, _all_ alleles should be counted.


### Frequencies Based on Arbitrary Groupings

If the report contains additional sample metadata, you can also generate stats based on this metadata.

To do so, you will need to know the column number of the data (starting from 1).

If you only added a single column (e.g. _Biogeographic Group_), it would have been added as column 17 in the
.report.tsv file.  To add frequencies based on this column, use the `-pc` parameter:

```console
# calc_allele_freqs -i <path_to_report_tsv> -o <directory> -pc 17
```  

This will generate one Excel file per gene that contains diplotype and allele frequencies.
