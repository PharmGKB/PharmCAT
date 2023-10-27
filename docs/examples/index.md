---
title: Examples
permalink: examples/
nav_order: 5
---

# Example 1 - All Reference

This is a collection of example files that demonstrate input and output files from a single run of PharmCAT. For more 
information about running PharmCAT, read the [documentation](/using/Running-PharmCAT).


## Input

First, a properly formatted [single-sample VCF file](/using/VCF-Requirements/) with all positions specified.

- [example input VCF](pharmcat.example.vcf)

This is an example of an optional file of [outside diplotype calls](/using/Outside-Call-Format/). This one specifies a CYP2D6 call.

- [example outside calls](pharmcat.example.outsideCall.tsv)


## Output

### Named Allele Matcher output

The `Named Allele Matcher` component generates both HTML and JSON files with detailed information about how data in the 
sample VCF matches up with haplotype definitions.

- [example matcher HTML](pharmcat.example.match.html)
- [example matcher JSON](pharmcat.example.match.json)

### Phenotyper output

The `Phenotyper` component takes data from the `Named Allele Matcher` and combines it with outside call data to assign 
function and metabolizer values.

- [example phenotype JSON](pharmcat.example.phenotype.json)

### Reporter output

The `Reporter` component takes data from the `Phenotyper` and matches phenotypes to information found in CPIC guideline
data. This data is visible in an HTML report and also in a JSON file for machine parsing.

- [example report HTML](pharmcat.example.report.html)
- [example report JSON](pharmcat.example.report.json)


# Example 2 - Some Non-Reference Alleles

This is a collection of example files with CYP2C19 *2/*2 and CYP2B6 *1/*6 diplotypes. One input VCF file that contains 
the source variant data and the corresponding output files that would result from running that input file.

## Input VCF
 
- [example2 input VCF](pharmcat.example2.vcf)

## Output

- [example2 matcher HTML](pharmcat.example2.match.html)
- [example2 matcher JSON](pharmcat.example2.match.json)
- [example2 phenotype JSON](pharmcat.example2.phenotype.json)
- [example2 report HTML](pharmcat.example2.report.html)
- [example2 report JSON](pharmcat.example2.report.json)
