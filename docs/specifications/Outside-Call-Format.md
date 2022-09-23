---
parent: Specifications
title: Outside Call Format
permalink: specifications/Outside-Call-Format/
nav_order: 4
---
# Outside Call Format

Typically, PharmCAT uses variant call data to match diplotypes used to find annotations. However, you can also give
diplotypes, phenotypes, or other allele calls to PharmCAT that were called by other tools.

These **outside call files** can be supplied to the [PharmCAT](/using/Running-PharmCAT#outside-calls) tool using the
`-po` flag.

Calls specified in this file will override results from the `Named Allele Matcher`.


## File format

The **outside call file** format is a tab-separated file. Lines starting with `#` will be ignored.

Each line has up to 4 fields, separated by tabs:

1. HGNC gene symbol (_required_)
2. Diplotype or single allele call (_required if third and fourth columns not specified_)
3. Phenotype or other gene result (_required if second and fourth columns not specified_)
4. Activity score (_required if second and third columns not specified)

The second, third, and fourth columns can be used individually or together. If you only supply diplotypes then PharmCAT
will attempt to assign phenotype and activity score (when applicable). If you specify both a diplotype and a phenotype
then PharmCAT will still match the diplotype to our known phenotype mapped to that diplotype. If our mapped phenotype
and the phenotype you supply do not match then PharmCAT will emit a warning.

Different genes on different lines can mix whether they give the diplotype, phenotype, or activity values.

Here's an example of an outside call file:

```text
CYP2D6	*1/*3
CYP2C9			2.0
G6PD	B (wildtype)/B (wildtype)
HLA-B		*57:01 positive
MT-RNR1	1555A>G
```

Notes:
* the HLA-B line has two tabs between the gene name and the gene result (`*57:01 positive`)
* the MT-RNR1 line specifies a single allele call since that gene is monoploid
