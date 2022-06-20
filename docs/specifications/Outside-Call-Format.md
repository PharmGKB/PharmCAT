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

## File format

The **outside call file** format is a tab-separated file. Lines starting with `#` will be ignored.

Each line has up to 3 fields, separated by tabs:

1. HGNC gene symbol (_required_)
2. Diplotype or single allele call (_required if third column not specified_)
3. Phenotype or other gene result (_required if second column not specified_)

The second and third column should not both be specified on the same line. Different genes on different lines can mix 
whether they give the diplotype and phenotype values.

Here's an example of an outside call file:

```text
CYP2D6	*1/*3
G6PD	B (wildtype)/B (wildtype)
HLA-B		*57:01 positive
MT-RNR1	1555A>G
```

Note that the HLA-B line has two tabs between the gene name and the allele presence. Also, the MT-RNR1 line specifies a 
since that gene is monoploid.
