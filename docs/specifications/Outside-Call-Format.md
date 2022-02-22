---
title: Outside Call Format
permalink: specifications/outside-call-format/
parent: Specifications
nav_order: 4
---
# Outside Call Format

Typically, PharmCAT uses variant call data to match diplotypes used to find annotations. However, you can also give 
diplotypes to PharmCAT that were called by other tools. This is especially useful for genes that PharmCAT will not match 
like CYP2D6.

These **outside call files** can be supplied to the overall [PharmCAT](/technical-docs/Running-PharmCAT) tool using the `-a` flag or to the 
[Phenotyper](/technical-docs/Running-PharmCAT#running-the-phenotyper) using the `-o` flag.

## File format

The **outside call file** format is a tab-separated file.

Each line starts with the HGNC gene symbol, then a tab, then the diplotype (e.g. `*1/*3`). Lines starting with `#` will 
be ignored.

Here's an example of an outside call file:

```text
CYP2D6	*1/*3
G6PD	B (wildtype)/B (wildtype)
MT-RNR1	1555A>G
```
