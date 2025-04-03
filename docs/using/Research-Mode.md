---
parent: Using PharmCAT
title: Research Mode
permalink: using/Research-Mode/
nav_order: 12
---
# Research Mode

PharmCAT has a couple of features which are not enabled by default and not recommended for typical runs but which are
included for particular use cases and should be fully understood before enabling.

These features are part of what we refer to as "Research Mode". As the name indicates, these features are not intended
for use outside a research context. They are enabled on the command line by using the `-research` flag and specifying
which research mode features you'd like to enable.

__If you specify any Research Mode feature PharmCAT will not output any data from the Reporter module__. Matcher and
Phenotyper data can still be written. Research features change the data in ways that are not compatible with making 
reliable prescribing recommendation data.

The following explains each research mode feature.


## CYP2D6 Named Allele Matching

As outlined on the [Calling CYP2D6](/using/Calling-CYP2D6) page, we **do NOT recommend calling CYP2D6 from VCF** due to
the large influence of SV and CNV on phenotype prediction.

PharmCAT can, however, call CYP2D6 star alleles that are defined based on SNPs and/or INDELs. If you are willing to
accept all these caveats, this functionality is part of Research Mode.

To get PharmCAT to call CYP2D6, use `-research cyp2d6`

```console
# java -jar pharmcat.jar -vcf patient_001.vcf -research cyp2d6
```


## Combination and Partial Allele Matches

If given the `-research combinations` flag, PharmCAT will try to call combination and partial alleles.
These are only called if an exact match to any single defined allele cannot be found.
Without this research flag these samples will yield a "not called" result from the `Named Allele Matcher`.

For more details, please see
[Combinations and Partial Alleles](/methods/NamedAlleleMatcher-101#combinations-and-partial-alleles).

To call combinations and partial alleles, use the `-research combinations` flag.

```console
# java -jar pharmcat.jar -vcf patient_001.vcf -research combinations
```


## Running multiple research mode options

Both research mode features can be run at the same time. Just write both `-research` flag values separated by a comma
to turn on both features.

```console
# java -jar pharmcat.jar -vcf patient_001.vcf -research combinations,cyp2d6
```
