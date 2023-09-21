---
parent: Using PharmCAT
title: Research Mode
permalink: using/Research-Mode/
nav_order: 10
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

If given the `-research combinations` flag, PharmCAT will try to call combination and partial alleles.  These are only called if an exact match to any single defined allele cannot be found.  Without this research flag these samples will yield a "not called" result from the `Named Allele Matcher`.

This option addresses variant combinations not catalogued by PharmVar or other nomenclature sites. It does not consider novel variants; it only considers variants included in existing allele definitions found in novel combinations.

A combination allele is when a sample matches a combination of 2 or more defined alleles.  For example, `[*6 + *14]` in the CYP2B6 `[*6 + *14]/*13` diplotype output.

PharmCAT's syntax for combination calls uses square brackets to reflect that it is a variation on one gene copy and to
distinguish it from gene duplications (e.g. tandem arrangements like CYP2D6 `*36+*10`).

A partial allele is when a sample matches all the (core) variants of a defined allele but also has additional variants.  For example, CYP2C19 `*2/[*17 + g.94781859G>A]`.  In the case where a partial call occurs off the reference allele, only the positions are listed (e.g. `*2/g.94781859G>A`).  A partial off the reference allele will only be called if the data is phased, or the unphased data only has 2 possible sequence combinations.

Note that PharmCAT only provides the match(es) with the highest score by default. Because PharmCAT scores on the number of matched positions in the definitions, the reference named allele (usually *1) will get the highest score. As such, scoring is biased towards grouping combinations together.  For example, CYP2B6 `*1/[*5 + *9 + *23]` will be the call with the highest score but permutations such as `*5/[*9 + *23]`, `*9/[*5 + *23]`, `*23/[*5 + *9]` are also valid.

For more details on combinations and partial alleles, please see [NamedAlleleMatcher 201](/methods/NamedAlleleMatcher-201#combinations-and-partial-alleles).

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
