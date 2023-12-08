---
parent: Methods
title: Named Allele Matcher 201
permalink: methods/NamedAlleleMatcher-201/
nav_order: 5
---
# Named Allele Matcher 201

This document attempts to cover more advanced details of the `Named Allele Matcher`.

This is really only applicable if you are defining your own named allele definitions.


## Scoring

The score for a given named allele is determined by the number of positions at which alleles have been provided.

Example:

|     | rs1 | rs2 | rs3 | rs4 | rs5 | score |
| --- | --- | --- | --- | --- | --- | ----- |
| *1  | C   | C   | T   | G   | A   | 5     |
| *2  | T   | T   |     | A   |     | 3     |

The default named allele definitions are designed to assume that missing alleles are the same as the reference named 
allele, which is defined as the first named allele in the definition file.

It is, however, possible to increase the score of a named allele by specifying the reference allele.  For example, this
gene definition table is effectively identical to the one above, but _*2_ has a different score.

|     | rs1 | rs2 | rs3 | rs4 | rs5 | score |
| --- | --- | --- | --- | --- | --- | ----- |
| *1  | C   | C   | T   | G   | A   | 5     |
| *2  | T   | T   | T   | A   | A   | 5     |



## Exemptions

`src/main/resources/org/pharmgkb/pharmcat/definition/alleles/exemptions.json` gives you a way to modify the behavior of
the `Named Allele Matcher`.


### Ignoring Named Alleles

If you are designing your own named allele definitions, you might need to define a named allele but not want it to be
considered by the `Named Allele Matcher`.

You can add an exemption for this using `ignoredAlleles` and `ignoredAllelesLc` (the latter is just a lower-cased
collection of the former).

```json
  {
    "gene": "XXX",
    "ignoredAlleles": [
      "*1S"
    ],
    "ignoredAllelesLc": [
      "*1s"
    ]
  }
```


## Combinations and Partial Alleles

{: .warn}
> Calling combination and partial alleles is intended for **research use only**.

A combination allele is when a sample matches a combination of 2 or more defined alleles.  For example, `[*6 + *14]` in
the CYP2B6 `[*6 + *14]/*13` diplotype output.

PharmCAT's syntax for combination calls uses square brackets to reflect that it is a variation on one gene copy and to
distinguish it from gene duplications (e.g. tandem arrangements like CYP2D6 `*36+*10`).

A partial allele is when a sample matches all the (core) variants of a defined allele but also has additional variants.
For example, CYP2C19 `"*2/[*17 + g.94781859G>A]"`.  In the case where a partial call occurs off the reference allele,
only the positions are listed (e.g. `*2/g.94781859G>A`).

When asked to find combination and partial alleles, the `Named Allele Matcher` will only attempt to do so if no viable 
call can be made. 

The `Named Allele Matcher` will only look for variant combinations not catalogued by PharmVar or other nomenclature 
sites. It does not consider novel variants; it only considers variants included in existing allele definitions found
in novel combinations.


### Phased vs. Unphased Data

When dealing with unphased data, note that the `Named Allele Matcher` will never produce a combination/partial allele
call if there is a viable non-combination/partial call.  For example, if the sample would be called CYP2B6
`*1/[*8 + *9]` with phased data, the `Named Allele Matcher` will only call `*8/*9` if the data is unphased because that
is a viable call.  It will not attempt to look for potential combination/partial alleles.

In addition, to limit the potential search space, a partial off the reference allele will only be called if the data is
phased or the unphased data only has 2 possible sequence combinations.


### Scoring

Because PharmCAT scores on the number of matched positions in the definitions, the reference named allele (usually *1) 
will get the highest score. As such, scoring is biased towards grouping combinations together.  For example, CYP2B6 
`*1/[*5 + *9 + *23]` will be the call with the highest score but permutations such as `*5/[*9 + *23]`, `*9/[*5 + *23]`,
`*23/[*5 + *9]` are also possible.


## Undocumented Variations

By default, only genetic variations that are defined in the allele definitions can be mapped to genotypes by the 
`Named Allele Matcher`. If the sample includes a variant call that is not included in the allele definitions at a
position that is part of the allele definitions, the `Named Allele Matcher` produces a "Not called" for the affected
gene since the sample matches neither the reference nor any defined variant.

A different approach is taken for genes for which the defined variants affect drug toxicity.
When these "undefined" variant calls are encountered, they will be treated as reference.
This applies to:

* CACNA1S
* DPYD
* G6PD
* NUDT15
* RYR1
* TPMT
