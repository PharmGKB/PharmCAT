---
parent: How It Works
title: Named Allele Matcher 201
permalink: methods/NamedAlleleMatcher-201/
nav_order: 3
---
# Named Allele Matcher 201

This document is _mainly for those interested in defining your own named allele definitions_ and covers some of the more
advanced details of the `Named Allele Matcher` .



## Scoring

The score for a given named allele is determined by the number of positions at which alleles have been provided.

Example:

|     | rs1 | rs2 | rs3 | rs4 | rs5 | score |
| --- | --- | --- | --- | --- | --- | ----- |
| \*1 | C   | C   | T   | G   | A   | 5     |
| \*2 | T   | T   |     | A   |     | 3     |

The default named allele definitions are designed to assume that missing alleles are the same as the reference named 
allele, which is defined as the first named allele in the definition file.

It is, however, possible to increase the score of a named allele by specifying the reference allele.  For example, this
gene definition table is effectively identical to the one above, but `*2` has a different score.

|     | rs1 | rs2 | rs3 | rs4 | rs5 | score |
| --- | --- | --- | --- | --- | --- | ----- |
| \*1 | C   | C   | T   | G   | A   | 5     |
| \*2 | T   | T   | T   | A   | A   | 5     |



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
