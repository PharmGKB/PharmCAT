---
parent: How It Works
title: Named Allele Matcher 101
permalink: methods/NamedAlleleMatcher-101/
nav_order: 2
---
# Named Allele Matcher 101

The `Named Allele Matcher` is responsible for calling diplotypes from variant call data.  While it is designed to be
used in the PharmCAT pipeline, it can also be run independently.

The basic process:

1. Read in all named allele definitions from the gene definition table.
   Each gene has a reference allele defined by the first definition row in the table (e.g. \*1).  By default, any
   non-reference named allele that does not contain a base call for a given position (i.e. blank spots in the definition
  table) will default to the reference row's base call.
2. Read in sample data (VCF file), ignoring positions that are not used in the gene definition tables.
3. For each gene:
    * If the data is phased:
        1. Attempt to match the genotypes in each strand to a named allele.
        2. Return matched diplotype
    * If the data is unphased:
        1. Generate all possible combinations of genotypes for the positions of interest.
        2. Attempt to match each combination to a named allele.<sup>[1](#notes)</sup>
        3. If there are matches, try to build diplotypes by making sure that the genotype combinations are possible.
        4. Return matched diplotype(s).  If multiple diplotypes are possible, they are scored and only the top-scoring
           diplotype(s) is returned.<sup>[2](#notes)</sup>

Note that some genes have to be handled differently (e.g. DPYD and SLCO1B1).
See [Gene Definition Exceptions](/methods/Gene-Definition-Exceptions) for details.


{: .info }
When given unphased data, the `Named Allele Matcher` may call multiple diplotypes. These are all equally possible and
valid in the absence of additional information. They are assigned scores not because one is better than another, but so
that we can consistently pick one over another when necessary since we are only reporting those with the highest score
by default.<sup>[2](#notes)</sup> It is as arbitrary as sorting them alphabetically and defaulting to the first one.
<br /><br />
_For reliable diplotype calls, phased data is best_.


### Structural Variants and Multi-allelic results 

The `Named Allele Matcher` does not currently support structural variants, including gene copy number.
If structural variants are detected in the VCF data, it will be ignored and a warning will be issued.

If it detects more than the expected number of alleles in the `GT` column of the VCF, only the first two alleles will be
used and a warning will be issued.  On haploid chromosomes, only the first allele will be used.


### Partial Genotypes

If the genotype field (GT) in the VCF lists a partial genotype (e.g. `0/.` or `.|1`), the `Named Allele Matcher` will
check if it's on a haploid or diploid chromosome. For diploid chromosomes, if the partial is the reference (`0`), the
position will be ignored and treated as a missing position.  If the partial is an alternate allele (ALT), then the
partial will be accepted as-is.  However, this means that the `Named Allele Matcher` will not be able to call
a diplotype (unless the combination research mode is enabled).


## Scoring

Each named allele is given a score based on the number of variant positions used to define the allele (non-blank cells
in that row).  This means that the reference allele will always have the maximum score because all positions are defined
for that allele.

Take a look at this sample gene definition table:

|     | rs1 | rs2 | rs3 | rs4 | rs5 | score |
| --- | --- | --- | --- | --- | --- | ----- |
| \*1 | C   | C   | T   | G   | A   | 5     |
| \*2 | T   | T   |     | A   |     | 3     |

Since the gene definition table contains 5 positions, the reference allele, \*1, gets a score of 5 while \*2 only has 3
positions defined and gets a score of 3. <sup>[3](#notes)</sup>

A diplotype's score is the combined score of its component named alleles.  A \*1/\*2 from the example above would have a
score of 8.


### Missing Positions

If the sample data has missing positions that are required by a named allele definition, the position will be dropped
from consideration.

This is the only reason the score for a diplotype might be different between two samples.


### Examples

Using the following gene definition table:

|     | rs1 | rs2 | rs3 | rs4 | rs5 | score |
| --- | --- | --- | --- | --- | --- | ----- |
| \*1 | C   | C   | T   | G   | A   | 5     |
| \*2 | T   |     |     | A   |     | 2     |
| \*3 |     |     |     | A   |     | 1     |
| \*4 |     | T   |     |     |     | 1     |
| \*5 | T   |     |     |     |     | 1     |

And the following (unphased) sample data:

| rs1 | rs2 | rs3 | rs4 | rs5 |
| --- | --- | --- | --- | --- |
| C/T | C/C | T/T | A/G | A/A |

The potential permutations of those genotypes will match \*1, \*2, \*3 and \*5.

From that, plausible diplotypes are:

| Diplotype | Score     |
| --------- | --------- |
| \*1/\*2   | 5 + 2 = 7 |
| \*3/\*5   | 1 + 1 = 2 |

Which results in \*1/\*2 being returned.

Note that \*1/\*3 is not a plausible diplotype because one chromosome must have a `C` and the other must have a `T` at
position rs1.  They can't both be `C`.  Similarly, \*2/\*3 is not a plausible diplotype either because it cannot be 
homozygous at rs4.


#### Missing rs5

If we use the same gene definition table as the example above and the following (unphased) sample data, _with no data available for rs5_:

| rs1 | rs2 | rs3 | rs4 | rs5 |
| --- | --- | --- | --- | --- |
| C/T | C/C | T/T | A/G |     |

The results would be the same, except the scores would be different:

| Diplotype | Score     |
| --------- | --------- |
| \*1/\*2   | 4 + 2 = 6 |
| \*3/\*5   | 1 + 1 = 2 |


#### Missing rs1

If we use the same gene definition table as the example above and the following (unphased) sample data, _with no data available for rs1_:

| rs1 | rs2 | rs3 | rs4 | rs5 |
| --- | --- | --- | --- | --- |
|     | C/C | T/T | A/G | A/A |

Then the results would be different:

| Diplotype | Score     |
| --------- | --------- |
| \*1/\*2   | 4 + 1 = 5 |
| \*1/\*3   | 4 + 1 = 5 |

As such, \*1/\*2 and \*1/\*3 would be returned.

Note that \*5 could never be matched in this scenario because its defining allele is missing.


## Undocumented Variations

By default, only genetic variations that are defined in the allele definitions can be mapped to genotypes by the
`Named Allele Matcher`. If the sample includes a variant call that is located at an allele-defining position but itself
not included in the allele definitions, the `Named Allele Matcher` produces a "Not called" for the affected gene since
the sample matches neither the reference nor any defined variant.

A "Not called" output cannot be connected to guideline recommendations, even if the sample has other defined, actionable
variants. The decision was made that, in the interest of providing recommendation guidance, these undocumented variants
are set to reference for genes for which the defined variants affect drug toxicity. This applies to:

* CACNA1S
* DPYD
* G6PD
* NUDT15
* RYR1
* TPMT


## Combinations and Partial Alleles

{: .warn}
> Calling combinations and partial alleles is intended for **research use only**.

A combination allele is when a sample matches a combination of 2 or more named alleles _in the same single gene copy_.
For example, `[*6 + *14]` in the CYP2B6 `[*6 + *14]/*13` diplotype output.

A partial allele is when a sample matches all the (core) variants of a defined allele but also has additional variants.
For example, CYP2C19 `*2/[*17 + g.94781859G>A]`.  In the case where a partial call occurs off the reference allele,
only the positions are listed (e.g. `*2/g.94781859G>A`).

When asked to find combinations and partial alleles, the `Named Allele Matcher`
* will only attempt to do so if a normal call cannot be made.
* will only look for combinations catalogued by PharmVar or other nomenclature sites.
It does not consider novel variants (i.e. variants at undocumented positions); it only considers variants included in
existing allele definitions found in novel combinations.

When looking for combinations, the variations of every named allele must not overlap – the same SNV cannot be counted 
for more than one star (*) allele. Examples:

1. A combination such as CYP2D6 `[*10 + *37]` could not exist because `100C>T` and `4181G>C` are present on both named
   alleles, and each SNV can only be assigned to one named allele in a single gene copy.
2. On the other hand CYP2D6 `[*10 + *25]` could theoretically exist because the core SNV in `*25` is `3199C>G` and the
   core SNVs in `*10` are `100C>T` and `4181G>C`.

{: .info}
> PharmCAT's syntax for combinations uses square brackets to surround the alleles and will always have spaces around the
> `+`.
>
> Examples:
> * `*2/[*17 + g.94781859G>A]`
> * `[*6 + *14]/*13`
> * `[*6 + *14 + g.94781859G>A]/*13`

PharmCAT's syntax for combinations uses square brackets to reflect that it is a variation on one gene copy and to
distinguish it from gene duplications (e.g. tandem arrangements like CYP2D6 `*36+*10`).


### Scoring

Because PharmCAT scores on the number of matched positions in the definitions, the reference named allele (usually \*1)
will get the highest score. As such, scoring is biased towards grouping combinations together.  For example, CYP2B6
`*1/[*5 + *9 + *23]` will be the call with the highest score but permutations such as `*5/[*9 + *23]`, `*9/[*5 + *23]`,
`*23/[*5 + *9]` are also possible.


### Unphased Data

When dealing with unphased data, note that the `Named Allele Matcher` will never produce a combination/partial allele
call if there is a viable non-combination/partial call.  For example, if the sample is CYP2B6 `*1/[*8 + *9]` with phased
data, the `Named Allele Matcher` will only call `*8/*9` if the data is unphased because that is a viable call.
It will not attempt to look for potential combination/partial alleles.

Once the `Named Allele Matcher` starts looking for combinations, it will look for all occurrences of named alleles in
all possible permutations of the unphased data.  It will only attempt to find partials based on the reference sequence
if there are zero or one named allele matches.  The matching diplotypes will then be scored like normal.   

Examples:

1. If the SNVs defining CYP2C9 `*2` and CYP2C9 `*3` are found alone in unphased data, the predicted diplotype will be 
   CYP2C9 `*2/*3`.
2. If the SNV defining CYP2C9 `*2` is found along with 1 SNV defining `*18` (`g.94986073A>C`, rs72558193), the predicted
   diplotype will be CYP2C9 `*1/[*2 + g.94986073A>C]` as the highest scoring result and CYP2C9 `*2/g.94986073A>C` as an
   equal possible option but with a lower matching score.
3. If the SNV defining CYP2C9 `*2` is found to be homozygous and the SNV defining CYP2C9 `*3` is also found, the 
   predicted diplotype will be CYP2C9 `*2/[*2 + *3]`.



## Notes

__1__: If sample data is not phased and we do not assume the reference for missing positions in the definition, it is
possible to have multiple matches for a single named allele.

__2__: This behavior can be modified to return all potential diplotype matches.

__3__: This score is the same regardless of whether we assume the reference for missing positions.
