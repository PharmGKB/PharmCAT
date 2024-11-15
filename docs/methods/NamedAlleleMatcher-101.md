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

A combination allele is when a sample matches a combination of 2 or more defined alleles.  For example, `[*6 + *14]` in
the CYP2B6 `[*6 + *14]/*13` diplotype output.

PharmCAT's syntax for combination calls uses square brackets to reflect that it is a variation on one gene copy and to
distinguish it from gene duplications (e.g. tandem arrangements like CYP2D6 `*36+*10`).

A partial allele is when a sample matches all the (core) variants of a defined allele but also has additional variants.
For example, CYP2C19 `*2/[*17 + g.94781859G>A]`.  In the case where a partial call occurs off the reference allele,
only the positions are listed (e.g. `*2/g.94781859G>A`).

When asked to find combinations and partial alleles, the `Named Allele Matcher` will only attempt to do so if no viable
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

Because PharmCAT scores on the number of matched positions in the definitions, the reference named allele (usually \*1)
will get the highest score. As such, scoring is biased towards grouping combinations together.  For example, CYP2B6
`*1/[*5 + *9 + *23]` will be the call with the highest score but permutations such as `*5/[*9 + *23]`, `*9/[*5 + *23]`,
`*23/[*5 + *9]` are also possible.





### Notes

__1__: If sample data is not phased and we do not assume the reference for missing positions in the definition, it is
possible to have multiple matches for a single named allele.

__2__: This behavior can be modified to return all potential diplotype matches.

__3__: This score is the same regardless of whether we assume the reference for missing positions.
