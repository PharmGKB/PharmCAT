---
title: Named Allele Matcher 101
permalink: methods/NamedAlleleMatcher-101/
parent: Methods
---
# NamedAlleleMatcher 101

The basic process:

1. Read in all named allele definitions from the gene definition table.
   Each gene has a reference allele defined by the first definition row in the table (eg. _*1_).  By default, any non-reference named allele that does not contain a base call for a given position (eg. blank spots in the definition table) will default to the reference row's base call.
2. Read in sample data (VCF file), ignoring positions that are not used in the gene definition tables.
3. For each gene:
    1. If data is unphased, generate all possible combinations of genotypes for the positions of interest.
    2. Attempt to match each combination to a named allele.<sup>[1](#notes)</sup>
    3. If there are matches, and data is unphased, try to build diplotypes by making sure that the genotype combinations are possible.
    4. If diplotypes can be found, they are scored and only the top-scoring diplotype(s) is returned.<sup>[2](#notes)</sup>




### Scoring

Each named allele is given a score based on the number of variant positions used to define the allele (non-blank cells in that row).  This means that the reference allele will always have the maximum score because all positions are defined for that allele.

Take a look at this sample gene definition table:

|     | rs1 | rs2 | rs3 | rs4 | rs5 | score |
| --- | --- | --- | --- | --- | --- | ----- |
| *1  | C   | C   | T   | G   | A   | 5     |
| *2  | T   | T   |     | A   |     | 3     |

Since the gene definition table contains 5 positions, the reference allele, _*1_, gets a score of 5 while _*2_ only has 3 positions defined and gets a score of 3. <sup>[3](#notes)</sup>

A diplotype's score is the combined score of its component named alleles.  A _*1_/_*2_ from the example above would have a score of 8.


### Missing Positions

If the sample data has missing positions that are required by a named allele definition, the position will be dropped from consideration.

This is the only reason the score for a diplotype might be different between two samples.


### Examples

Using the following gene definition table:

|     | rs1 | rs2 | rs3 | rs4 | rs5 | score |
| --- | --- | --- | --- | --- | --- | ----- |
| *1  | C   | C   | T   | G   | A   | 5     |
| *2  | T   |     |     | A   |     | 2     |
| *3  |     |     |     | A   |     | 1     |
| *4  |     | T   |     |     |     | 1     |
| *5  | T   |     |     |     |     | 1     |

And the following (unphased) sample data:

| rs1 | rs2 | rs3 | rs4 | rs5 |
| --- | --- | --- | --- | --- |
| C/T | C/C | T/T | A/G | A/A |

The potential permutations of those genotypes will match _*1_, _*2_, _*3_ and _*5_.

From that, plausible diplotypes are:

| Diplotype | Score     |
| --------- | --------- |
| *1/*2     | 5 + 2 = 7 |
| *3/*5     | 1 + 1 = 2 |

Which results in _*1_/_*2_ being returned.

Note that _*1_/_*3_ is not a plausible diplotype because one chromosome must have a `C` and the other must have a `T` at position rs1.  They can't both be `C`.  Similarly, _*2_/_*3_ is not a plausible diplotype either because it cannot be homozygous at rs4.


#### Missing rs5

If we use the same gene definition table as the example above and the following (unphased) sample data, _with no data available for rs5_:

| rs1 | rs2 | rs3 | rs4 | rs5 |
| --- | --- | --- | --- | --- |
| C/T | C/C | T/T | A/G |     |

The results would be the same, except the scores would be different:

| Diplotype | Score     |
| --------- | --------- |
| *1/*2     | 4 + 2 = 6 |
| *3/*5     | 1 + 1 = 2 |


#### Missing rs1

If we use the same gene definition table as the example above and the following (unphased) sample data, _with no data available for rs1_:

| rs1 | rs2 | rs3 | rs4 | rs5 |
| --- | --- | --- | --- | --- |
|     | C/C | T/T | A/G | A/A |

Then the results would be different:

| Diplotype | Score     |
| --------- | --------- |
| *1/*2     | 4 + 1 = 5 |
| *1/*3     | 4 + 1 = 5 |

As such, _*1_/_*2_ and _*1_/_*3_ would be returned.

Note that *5 could never be matched in this scenario because it's definining allele is missing.


### Notes

__1__:  If sample data is not phased and we do not assume the reference for missing positions in the definition, it is possible to have multiple matches for a single named allele.

__2__: This behavior can be modified to return all potential diplotype matches.

__3__: This score is the same irregardless of whether we assume the reference for missing positions.
