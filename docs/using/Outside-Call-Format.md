---
parent: Using PharmCAT
title: Outside Call Format
permalink: using/Outside-Call-Format/
nav_order: 4
---
# Outside Call Format

Typically, PharmCAT uses variant call data to match diplotypes used to find annotations. However, you can also give
diplotypes, phenotypes or other allele calls to PharmCAT that were called by other tools.

These **outside call files** can be supplied to the [PharmCAT](/using/Running-PharmCAT#outside-calls) tool using the
`-po` flag.

{: .warn}
> Calls specified in this file will override results from the `Named Allele Matcher`.


## File format

The **outside call file** format is a tab-separated file. Lines starting with `#` will be ignored.

Each line has up to 4 fields, separated by tabs:

1. HGNC gene symbol (_required_)
2. Diplotype or single allele call (_required if third and fourth columns are not specified_)
3. Phenotype or another gene result (_required if second and fourth columns are not specified_)
4. Activity score (_required if second and third columns are not specified_)

The second, third and fourth fields can be used individually or together.

Different genes on different lines can mix whether they give the diplotype, phenotype or activity values.


### Activity score genes

For [activity score genes](/methods/Matching-Recommendations/#method-2-activity-score-lookup), a provided activity score
value will trump other values.  If you also provide diplotype and/or phenotype value, a warning will be issued if they
do not match the expected values.  For example, if you specify a `*1/*6` diplotype and `4.0` activity score, a warning
will be issued because CPIC expects an activity score of `1.0` for `*1/*6`.

If you don't specify an activity score, PharmCAT will look up the activity score based on the phenotype (if provided)
and then diplotype (if phenotype is not provided).  Once again, PharmCAT will issue a warning if both phenotype and 
diplotype are provided, and they do not match the expected values. 

### Non-activity score genes

For non-activity score genes, the activity score is ignored (although it may be displayed in the final report).

If you specify both a diplotype and a phenotype, then PharmCAT will _rely on your phenotype_.
PharmCAT will issue a warning if your phenotype does not match the expected phenotype for the given diplotype.


### Example

Here's an example of an outside call file:

```text
CYP2D6	*1/*3
CYP2C9			2.0
HLA-B		*57:01 positive
MT-RNR1	1555A>G
```

Notes:
* the HLA-B line has two tabs between the gene name and the gene result (`*57:01 positive`)
* the MT-RNR1 line specifies a single allele call since that gene is monoploid


## Details

We rely on string matching to match outside calls to recommendations.

Consult the [Phenotypes List](/Phenotypes-List) for a complete list of named alleles, phenotypes and activity scores.

Prefixing allele names with the gene symbol in the second field (e.g. `CYP2C9*1/CYP2C9*3`) is not necessary. The gene is
specified in the first field, so repeating it in the second field is not necessary. Prefixed gene symbols will be
stripped from the allele names.

If there is an outside call for a gene that also has data from the VCF, the outside call will trump the VCF data.

{: .info}
> To avoid potential problems, please use [UTF-8 encoding](https://en.wikipedia.org/wiki/UTF-8) for outside call files.
>
> Pay attention to the use of special characters, especially "&ge;", which is used in CYP2D6 diplotypes and activity
> scores.


### Diplotypes

Named allele matching in diplotypes/single allele calls should be fairly straightforward.

PharmCAT will automatically convert named alleles into a format usable by PharmCAT based on existing conventions for
CYP2D6, HLA-A and HLA-B. For example:

* `HLA-B *07:02:01:127` will be truncated to `*07:02`
* `CYP2D6 *4.024` will be truncated to `*4`


#### Gene duplications and copy numbers

{: .info}
Gene duplication occurs when there are two or more copies of the same gene present on the same chromosome (in cis).
These gene copies may contain genetic variation defined by named alleles.
The gene copies may be identical (i.e. the same named allele) or non-identical (i.e. different named alleles).

PharmCAT relies on PharmVar gene definitions. CYP2D6 is the only gene with gene duplications defined.
More information about gene duplications, terminology and notation can be found in the Structural Variation document for
CYP2D6, accessible at the top of the [PharmVar CYP2D6 page](https://www.pharmvar.org/gene/CYP2D6).

Highlights include:
1. If the gene copy number is known, duplications/multiplications are recommended to be annotated and reported as `x2`,
   `x3`, etc., and if the number is unknown as `xN`.
2. To date, only CYP2D6 `*1`, `*2`, `*4`, and `*41` have been described to have 3 or more copies, whereas an increasing
   number of other star alleles have been described in the duplicated state.
3. Over 90% of CYP2D6 SVs (structural variants) involve identical copies of the gene.
4. In the past, "tandem" was used to distinguish allelic variants with two or more gene units that are not identical
   from those with identical units that are duplicated or multiplied. PharmVar no longer recommends using this term.
5. PharmVar maintains a table of non-identical CYP2D6 duplications that have been described in the literature and/or
   submitted to PharmVar.
6. PharmVar maintains a table of recommended ways to report structural variations, including gene duplications.

Syntax:
* PharmCAT follows PharmVar recommended notation for gene duplications and copy numbers.
* For identical gene copies on the same allele (in cis), the star allele is followed by an “x” and the number of
  gene copies. Gene copies may vary at the suballele level. If the number of gene copies is unknown, the star number
  is followed by "xN". Note that alleles with an unknown number of gene copies cannot be linked to prescribing
  recommendations. Examples:
    * CYP2D6 `*2x2/*4`
    * CYP2D6 `*1/*41x3`
* For non-identical gene copies on the same allele (in cis), or multiple SVs/CNVs on the same allele, the upstream
  gene copy is written first, followed by a "+" and the downstream gene copy. Although the order of the gene copies
  may not be experimentally determined in routine clinical testing, they should be displayed in their most likely order
  (the order found in PharmVar’s "Structural Variation for CYP2D6" document referred to above) for consistency.
  Examples:
    * CYP2D6 `*68+*4/*10`
    * CYP2D6 `*2/*36+*10`

In addition to only recognizing gene duplications defined in PharmVar, PharmCAT also only recognizes copy number
variations that have a function assignment from CPIC.
Consult the [CYP2D6 phenotypes list](/Phenotypes-List#cyp2d6) for the full list.
These alleles are part of the CPIC diplotype to phenotype translation and can be connected to a corresponding
recommendation.  Some copy number variations (e.g. for `*1`, `*2` and `*4`) over 3 are combined in a single bin (≥3).
So if you have `*1x3` or `*1x5`, you will need to translate that to `*1≥3`.

PharmCAT will automatically attempt to translate your CYP2D6 copy number variation into a matching CPIC copy number
variation if possible.

**IMPORTANT**: PharmCAT expects files encoded in UTF-8.
This is particularly important when it comes to the "≥" signs that are used in copy number names.


#### Combination Calls

If the call is a [combination call](/methods/NamedAlleleMatcher-101#combinations-and-partial-alleles)
(e.g. `[*2 + *3]`), it needs to use PharmCAT's combination syntax: it has to be wrapped in square brackets (`[` and `]`)
and each named allele must be separated with a plus sign with a space on either side (` + `).
More examples: `*1/[*6 + *8]`, `[*3 + *4 + *5]/[*18 + *37]`.

Please note the differences in syntax between combination calls and gene duplication calls described above!


### Phenotypes

When providing phenotypes, you will need to use
[CPIC standardized terms](https://cpicpgx.org/resources/term-standardization/), although we do provide some
interpretation:

1. We automatically normalize spelling (e.g. "metabolizer" instead of "metaboliser") and capitalization.
2. We will translate common synonyms:
    * PM = Poor Metabolizer
    * IM = Intermediate Metabolizer
    * NM = Normal Metabolizer
    * EM = Normal Metabolizer
    * UM = Ultrarapid Metabolizer

   Note that CPIC uses "normal" instead of "extensive", so the translation from "EM" above to "Normal Metabolizer" is not a mistake. 
3. We will try to extract the main phenotypes above if possible.  For example, "CYP2D6 Ultrarapid Metabolizers (UM)" becomes "Ultrarapid Metabolizer".
4. Some CPIC standardized terms for phenotypes include modifiers such as "likely" or "possible".  We retain these modifiers.  For example, "likely cyp2c19 poor metaboliser" becomes "Likely Poor Metabolizer".

### Activity Scores

String matching also applies to activity scores.  PharmCAT only recognizes CPIC assigned
activity scores.  For example, one possible CPIC activity score for CYP2D6 is "≥6.0".  If you provide "7.0", this will result in a no call.  Similarly, if CPIC defines activity scores of "0.0" and "0.25", and you provide "0.1", this will also result in a no call. 
