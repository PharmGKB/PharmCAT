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

Calls specified in this file will override results from the `Named Allele Matcher`.


## File format

The **outside call file** format is a tab-separated file. Lines starting with `#` will be ignored.

Each line has up to 4 fields, separated by tabs:

1. HGNC gene symbol (_required_)
2. Diplotype or single allele call (_required if third and fourth columns not specified_)
3. Phenotype or other gene result (_required if second and fourth columns not specified_)
4. Activity score (_required if second and third columns not specified_)

The second, third, and fourth fields can be used individually or together. If you only supply diplotypes then PharmCAT
will attempt to assign phenotype and activity score where applicable. If you specify both a diplotype and a phenotype
then PharmCAT will rely on your phenotype, although it will emit a warning if your phenotype does not match the expected
phenotype for the given diplotype.

Different genes on different lines can mix whether they give the diplotype, phenotype, or activity values.

Here's an example of an outside call file:

```text
CYP2D6	*1/*3
CYP2C9			2.0
G6PD	B (wildtype)/B (wildtype)
HLA-B		*57:01 positive
MT-RNR1	1555A>G
```

Notes:
* the HLA-B line has two tabs between the gene name and the gene result (`*57:01 positive`)
* the MT-RNR1 line specifies a single allele call since that gene is monoploid


## Caveats

We rely on string matching to match outside calls to recommendations.

Consult the [Phenotypes List](/Phenotypes-List) for a complete list of named alleles, phenotypes and activity scores.


#### Diplotypes

Named allele matching in diplotypes/single allele calls should be fairly straightforward.  The biggest potential problem comes when dealing with copy numbers.  Note that PharmCAT is only aware of copy numbers in CYP2D6, and only for copy number variations that have a function assignment through CPIC: `*1x2`, `*1x≥3`, `*2x2`, `*2x≥3`, `*3x2`, `*4x2`, `*4x≥3`, `*6x2`, `*9x2`, `*10x2`, `*17x2`, `*29x2`, `*35x2`, `*36x2`, `*41x2`, `*41x3`, `*43x2`, `*45x2`. These alleles are part of the CPIC diplotype to phenotype translation and can be connected to a corresponding recommendation.  For `*1`, `*2`, and `*4`, copy numbers over 3 are combined in a single bin (≥3).  So if you have `*1x3` or `*1x5`, you will need to translate that to `*1≥3`.

**IMPORTANT**: PharmCAT expects files encoded in UTF-8.  This is particularly important when it comes to the "≥" signs that are used in copy number names.

#### Phenotypes

When providing phenotypes, you will need to match our values, although we do some interpolation.

1. We automatically normalize spelling (e.g. "metabolizer" instead of "metaboliser") and capitalization.
2. We will translate common synonyms:
    * PM = Poor Metabolizer
    * IM = Intermediate Metabolizer
    * NM = Normal Metabolizer
    * EM = Normal Metabolizer
    * UM = Ultrarapid Metabolizer

   Note that we treat "extensive" as "normal", so the translation from "EM" above to "Normal Metabolizer" is not a mistake. 
3. We will try to extract the main phenotypes above if possible.  For example, "CYP2D6 Ultrarapid Metabolizers (NM)" becomes "Ultrarapid Metabolizer".  When this happens, we keep retain "likely" and "possible" modifiers.  For example, "likely cyp2d6 extensive metaboliser" becomes "Likely Normal Metabolizer".

#### Activity Scores

String matching also applies to activity scores.  For example, one possible activity score for CYP2D6 is "≥6.0".  If you provide "7.0", this will result in a no call.  Simlarly, if we look for an activity score of "0.0" and "0.25", and you provide "0.1", this will also result in a no call. 
