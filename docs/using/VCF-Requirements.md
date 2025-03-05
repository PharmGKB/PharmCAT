---
parent: Using PharmCAT
title: VCF Requirements
permalink: using/VCF-Requirements/
nav_order: 3
---
# VCF Requirements

PharmCAT expects incoming VCF files to follow the [official VCF spec](https://samtools.github.io/hts-specs/VCFv4.3.pdf).
You can find the variants PharmCAT requires and their expected representation in [pharmcat_positions.vcf](https://github.com/PharmGKB/PharmCAT/releases/latest).

PharmCAT only requires `CHROM`, `POS`, `REF`, `ALT`, and `FORMAT/GT` columns in a VCF file.
PharmCAT supports phase sets `FORMAT/PS` as of version 3.0.
Optional `FORMAT/AD` and `INFO/END` fields, when present, will be used for simple validation checks.

VCF should only contain _high-quality genotypes_.
It is the user's responsibility to remove data not meeting quality criteria _before_ passing it to PharmCAT.
As such, the `QUAL` and `FILTER` columns are not interpreted.<sup>[1](#1-the-qual-and-filter-columns-are-not-interpreted)</sup>

In addition, PharmCAT expects incoming VCF to meet the following requirements:

1. Build version must be [aligned to the __GRCh38 assembly__](#requirement-1---alignment-to-grch38-assembly) (aka `b38`, `hg38`, etc.).
2. [Specify all required positions](#requirement-2---specify-all-allele-defining-positions) you want to be considered.
    __Any position not in the input VCF is assumed to be a "no call"__.
    Missing positions will _not_ be interpreted as reference. 
3. [Use a parsimonious, left-aligned variant representation format](#requirement-3---use-parsimonious-left-aligned-variant-representation).
4. [Normalize insertions and deletions](#requirement-4---normalize-insertions-and-deletions) to the expected left-aligned representation.
5. [The `CHROM` field must be in the format "chr##"](#requirement-5---the-chrom-field-must-be-in-the-format-chr).

{: .attention}
We highly recommend that you use PharmCAT's [VCF Preprocessor](/using/VCF-Preprocessor) to prepare your VCF files for
use by PharmCAT.

The PharmCAT VCF Preprocessor automates the process of making your VCF PharmCAT-ready as much as possible.
It will handle requirements #2 - #5 in most cases.  However, it will always err on the side of caution and not make
any assumptions about your data.  If it cannot automate something, it's most likely because there is no obvious way for
it to do so, and you will have to manage the idiosyncrasies in your data yourself.

For real-world examples of changes the VCF Preprocessor will make, we have a page collecting
[before and after examples](/using/Variant-Normalization) of these normalization changes.

If you believe there is a tool that can do a better job with variant normalization other than bcftools
(or if there is a better way to use bcftools), please [let us know](mailto:pharmcat@pharmgkb.org).


## Deep Dive

The rest of this document explores the reasoning behind our requirements, some known issues, and specific examples of
fulfilling them.


### Requirement #1 - Alignment to GRCh38 Assembly

PharmCAT requires VCF files aligned to GRCh38. The VCF Preprocessor is unable to handle this automatically because it
is a non-trivial task, and you will need to make decisions that are dependent on your data.

While you can translate positions between assemblies, this can easily introduce errors into your data. For example,
LiftOver can be problematic due to ambiguity in positions. In addition, many tools won't update the genotype.

For example, if you use [CrossMap](http://crossmap.sourceforge.net/) it will update the position and reference
nucleotides. However, if the reference is now the same as the alternate, it will remove the variant from the file.
As a result, you will have to check these sites carefully by hand.

#### Example: LiftOver using the GATK

Due to the length of the example and its detailed explanations, it is being hosted
[on Google Docs](https://docs.google.com/document/d/15rxe0iG2kruEWsvBCLyNGof-YRo5T10zuQiBJkUbyJ0).

This is a GRCh37-toGRCh38 LiftOver example that uses the GATK LiftoverVcf tool on the UK Biobank genotype data in the
GRCh37 coordinates. The GATK LiftoverVcf tool can update genomic coordinates. If the reference allele and alternate
allele are swapped in the new genome build for a single nucleotide polymorphism (SNP) locus, the GATK LiftoverVcf can
reverse-complement the SNP, update the relevant genotypes (GT), and correct AF-like INFO fields.

Note: this example is not meant to be a comprehensive documentation of solutions to all LiftOver issues.
LiftOver may require additional data cleaning or preparation steps that are specific to your genomic data.


### Requirement #2 - Specify all allele-defining positions

All positions that are used to define alleles in PharmCAT must be present in the VCF file.
You can find the variants PharmCAT requires and their expected representation in
[pharmcat_positions.vcf](https://github.com/PharmGKB/PharmCAT/releases/latest).

The VCF Preprocessor requires even positions that are reference (`0/0`) and missing (`./.`) to be specified as such.
This is different from a typical VCF file, which usually only contains variant sites. A missing entry can mean that 
the reference base was detected _OR_ the base was not assayed or has no call.  PharmCAT will not make any assumptions,
so we ask that you declare each required position for PharmCAT so that there is no confusion.

In addition, we do not know what "reference" is because it can vary based on your reference sequence.
Did you convert it from GRCh37 to GRCh38? If so, the "reference" from the two may not have been updated correctly,
and your VCF would not provide any indication that this is the case.

You have to decide on how accurate you want the data you provide to PharmCAT should be, especially if you're making any
clinical decisions based on PharmCAT's results. If you wish to make assumptions on your data, you are welcome to do so.

If you have .bam files to work from, one possible approach is to use the following GATK command to generate the VCF: 

```console
# gatk --java-options "-Xmx4g" HaplotypeCaller \
     --alleles pharmcat_positions.vcf -R grc38.reference.fasta -I input.bam -O output.vcf \
     -L pharmcat_positions.vcf -ip 20 --max-mnp-distance 1 --output-mode EMIT_ALL_ACTIVE_SITES
```

This is expected to, to our best knowledge at the time of writing, call all the genetic variants specified by the 
`--alleles pharmcat_positions.vcf` from your input bam file. Please refer to the [GATK documentation](https://gatk.broadinstitute.org/hc/en-us/categories/360002369672-Tool-Index) for the latest updates.

The `pharmcat_positions.vcf` file is available on the [PharmCAT release page](https://github.com/PharmGKB/PharmCAT/releases) on GitHub.

Please refer to the [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/27007962724507-HaplotypeCaller) documentation for details and tuning options.


### Requirement #3 - Use parsimonious, left-aligned variant representation

#### The problem

Variant representation is an ongoing problem in NGS ([related article](https://macarthurlab.org/2014/04/28/converting-genetic-variants-to-their-minimal-representation/)).
For example, the following VCF lines all specify the same variant with different formats:

```text
chr7    117548628    .    GTTTTTTTA    GTTTTTA    .    PASS    CFTR:5T    GT    0/1
chr7    117548628    .    GTT     G    .    PASS    CFTR:5T    GT    0/1
chr7    117548628    .    G    .    .    PASS    CFTR:5T    GT    0/1
chr7    117548628    .    G(T)7A    G(T)5A    .    PASS    CFTR:5T    GT    0/1
```

Different NGS pipelines, the way VCF files are created (e.g. if a multi-sample file is split or not), and
post-processing software tools all lead to these differences. This can cause problems since PharmCAT is directly
matching these strings to what is in the definition files.

For example, PharmCAT expects to find deletions where the ref="ATCT" and alt="A", rather than ref="TCT" and alt=".".

Therefore, you will need to replace all the deletions within in the file. If the ref is a single letter, it means no
variant was found, so it's safe to replace it with the appropriate nucleotide string.


#### The solution

To avoid ambiguity in variant representation, PharmCAT uses a parsimonious, left-aligned variant representation format
(as discussed in [Unified Representation of Genetic Variants](https://doi.org/10.1093/bioinformatics/btv112) by Tan, Abecasis, and Kang).

PharmCAT's [VCF Preprocessor](/using/VCF-Preprocessor) takes care of this automatically.

Under the hood, we are using [bcftools](http://samtools.github.io/bcftools/bcftools.html).
You can perform this step yourself with:

```console
# bcftools norm -m+ -c ws -Oz -o output.vcf -f grc38.reference.fasta input.vcf
```

* `-m+` joins biallelic sites into multiallelic records (+)
* `-f <ref_seq_fasta>` is the reference sequence. This flag is required for normalization
* `-c ws` when incorrect or missing REF allele is encountered, warn (w) and set/fix(s) bad sites

Please consult the [bcftools documentation](http://samtools.github.io/bcftools/bcftools.html) for details. It is highly
recommended that you always check the output files from these tools manually to make sure the correct format
normalizations have been made.


### Requirement #4 - Normalize insertions and deletions

This is really a recapitulation of requirement #3 - to use a parsimonious, left-aligned variant representation format
and make sure it's applied to insertions and deletions as well as SNPs.

#### Deletions

PharmCAT expects deletions to be represented with an "anchoring" base at the beginning of the `REF` sequence and then
the anchoring base to also appear in the `ALT` sequence. For example, the following shows a deletion of `AGAAATGGAA`:

```text
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr10	94942212	.	AAGAAATGGAA	A	.	PASS	desired-deletion-format	GT	0/1
```

as opposed to the unwanted format:

```text
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr10	94942212	.	AGAAATGGAA	.	.	PASS	do-not-want	GT	0/1
```

If the REF is a single letter, it means no variant was found, so it's safe to replace it with the appropriate nucleotide
string.

#### Insertions

Similarly, PharmCAT expects to find insertions with a reference base `REF="A" ALT="ATCT"`. For example, here's an
insertion of `A`:

```text
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr7	99652770	rs41303343	T	TA	.	PASS	desired-insertion-format	GT	0/1
```

#### INDELs

##### Consecutive homozygous reference genotypes will not be translated to a homozygous reference INDEL

We do not consider consecutive homozygous reference genotypes as evidence for homozygous reference INDELs.
For example, one might equate the following two VCF entries:

```text
# example 1
chr2    233760233       .       C       <NON_REF>       .       .       .       GT  0/0
chr2    233760234       .       A       <NON_REF>       .       .       .       GT  0/0
chr2    233760235       .       T       <NON_REF>       .       .       .       GT  0/0

# example 2
chr2    233760233       .       CAT       C       .       .       .       GT  0/0
chr2    233760233       .       C       CAT       .       .       .       GT  0/0
```

However, example 1 does not rule out the possibility of a small deletion (i.e. `CAT>CATAT`), while example 2 explicitly
states the lack of a small AT insertion or deletion at this genomic position.

You can find more examples in [bcftools issue #2163](https://github.com/samtools/bcftools/issues/2163).

It is hard to accurately infer INDELs for all cases. And as PharmCAT expects high-quality genotypes, it is beyond
the scope of PharmCAT's VCF Preprocessor to perform in-depth genotype-calling tasks.

If you understand the risks, you can make the call on your data before passing it on to PharmCAT.

If you have access to the raw sequencing data, you can

* try the [force-calling feature in DeepVariant](https://github.com/google/deepvariant/issues/433)
* the `--alleles` feature in the [GATK HaplotypeCaller > 4.6.0.0](https://gatk.broadinstitute.org/hc/en-us/articles/27007962724507-HaplotypeCaller). See [Specify all allele-defining positions](#requirement-2---specify-all-allele-defining-positions) for a sample GATK command.
* use bcftools to merge a gVCF with pharmcat_positions.vcf to force the representation of INDELs before
any downstream file normalization. Again, please note that PharmCAT expects high-quality VCF as the input:
```shell
# bcftools merge -m both <input_vcf.gz> pharmcat_position.vcf.bgz | bcftools view -s ^PharmCAT -Oz -o <output.vcf.gz>
```

We thank the PharmCAT users who reported the issue and shared solutions in the [original GitHub issue](https://github.com/PharmGKB/PharmCAT/issues/128).




### Requirement #5 - The `CHROM` field must be in the format "chr##"

PharmCAT expects the `CHROM` field to have entries that begin with "chr" (e.g. `chr1` instead of just `1`).

PharmCAT's [VCF Preprocessor](/using/VCF-Preprocessor) takes care of this issue by automatically detecting and updating the `CHROM` field format.

Alternatively,`CHROM` field format can be updated using the following command:

```console
# perl -pe '/^((?!^chr).)*$/ && s/^([^#])/chr$1/gsi' merged_output.vcf > merged_output.chrfixed.vcf
```


### Miscellaneous

#### <sup>1</sup> The `QUAL` and `FILTER` columns are not interpreted

PharmCAT considers all variants in your file, even if they fail filters.

If you use a tool like VSQR you can use the following Perl one-liner to change the genotype to 0/0 if necessary:

```console
# perl -pe '/^((?!PASS).)*$/ && /^((?!0\/0).)*$/ && /^((?!0\|0).)*$/ && s/[0-9][\/|][0-9]/0|0/' merged_output.vcf > merged_output_pass.vcf
```
