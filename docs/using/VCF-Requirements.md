---
parent: Using PharmCAT
title: VCF Requirements
permalink: using/VCF-Requirements/
nav_order: 2
---
# VCF Requirements
{: .no_toc }

1. TOC
{:toc}

If you believe there is a tool out there other than bcftools that can do a better job with variant normalization (or if there is a better way to use bcftools), please [let us know](mailto:pharmcat@pharmgkb.org){: .btn .btn-blue }.

## Expected formats

PharmCAT expects incoming VCF files to follow the [official VCF spec](https://samtools.github.io/hts-specs/VCFv4.3.pdf). You can find the expected, exemplary variant representation formats in [pharmcat_positions.vcf](https://github.com/PharmGKB/PharmCAT/releases/latest).

### Requirements and expectations

PharmCAT expects incoming VCF to meet the following requirements:

1. Build version must be aligned to the __GRCh38 assembly__ (aka `b38`, `hg38`, etc.).
2. __Any position not in the input VCF is assumed to be a "no call"__. Missing positions will _not_ be interpreted as reference. You must specify all positions in the input VCF that you want to be considered.
3. Use a parsimonious, left-aligned variant representation format.
4. Have insertions and deletions normalized to the expected left-aligned representation.
5. The `CHROM` field must be in the format of __chr##__.

PharmCAT assumes the incoming VCF have the following properties:

1. The input VCF has __high-quality genotypes__. It is the user's responsibility to remove data not meeting quality criteria _before_ passing it to PharmCAT.
2. PharmCAT only requires _CHROM_, _POS_, _REF_, _ALT_, and _FORMAT/GT_ columns in a VCF file. Optional _FORMAT/AD_ and _INFO/END_ fields, when present, will be used for quality assurance checks.
3. The `QUAL` and `FILTER` columns are __not interpreted__. 


### More on parsimonious left-alignment

To avoid ambiguity in variant representation, PharmCAT uses a parsimonious, left-aligned variant representation format (as discussed in [Unified Representation of Genetic Variants](https://doi.org/10.1093/bioinformatics/btv112) by Tan, Abecasis, and Kang).

#### Deletions

PharmCAT expects deletions to be represented with an "anchoring" base at the beginning of the `REF` sequence and then the anchoring base to also appear in the `ALT` sequence. For example, the following shows a deletion of `AGAAATGGAA`:

```text
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr10	94942212	.	AAGAAATGGAA	A	.	PASS	desired-deletion-format	GT	0/1
```

as opposed to the unwanted format:

```text
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr10	94942212	.	AGAAATGGAA	.	.	PASS	do-not-want	GT	0/1
```

If the REF is a single letter it means no variant was found, so it's safe to replace it with the appropriate nucleotide string.

#### Insertions

Similarly, PharmCAT expects to find insertions with a reference base `REF="A" ALT="ATCT"`. For example, here's an insertion of `A`:

```text
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr7	99652770	rs41303343	T	TA	.	PASS	desired-insertion-format	GT	0/1
```

---

## Preparing VCF using PharmCAT's VCF Preprocessor

{: .info}
> We highly recommend that you use [PharmCAT's VCF Preprocessor](/using/VCF-Preprocessor) to prepare your VCF files for use by PharmCAT.

If you'd like to prepare your VCF yourself, the rest of this document explores the reasoning behind our requirements, some known issues, and specific examples of fulfilling them.


## Known issues and alternative solutions

### Alignment to GRCh38 Assembly

PharmCAT requires VCF files aligned to GRCh38.

While you can translate positions between assemblies, this can easily introduce errors into your data. For example, LiftOver can be problematic due to ambiguity in positions. In addition, many tools won't update the genotype.

#### Known issue with remapping

If you use [CrossMap](http://crossmap.sourceforge.net/) it will update the position and reference nucleotides. However, if the reference is now the same as the alternate, it will remove the variant from the file.

Likewise, [NCBI Remap](http://crossmap.sourceforge.net/) will not update genotypes, so even if the updated position means the variant should now be called as 0/0 it will not be updated from 1/1. Presently there are only a few positions where this matters, but this will need to be very carefully checked.

Another issue with NCBI remap is that some alternates can be lost following remap. For instance:

```text
 Pre LiftOver: 2	234668879	rs57191451	CAT	C,CATAT,CATATAT	0/3
Post LiftOver: 2	233760233	rs57191451	CAT	CATAT,CATATAT	0/3
```

We currently do not have a fix for this (NCBI has been contacted), so again you will have to check these sites carefully by hand.

#### Example: LiftOver using the GATK

This is a GRCh37-toGRCh38 LiftOver example which uses the GATK LiftoverVcf tool on the UK Biobank genotype data in the GRCh37 coordinates. The GATK LiftoverVcf tool can update genomic coordinates. If the reference allele and alternate allele are swapped in the new genome build for a single nucleotide polymorphism (SNP) locus, the GATK LiftoverVcf can reverse-complement the SNP, update the relevant genotypes (GT), and correct AF-like INFO fields. Due to the length of the example, we host the [example with detailed explanations on Google Drive](https://docs.google.com/document/d/15rxe0iG2kruEWsvBCLyNGof-YRo5T10zuQiBJkUbyJ0/edit?usp=sharing).

Note: this example is not meant to be a comprehensive documentation of solutions to all LiftOver issues. LiftOver may require additional data cleaning or preparation steps that are specific to your genomic data.


### Must have all allele-defining positions

All positions that are used to define alleles in PharmCAT must be present in the VCF file you want to run through PharmCAT, _even if they are `0/0` (reference) or `./.` (missing)_. This is different from a typical VCF file which usually only contains variant sites.

PharmCAT needs this level specificity because reference and missing alleles are interpreted differently by the allele matcher. When positions are missing from the input PharmCAT does not make assumptions about whether this is because they are reference or missing.

For one thing, we do not know what "reference" is because it can vary based on your reference sequence. Did you convert it from GRCh37 to GRCh38? If so, the "reference" from the two may not have been updated correctly and your VCF would not provide any indications that this is the case. 

Secondly, a missing entry can mean that the reference base was detected _OR_ the base was not assayed or has no call. We cannot distinguish between uncalled positions and reference in a VCF file. So we ask that you declare each required position for PharmCAT to be clear about the input.

You have to decide on how accurate you want the data you provide to PharmCAT should be, especially if you're making any clinical decisions based on PharmCAT's results. If you wish to make assumptions of your data, you are welcome to do so following the alternative solutions in this section.

The following is an exemplary GATK command for generating VCF:

```console
# gatk --java-options "-Xmx4g" HaplotypeCaller \
     --alleles pharmcat_positions.vcf -R grc38.reference.fasta -I input.bam -O output.vcf \
     -L pharmcat_positions.vcf -ip 20 --max-mnp-distance 1 --output-mode EMIT_ALL_ACTIVE_SITES
```

The `pharmcat_positions.vcf` file is available on the [PharmCAT release page](https://github.com/PharmGKB/PharmCAT/releases) on GitHub.

Please refer to the [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/27007962724507-HaplotypeCaller) documentation for details and tuning options.


### Must normalize variant representation

Variant representation is an ongoing problem in NGS ([related article](https://macarthurlab.org/2014/04/28/converting-genetic-variants-to-their-minimal-representation/)). For example, the following vcf lines all specify the same variant with different formats:

```text
chr7    117548628    .    GTTTTTTTA    GTTTTTA    .    PASS    CFTR:5T    GT    0/1
chr7    117548628    .    GTT     G    .    PASS    CFTR:5T    GT    0/1
chr7    117548628    .    G    .    .    PASS    CFTR:5T    GT    0/1
chr7    117548628    .    G(T)7A    G(T)5A    .    PASS    CFTR:5T    GT    0/1
```

Different NGS pipelines, the way files are created (for instance if a multisample file is split), and post-processing software tools all lead to these differences. As PharmCAT is directly matching these strings to what is in the definition files this can cause problems. For example, PharmCAT expects to find deletions as the ref="ATCT"  alt="A", rather than ref="TCT" alt=".". Therefore, you will need to replace all the deletions within in the file. If the ref is a single letter it means no variant was found, so it's safe to replace it with the appropriate nucleotide string.

To avoid ambiguity in variant representation, PharmCAT uses a parsimonious, left-aligned variant representation format (as discussed in [Unified Representation of Genetic Variants](https://doi.org/10.1093/bioinformatics/btv112) by Tan, Abecasis, and Kang).

We recommend performing this normalization with [bcftools](http://samtools.github.io/bcftools/bcftools.html):

```console
# bcftools norm -m+ -c ws -Oz -o output.vcf -f grc38.reference.fasta input.vcf
```

* `-m+` joins biallelic sites into multiallelic records (+)
* `-f <ref_seq_fasta>` is the reference sequence. This flag is required for normalization
* `-c ws` when incorrect or missing REF allele is encountered, warn (w) and set/fix(s) bad sites

Please consult the [bcftools documentation](http://samtools.github.io/bcftools/bcftools.html) for details.

Alternatively, you can use PharmCAT's [VCF Preprocessor](/using/VCF-Preprocessor), which also relies on bcftools for this.

It is highly recommended that you always check the output files from these tools manually to make sure the correct format normalizations have been made.


### The `CHROM` field must be in the format __chr##__

PharmCAT expects the `CHROM` field to have entries that begin with "chr" (e.g. `chr1` instead of just `1`).

[PharmCAT's VCF Preprocessor](/using/VCF-Preprocessor) takes care of this issue by automatically detecting and updating the `CHROM` field format.

Alternatively,`CHROM` field format can be updated using the following command:

```console
# perl -pe '/^((?!^chr).)*$/ && s/^([^#])/chr$1/gsi' merged_output.vcf > merged_output.chrfixed.vcf
```


### The `QUAL` and `FILTER` columns are __not interpreted__

PharmCAT considers all variants in your file, even if they fail filters.

If you use a tool like VSQR you can use the following Perl one-liner to change the genotype to 0/0 if necessary:

```console
# perl -pe '/^((?!PASS).)*$/ && /^((?!0\/0).)*$/ && /^((?!0\|0).)*$/ && s/[0-9][\/|][0-9]/0|0/' merged_output.vcf > merged_output_pass.vcf
```


### Consecutive homozygous reference genotypes will not be translated to a homozygous reference INDEL

We do not consider consecutive homozygous reference genotypes as evidence for homozygous reference INDELs. For example, one might equate the following two VCF entries

```text
# example 1
chr2    233760233       .       C       <NON_REF>       .       .       .       GT  0/0
chr2    233760234       .       A       <NON_REF>       .       .       .       GT  0/0
chr2    233760235       .       T       <NON_REF>       .       .       .       GT  0/0

# example 2
chr2    233760233       .       CAT       C       .       .       .       GT  0/0
chr2    233760233       .       C       CAT       .       .       .       GT  0/0
```

However, example 1 does not rule out the possibility of a small deletion, _e.g._, `CAT>CATAT`, while example 2 explicitly states the lack of a small AT insertion or deletion at this genomic position.

You can find more examples on [this GitHub issue page](https://github.com/samtools/bcftools/issues/2163).

It is hard to accurately infer INDELs for various cases. And as PharmCAT expects high-quality genotypes, it is beyond the scope of PharmCAT's VCF Preprocessor to perform in-depth genotype-calling tasks.

If you understand the risks, you can make the call on your data before passing it on to PharmCAT.

If you have access to the raw sequencing data, one option is to try the [force-calling feature in DeepVariant](https://github.com/google/deepvariant/issues/433). Or try out the `--alleles` feature in the [GATK HaplotypeCaller > 4.6.0.0](https://gatk.broadinstitute.org/hc/en-us/articles/27007962724507-HaplotypeCaller). See [Must have all allele-defining positions](/using/VCF-Preprocessor/#Must-have-all-allele-defining-positions) for the toy GATK command.

Alternatively, use bcftools to merge a gVCF with the PharmCAT position VCF to force the representation of INDELs before any downstream file normalization. Again, please note that PharmCAT expects high-quality VCF as the input.
```shell
# bcftools merge -m both <input_vcf.gz> pharmcat_position.vcf.bgz | bcftools view -s ^PharmCAT -Oz -o <output.vcf.gz>
```

We thank the PharmCAT users who reported the issue and shared solutions in the [original GitHub issue](https://github.com/PharmGKB/PharmCAT/issues/128)
