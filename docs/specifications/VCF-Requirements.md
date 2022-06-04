---
title: VCF Requirements
parent: Specifications
nav_order: 2
---
# PharmCAT VCF Requirements

PharmCAT expects the incoming VCF files to follow the [official VCF spec](https://samtools.github.io/hts-specs/VCFv4.3.pdf).

In addition, PharmCAT expects incoming VCF to have the following properties:

1. Build version must be aligned to the __GRCh38 assembly__ (aka `b38`, `hg38`, etc.).
1. __Any position not in the input VCF is assumed to be a "no call"__. Missing positions will _not_ be interpreted as reference. You must specify all positions in the input VCF that you want to be considered.
1. Use a parsimonious, left aligned variant representation format.
1. Have insertions and deletions normalized to the expected representation.
1. The `CHROM` field must be in the format __chr##__.
1. The `QUAL` and `FILTER` columns are __not interpreted__. It is left to the user to remove data not meeting quality criteria _before_ passing it to PharmCAT.
1. Should only have data for a __single sample__.  If it's a multi-sample VCF file, __only the first sample is used__.



### Variant Representation Format

To avoid ambiguity in variant representation, PharmCAT is using a parsimonious, left-aligned variant representation format (as discussed in [Unified Representation of Genetic Variants](https://doi.org/10.1093/bioinformatics/btv112) by Tan, Abecasis, and Kang).


### Insertions & Deletions

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


## More Information

Every PharmCAT [release](https://github.com/PharmGKB/PharmCAT/releases) includes a `pharmcat_positions.vcf` VCF file that contains all positions of interest to PharmCAT.

For more details about fulfilling these requirements for PharmCAT read the [Preparing VCF Files](/specifications/Preparing-VCF-Files) page.

See [PharmCAT's VCF Preprocessor](/using/VCF-Preprocessor) for a script to automate some of these steps.
