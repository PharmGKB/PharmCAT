---
parent: FAQs
title: Instances of Normalizing VCF Variant Representation Formats
permalink: Variants-Normalization/
---

This page presents a variety of instances of VCF variant representation formats that should be normalized and standardized for PharmCAT.

You should be able to concatenate different instances of VCF records into a single VCF file to test the robustness of a VCF preprocessor tool.


# VCF Header

This is the header you will need for every VCF file.
 ```VCF header
##fileformat=VCFv4.3
##source=PharmCAT allele definitions
##reference=hg38
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=FT,Number=.,Type=String,Description="Genotype filters.">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype Likelihoods">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block.">
##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions.">
##contig=<ID=chr1,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr2,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr3,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr4,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr5,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr6,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr7,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr8,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr9,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr10,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr11,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr12,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr13,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr14,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr15,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr16,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr17,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr18,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr19,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr20,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr21,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chr22,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chrX,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chrY,assembly=GRCh38.p13,species="Homo sapiens">
##contig=<ID=chrM,assembly=GRCh38.p13,species="Homo sapiens">
##INFO=<ID=PX,Number=.,Type=String,Description="Gene">
##INFO=<ID=POI,Number=0,Type=Flag,Description="Position of Interest but not part of an allele definition">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=PCATxREF,Description="Reference allele does not match PharmCAT reference alleles">
##FILTER=<ID=PCATxALT,Description="Alternate alleles do not match PharmCAT alternate alleles">
##FILTER=<ID=PCATxINDEL,Description="Unexpected format for INDELs">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample_1
 ```

# instance 1: chromosome names

```vcf deficient format
1	97078987	rs114096998	G	T	.	PASS	.	GT	0/0
```

```vcf accurate format
chr1	97078987	rs114096998	G	T	.	PASS	.	GT	0/0
```

# instance 2: homozygous reference SNPs with unspecified ALT

```vcf deficient format
chr1	97078993	rs148799944	C	.	.	PASS	.	GT	0/0
chr1	97079005	rs140114515	C	<*>	.	PASS	.	GT	0/0
```

```vcf accurate format
chr1	97078993	rs148799944	C	G	.	PASS	.	GT	0/0
chr1	97079005	rs140114515	C	T	.	PASS	.	GT	0/0
```

# instance 3: left-alignment for a tandem repeat INDELs

Note that the ALT allele is updated after left-alignment.

```vcf deficient format
chr1	97740414	rs72549309	AATGA	A	.	PASS	.	GT	1/0
```

```vcf accurate format
chr1	97740410	rs72549309	GATGA	G	.	PASS	.	GT	1/0
```

# instance 3: multiallelic INDELs

For a genomic position that harbors multiple pharmacogenetic INDELs, all pharmacogenetic INDEL variants should be combined into a multiallelic record. 

```vcf deficient format
chr2	233760233	rs3064744	C	CAT	.	PASS	.	GT	1/0
chr2	233760233	rs3064744	CAT	C	.	PASS	.	GT	0/0
chr2	233760233	rs3064744	C	CATAT	.	PASS	.	GT	0/1
```

```vcf accurate format
chr2	233760233	rs3064744	CAT	CATAT,CATATAT,C	.	PASS	.	GT	1/2
```

# instance 4: prioritizing pharmacogenetic variant(s) over non-pharmacogenetic ones at a multiallelic locus

Note that the pharmacogenetic SNP is reordered before the INDEL that is located at the same genomic position.

```vcf deficient format
chr7	117509035	.	GA	G	.	PASS	.	GT	0/0
chr7	117509035	.	G	A	.	PASS	.	GT	0/1
```

```vcf accurate format
chr7	117509035	rs397508256	G	A	.	PASS	.	GT	0/1
chr7	117509035	.	GA	G	.	PCATxREF	.	GT	0/0
```

# instance 5: left-alignment for long INDELs

Note that the ALT allele is updated after left-alignment.

```vcf deficient format
chr10	94942212	rs1304490498	AAGAAATGGAA	A	.	PASS	.	GT	1/0
```

```vcf accurate format
chr10	94942205	rs1304490498	CAATGGAAAGA	C	.	PASS	.	GT	1/0
```

# instance 6: warning for unexpected INDEL format

In this instance, chr10:94949281 is homozygous reference. There is no sufficient information to infer genotypes for a pharmacogenetic INDEL `chr10:94949281:GA:G` at this position. A warning should be listed. And the pharmacogenetic variant should be reported as missing.

```vcf deficient format
chr10	94949281	.	G	.	.	PASS	.	GT	0/0
```

```vcf accurate format
chr10	94949281	.	G	.	.	PCATxINDEL	.	GT	0/0
```

# instance 7: prioritizing pharmacogenetic variants at positions that have both SNP(s) and INDEL(s)
!! need to double check

Note that pharmacogenetic INDELs should be listed ahead of the non-pharmacogenetic SNPs. Information of other alternative alleles need to be completed so that the variant representation format matches the format that PharmCAT expects.

```vcf deficient format
chr13	48037782	.	A	C,<*>	.	PASS	.	GT	0/0
chr13	48037782	rs746071566	AGGAGTC	A,<*>	.	PASS	.	GT	0/0
```

```vcf accurate format
chr13	48037782	rs746071566	AGGAGTC	AGGAGTCGGAGTC,A	.	PASS	.	GT	0/0
chr13	48037782	rs746071566	A	C	.	PASS	.	GT	0/0
chr13	48037782	rs746071566	AGGAGTC	<*>	.	PCATxINDEL	.	GT	0/0
```

# instance 8: filling up multiallelic SNPs

If a SNP position is present in the VCF while the exact pharmacogenetic allele is not present, there is sufficient information to infer that the specific pharmacogenetic alleles are not present in the VCF. It should be safe to add back the pharmacogenetic alleles back to the VCF. Bcftools should probably adjust the genotypes after the alleles are added back.

```vcf deficient format
chr19	38448712	rs121918592	G	A,<*>	.	PASS	.	GT	1/0
chr19	40991381	rs33973337	A	T	.	PASS	.	GT	1/0
```

```vcf accurate format
chr19	38448712	rs121918592	G	A,C	.	PASS	.	GT	1/0
chr19	40991381	rs33973337	A	C,T	.	PASS	.	GT	2/0
```

# instance 9: prioritizing pharmacogenetic variants over non-pharmacogenetic variants 

This is similar to instance 4 but observed in actual UK Biobank dataset

```vcf deficient format
chrX	154532608	.	C	CG,<*>	0	RefCall	.	0/0
chrX	154532608	.	C	T,<*>	0	PASS	.	0/0
chrX	154532608	.	CG	C,<*>	0	RefCall	.	0/0
```

```vcf accurate format
chrX	154532608	.	C	T	0	RefCall	.	0/0
chrX	154532608	.	C	CG	0	RefCall	.	0/0
chrX	154532608	.	CG	<*>	0	PCATxREF	.	0/0
```
