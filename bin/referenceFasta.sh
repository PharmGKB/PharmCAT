#!/bin/bash
#
# Builds reference FASTA file for use by VCF preprocessor.
#

echo "Removing old files..."
rm -f GRCh38_reference_fasta.tar reference.fna.bgz reference.fna.bgz.fai reference.fna.bgz.gzi genomic.fna genomic.fna.gz chrfix.fna

echo ""
echo "Downloading from NCBI..."
curl -#fSL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz -o genomic.fna.gz

echo ""
echo "Uncompressing..."
gunzip genomic.fna.gz

echo ""
echo "Fixing chromosomes..."
cat genomic.fna | sed -r 's/^>(NC.*Homo sapiens chromosome ([0-9XY]+),.*)/>chr\2 \1/g' | sed -r 's/^>(NC.*Homo sapiens mitochondrion,.*)/>chrM \1/g' > chrfix.fna

echo ""
echo "Removing non-primary sequences"
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' chrfix.fna | grep ">chr" - | tr "\t" "\n" > chrfix.short.fna

echo ""
echo "Compressing with bgzip..."
bgzip -c chrfix.short.fna > reference.fna.bgz

echo ""
echo "Building index..."
samtools faidx reference.fna.bgz

echo ""
echo "Creating tarball..."
tar -czvf GRCh38_reference_fasta.tar reference.fna.bgz reference.fna.bgz.fai reference.fna.bgz.gzi

echo ""
echo "Done."
