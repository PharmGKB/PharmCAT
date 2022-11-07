#!/bin/bash
#
# Builds reference FASTA file for use by VCF preprocessor.
# Download from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/
# * The no_alt_analysis_set contains the sequences, in FASTA format, of the chromosomes, mitochondrial genome,
#   un-localized scaffolds, and unplaced scaffolds.
# * The definition line has UCSC-style sequence identifier and contains metadata in a series of space-separated
#   tag-value pairs.
# * Soft-masking (low complexity sequence aka repetitive regions) is converted to uppercase in reference genome
#   sequences for alignment pipelines.
#
set -e
set -u
set -o pipefail

echo "Removing old files..."
rm -f GRCh38_reference_fasta.tar reference.fna.bgz reference.fna.bgz.fai reference.fna.bgz.gzi genomic.fna genomic.fna.gz chrfix.fna

echo ""
echo "Downloading from NCBI..."
curl -#fSL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz -o genomic.fna.gz

echo ""
echo "Uncompressing..."
gunzip genomic.fna.gz

# 'AS' means assembly-name. Only keep the main GRCh38 assembly sequence
echo ""
echo "Removing non-primary sequences"
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' genomic.fna | grep -v "^>chr\S*_" - | tr "\t" "\n" > genomic.short.fna

echo ""
echo "Compressing with bgzip..."
bgzip -c genomic.short.fna > reference.fna.bgz

echo ""
echo "Building index..."
samtools faidx reference.fna.bgz

echo ""
echo "Creating tarball..."
tar -czvf GRCh38_reference_fasta.tar reference.fna.bgz reference.fna.bgz.fai reference.fna.bgz.gzi

echo ""
echo "Done."
