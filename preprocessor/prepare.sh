#!/bin/bash

file=$1

echo "Starting ${file}"
# use bcftools to split a multipart vcf if present
# https://samtools.github.io/bcftools/bcftools.html
# Suggest downloading brew or bioconda
for sample in `bcftools view -h $file | grep "^#CHROM" | cut -f10-`; do
    dir="$(dirname ${file})"
    echo "START"
    echo "$sample"
    bcftools view -I -s $sample $file > $sample.vcf

    echo "Remove non variants"
    awk '$0 ~ "^#" || ($5 != "." && $5 != "<*>")' $sample.vcf > $sample.filtered.vcf

    # add chr to vcf file (required by PharmCAT) 
    echo "fixing chr"
    perl -pe '/^((?!^chr).)*$/ && s/^([^#])/chr$1/gsi' $sample.filtered.vcf > $sample.chr.vcf

    # bzip and tabix the data - required by vcf-merge
    echo "bzip and tabix results"
    vcf-sort $sample.chr.vcf > $sample.chr.sorted.vcf
    bgzip -c $sample.chr.sorted.vcf > $sample.chr.vcf.gz
    tabix -p vcf $sample.chr.vcf.gz

    # combine with PharmCAT vcf
    echo "Combine with PharmCAT vcf"
    vcf-merge $sample.chr.vcf.gz pgx.vcf.gz > $sample.chr.combined.vcf

    # remove variants not in PharmCAT
    echo "Remove no PGX variants"
    grep -i "PX\|#" $sample.chr.combined.vcf > $sample.chr.combined_removed.vcf  

    # replace empty positions
    echo "Replace unknown positions with 0/0 (not ideal - better to extract from .bam or a .bed file)"
    awk '{gsub("\t.\t0/0","\t0/0\t0/0");print $0}' $sample.chr.combined_removed.vcf > $sample.chr.combined_removed2.vcf

    # removed second column of sample data (from the PharmCAT vcf)
    cut -f1-10 $sample.chr.combined_removed2.vcf > ${file/.vcf*/.$sample.final.vcf}

    # Iniital run of PharmCAT as a sanity check.  Can be re-run with more parameters etc.
    echo "Run PharmCAT initial check"
    java -jar /PharmCAT/build/libs/pharmcat-*.jar  -o $dir -vcf ${file/.vcf*/.$sample.final.vcf}  -k -j


    # parse the output json to produce a summary line
    echo "Parse json"
    python parse.py ${file/.vcf*/.$sample.final.call.json}
done