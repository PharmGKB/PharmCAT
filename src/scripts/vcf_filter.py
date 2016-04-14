# This is a Python script for dealing with variant position re-mapping
# between human genome assemblies
#
# Author: Chunlei Wu
#
# This script requires pyVCF module, which can be installed easily as:
#     pip install pyVCF
#

from __future__ import print_function
import sys

import vcf.filters


CPIC_POS_MAPPING_FILE = './cpic_pos_mapping.txt'
SUPPORTED_ASSEMBLIES = ['hg19', 'hg38']


def get_pos_mapping():
    hg38_pos_li = []
    hg19_pos_li = []
    with open(CPIC_POS_MAPPING_FILE) as pos_f:
        for line in pos_f:
            ld = line.strip().split('\t')
            if len(ld) < 5:
                continue
            chr = ld[3].split(':')[0]
            start_hg38 = int(ld[3].split(':')[1].split('-')[0])
            start_hg19 = int(ld[4].split(':')[1].split('-')[0])
            hg38_pos_li.append((chr, start_hg38))
            hg19_pos_li.append((chr, start_hg19))
    return {'hg38': hg38_pos_li, 'hg19': hg19_pos_li}


def get_cpic_pos_list(assembly='hg38'):
    assert assembly in ['hg38', 'hg19'], "Unsupported assembly."
    cpic_pos_li = []
    with open(CPIC_POS_MAPPING_FILE) as pos_f:
        for line in pos_f:
            ld = line.split('\t')
            if assembly == 'hg38':
                start = ld[3].split(':')[1].split('-')[0]
            elif assembly == 'hg19':
                start = ld[4].split(':')[1].split('-')[0]
            cpic_pos_li.append(start)

    return cpic_pos_li


def vcf_filter(input_vcf, output_vcf=None):
    '''input_vcf is the path to the input vcf file (gzipped or not)'''

    pos_mapping_d = get_pos_mapping()
    if output_vcf:
        output_f = open(output_vcf, 'w')
    else:
        output = []

    cnt = 0
    with open(input_vcf) as vcf_f:
        reader = vcf.Reader(vcf_f)
        assembly = reader.metadata['reference']
        assert assembly in SUPPORTED_ASSEMBLIES, "Assembly \"{}\" not supported!".format(assembly)
        cpic_pos_li = pos_mapping_d[assembly]
        if output_vcf:
            output = vcf.Writer(output_f, reader)
        for record in reader:
            if (record.CHROM, record.POS) in cpic_pos_li:
                cnt += 1
                if output_vcf:
                    output.write_record(record)
                else:
                    output.append(record)
        print("Found {} qualified VCF lines".format(cnt))
        if output_vcf:
            output.close()
        else:
            return output


def main():
    if len(sys.argv) != 3:
        print("Usage: python vcf_filter.py <input_vcf_file> <output_vcf_file>")
    else:
        input_vcf, output_vcf = sys.argv[1], sys.argv[2]
        vcf_filter(input_vcf, output_vcf)


if __name__ == '__main__':
    main()
