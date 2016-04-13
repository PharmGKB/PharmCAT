# This is a Python script for dealing with variant position re-mapping
# between human genome assemblies
#
# Author: Chunlei Wu

import re

chr_mapping_d = {
    'NC_000001': 'chr1',
    'NC_000002': 'chr2',
    'NC_000006': 'chr6',
    'NC_000010': 'chr10',
    'NC_000012': 'chr12'
}

def get_variant_list(input_tsv_files):
    variant_li = []
    for input_tsv in input_tsv_files:
        with open(input_tsv) as input_f:
            gene_line = input_f.readline()
            gene = gene_line.split('\t')[0].split(':')[1].strip()
            #skip next two lines
            input_f.readline()
            input_f.readline()
            #get the position line:
            pos_line = input_f.readline().strip().split('\t')
            match = re.search(r'(NC_\d+\.\d+)', pos_line[0])
            if match:
                chr_seq = match.group(0)
            else:
                chr_seq = None
            # make sure all position is on "GRCh38.p2"
            assert pos_line[0].find('GRCh38.p2') != -1
            for _var in pos_line[1:]:
                for _v in _var.split(';'):
                    # deal with multiallelic variants
                    variant_li.append((gene, chr_seq, _v.strip()))
    return variant_li
