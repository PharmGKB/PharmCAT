# This is a Python script for dealing with variant position re-mapping
# between human genome assemblies
#
# Author: Chunlei Wu
from __future__ import print_function
import re


chr_mapping_d = {
    'NC_000001': 'chr1',
    'NC_000002': 'chr2',
    'NC_000003': 'chr3',
    'NC_000004': 'chr4',
    'NC_000005': 'chr5',
    'NC_000006': 'chr6',
    'NC_000007': 'chr7',
    'NC_000008': 'chr8',
    'NC_000009': 'chr9',
    'NC_000010': 'chr10',
    'NC_000011': 'chr11',
    'NC_000012': 'chr12',
    'NC_000013': 'chr13',
    'NC_000014': 'chr14',
    'NC_000015': 'chr15',
    'NC_000016': 'chr16',
    'NC_000017': 'chr17',
    'NC_000018': 'chr18',
    'NC_000019': 'chr19',
    'NC_000020': 'chr20',
    'NC_000021': 'chr21',
    'NC_000022': 'chr22',
    'NC_000023': 'chrX',
    'NC_000024': 'chrY'
}

def get_variant_list(input_tsv_files, as_ucsc_region=True):
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
    variant_li = [_get_variant_region(var, as_ucsc_region=as_ucsc_region) for var in variant_li]
    return variant_li


def _get_variant_region(var, as_ucsc_region=False):
    '''A helper function, input is a tuple from above get_variant_list function.
       Note that the returned ucsc region has the start position consistent with VCF style,
       that is, for del/ins, start is 1-based ahead.
    '''

    assert isinstance(var, (tuple, list)) and len(var) == 3, "Should be a tuple/list of three items"
    chr = chr_mapping_d.get(var[1].split('.')[0], None)
    if chr is None:
        print("unknown NC seq accessioon: {}".format(var[1]))
    hgvs = var[2]
    # check if it is a snp
    mat = re.match('g\.(\d+)\w{1}\>\w{1}', hgvs)
    if mat:
        start, end = mat.group(1), mat.group(1)
    else:
        # check if it is a del or ins
        mat = re.match('g\.([\d_]+)(del|ins)', hgvs)
        if mat:
            pos = mat.group(1)
            if pos.find('_') != -1:
                # so it's a range
                start, end = pos.split('_')
            else:
                # it is a position
                start, end = pos, pos
            # move start pos 1-base ahead to match pos in VCF file
            start = str(int(start) - 1)
        else:
            # skipping other types for now
            print("\"{}\" is skipped".format(var))
            start, end = None, None
    if as_ucsc_region:
        return var + ("{}:{}-{}".format(chr, start, end), )
    return var + (chr, start, end)

# remap_api cmd line for converting grch38.p2 to grch37.p13 positions
# ./remap_api.pl --mode asm-asm --from GCF_000001405.28 --dest GCF_000001405.25 --annotation cpic_input.txt --annot_out cpic_output.txt --report_out report.txt --in_format region
#
