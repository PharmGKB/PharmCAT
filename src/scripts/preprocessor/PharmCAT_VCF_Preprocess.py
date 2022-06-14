#!/usr/bin/env python3

__author__ = 'BinglanLi'

import os
import subprocess
import sys
from pathlib import Path
from timeit import default_timer as timer

import vcf_preprocess_utilities as util


def run(args):
    ''' normalize and prepare the input VCF for PharmCAT'''
    start = timer()

    """
    validate arguments
    """
    # validate bcftools
    bcftools_path = args.path_to_bcftools if args.path_to_bcftools else 'bcftools'
    try:
        subprocess.run([bcftools_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except:
        print('Error: %s not found or not executable' % bcftools_path)
        sys.exit(1)
    # validate bgzip
    bgzip_path = args.path_to_bgzip if args.path_to_bgzip else 'bgzip'
    try:
        subprocess.run([bgzip_path, '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except:
        print('Error: %s not found or not executable' % bgzip_path)
        sys.exit(1)

    # validate input
    input_list = args.input_list
    input_vcf = args.input_vcf
    # check input arguments. can only take one of the two input options
    if bool(input_vcf) == bool(input_list):
        print("Invalid input file(s). Provide either [--input_vcf] or [--input_list].")
        sys.exit(1)
    # if single input vcf, validate and bgzip
    if input_vcf:
        if not os.path.exists(input_vcf):
            print("Cannot find", input_vcf)
            sys.exit(1)
        # compress if the input vcf is not bgzipped
        input_vcf = util.bgzipped_vcf(bgzip_path, input_vcf)
    # if a list of vcfs, validate the list file; and bgzip later
    if input_list:
        if not os.path.exists(input_list):
            print("Cannot find", input_list)
            sys.exit(1)

    # validate the reference vcf of PharmCAT PGx positions
    ref_pgx = args.ref_pgx_vcf
    if not os.path.exists(ref_pgx):
        print('Error: VCF of the reference PGx positions was not found at: %s' % ref_pgx)
        sys.exit(1)

    """
    organize the rest of the arguments
    """
    sample_file = args.sample_file
    output_prefix = args.output_prefix
    keep_intermediate_files = args.keep_intermediate_files
    missing_to_ref = args.missing_to_ref
    # define working directory, default the directory of the first input VCF
    if args.output_folder:
        output_dir = args.output_folder
        # create the output folder
        Path(output_dir).mkdir(parents=True, exist_ok=True)
    elif input_vcf:
        output_dir = os.path.split(input_vcf)[0]
    else:
        output_dir = os.path.split(input_list)[0]
    print("Saving output to", output_dir)

    # download the human reference sequence if not provided
    if args.ref_seq:
        ref_seq = args.ref_seq
    else:
        if os.path.exists(os.path.join(output_dir, 'reference.fna.bgz')):
            ref_seq = os.path.join(output_dir, 'reference.fna.bgz')
            print("Using default FASTA reference at ", ref_seq)
        elif os.path.exists(os.path.join(os.getcwd(), 'reference.fna.bgz')):
            ref_seq = os.path.join(os.getcwd(), 'reference.fna.bgz')
            print("Using default FASTA reference at ", ref_seq)
        else:
            ref_seq = util.get_default_grch38_ref_fasta_and_index(output_dir)
            print("Downloaded to %s" % ref_seq)

    # index ref_pgx if not already so
    if not os.path.exists(ref_pgx + '.csi'):
        util.index_vcf(bcftools_path, ref_pgx)

    # read the sample list
    sample_list = []
    if sample_file:
        with open(sample_file, 'r') as file:
            for line in file:
                line = line.strip()
                sample_list.append(line)
        file.close()
    elif input_vcf:
        sample_list = util.obtain_vcf_sample_list(bcftools_path, input_vcf)
    else:
        with open(input_list, 'r') as file:
            for line in file:
                line = line.strip()
                if os.path.isfile(line):
                    sample_list = util.obtain_vcf_sample_list(bcftools_path, line)
                    break
        file.close()

    # check if sample name violates bcftools sample name convention
    if any(',' in sample_name for sample_name in sample_list):
        print('Please remove comma \',\' from sample names, which violates bcftools sample name convention')
        sys.exit(1)

    # list of files to be deleted
    tmp_files_to_be_removed = []

    """
    normalize and prepare vcf for PharmCAT
    """
    # shrink input VCF down to PGx allele defining regions and selected samples
    # modify input VCF chromosomes naming format to <chr##>
    if input_list:
        vcf_pgx_regions = util.extract_regions_from_multiple_files(bcftools_path, bgzip_path, input_list,
                                                                   ref_pgx, output_dir, output_prefix, sample_list)
    else:
        vcf_pgx_regions = util.extract_regions_from_single_file(bcftools_path, input_vcf,
                                                                ref_pgx, output_dir, output_prefix, sample_list)
    tmp_files_to_be_removed.append(vcf_pgx_regions)

    # normalize the input VCF
    vcf_normalized = util.normalize_vcf(bcftools_path, vcf_pgx_regions, ref_seq)
    tmp_files_to_be_removed.append(vcf_normalized)

    # extract the specific PGx genetic variants in the reference PGx VCF
    # this step also generates a report of missing PGx positions in the input VCF
    vcf_normalized_pgx_only = util.filter_pgx_variants(bcftools_path, bgzip_path, vcf_normalized, ref_seq,
                                                       ref_pgx, missing_to_ref, output_dir, output_prefix)
    tmp_files_to_be_removed.append(vcf_normalized_pgx_only)

    # output PharmCAT-ready single-sample VCF
    # retain only the PharmCAT allele defining positions in the output VCF file
    util.output_pharmcat_ready_vcf(bcftools_path, vcf_normalized_pgx_only, output_dir, output_prefix, sample_list)

    # remove intermediate files
    if not keep_intermediate_files:
        for single_path in tmp_files_to_be_removed:
            util.remove_vcf_and_index(single_path)

    end = timer()
    print("Successfully preprocessed input VCF in %s seconds" % (str(end - start)))


if __name__ == "__main__":
    import argparse

    # describe the tool
    parser = argparse.ArgumentParser(description='Prepare an input VCF for the PharmCAT')

    # list arguments
    parser.add_argument("-i", "--input_vcf", type=str, help="Path to a VCF file.")
    parser.add_argument("--input_list", type=str,
                        help="A file containing paths to VCF files (one file per line), sorted by chromosome position.")
    parser.add_argument("-p", "--ref_pgx_vcf", type=str,
                        default=os.path.join(os.getcwd(), "pharmcat_positions.vcf.bgz"),
                        help="A sorted VCF of PharmCAT PGx variants, gzipped with preprocessor scripts. Default = "
                             "\'pharmcat_positions.vcf.bgz\' in the current working directory.")
    parser.add_argument("-g", "--ref_seq",
                        help="(Optional) the Human Reference Genome GRCh38/hg38 in the fasta format.")
    parser.add_argument("-S", "--sample_file",
                        help="(Optional) a file of samples to be prepared for the PharmCAT, one sample at a line.")
    parser.add_argument("--path_to_bcftools",
                        help="(Optional) an alternative path to the executable bcftools.")
    parser.add_argument("--path_to_bgzip",
                        help="(Optional) an alternative path to the executable bgzip.")
    parser.add_argument("--output_folder", type=str,
                        help="(Optional) directory for outputs, by default, directory of the first input VCF.")
    parser.add_argument("--output_prefix", default='pharmcat_ready_vcf', type=str,
                        help="(Optional) prefix of the output VCF, default = \'pharmcat_ready_vcf\'.")
    parser.add_argument("--keep_intermediate_files", action='store_true',
                        help="(Optional) keep intermediate files, false by default.")
    parser.add_argument("-0", "--missing_to_ref", action='store_true',
                        help="(Optional) assume genotypes at missing PGx sites are 0/0.  DANGEROUS!.")

    # parse arguments
    args = parser.parse_args()

    # print warnings here
    # alternatively, could use the "warnings" module
    if args.missing_to_ref:
        print('=============================================================\n'
              'Warning: Argument "-0"/"--missing_to_ref" supplied\n'
              '\n'
              'THIS SHOULD ONLY BE USED IF: you sure your data is reference\n'
              'at the missing positions instead of unreadable/uncallable at\n'
              'those positions.\n'
              '\n'
              'Running PharmCAT with positions as missing vs reference can\n'
              'lead to different results.\n'
              '=============================================================\n')

    # normalize variant representations and reconstruct multi-allelic variants in the input VCF
    run(args)
