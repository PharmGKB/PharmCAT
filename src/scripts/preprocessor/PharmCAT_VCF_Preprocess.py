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
    # validate input
    if bool(args.input_vcf) == bool(args.input_list):
        print("Invalid input file(s). Provide either [--input_vcf] or [--input_list].")
        sys.exit(1)

    # validate bcftools
    bcftools_path = args.path_to_bcftools if args.path_to_bcftools else 'bcftools'
    try:
        subprocess.run([bcftools_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except:
        print('Error: %s not found or not executable' % bcftools_path)
        sys.exit(1)
    # validate tabix
    tabix_path = args.path_to_tabix if args.path_to_tabix else 'tabix'
    try:
        subprocess.run([tabix_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except:
        print('Error: %s not found or not executable' % tabix_path)
        sys.exit(1)

    """
    organize arguments
    """
    input_list = args.input_list
    input_vcf = args.input_vcf
    pgx_vcf = args.ref_pgx_vcf
    sample_file = args.sample_file
    output_prefix = args.output_prefix
    keep_intermediate_files = args.keep_intermediate_files
    # define working directory, default the directory of the first input VCF
    if args.output_folder:
        output_dir = args.output_folder
    elif input_vcf:
        output_dir = os.path.split(input_vcf)[0]
    else:
        with open(input_list, 'r') as file:
            for line in file:
                if os.path.isfile(line):
                    output_dir = os.path.split(line)[0]
                    break
        file.close()
    # create the output folder
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    # download the human reference sequence if not provided
    ref_seq = util.download_grch38_ref_fasta_and_index(output_dir) if not args.ref_seq else args.ref_seq

    # index ref_pgx_vcf if not indexed
    if not os.path.exists(pgx_vcf + '.tbi'):
        util.tabix_index_vcf(tabix_path, pgx_vcf)

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
                if os.path.isfile(line):
                    sample_list = util.obtain_vcf_sample_list(bcftools_path, line)
                    break
        file.close()

    # list of files to be deleted
    tmp_files_to_be_removed = []

    """
    normalize and prepare vcf for PharmCAT
    """
    # shrink input VCF down to PGx allele defining regions and selected samples
    # modify input VCF chromosomes naming format to <chr##>
    if input_list:
        vcf_pgx_regions = util.extract_regions_from_multiple_files(bcftools_path, tabix_path, input_list, ref_seq,
                                                                   pgx_vcf, output_dir, output_prefix, sample_list)
    else:
        vcf_pgx_regions = util.extract_regions_from_single_file(bcftools_path, tabix_path, input_vcf,
                                                                pgx_vcf, output_dir, output_prefix, sample_list)
    tmp_files_to_be_removed.append(vcf_pgx_regions)

    # normalize the input VCF
    vcf_normalized = util.normalize_vcf(bcftools_path, tabix_path, vcf_pgx_regions, ref_seq)
    tmp_files_to_be_removed.append(vcf_normalized)

    # extract the specific PGx genetic variants in the reference PGx VCF
    # this step also generates a report of missing PGx positions in the input VCF
    vcf_normalized_pgx_only = util.filter_pgx_variants(bcftools_path, tabix_path, vcf_normalized, ref_seq, pgx_vcf,
                                                       output_dir, output_prefix)
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
    parser.add_argument("--input_vcf", type=str, help="Load a compressed VCF file.")
    parser.add_argument("--input_list", type=str,
                        help="A sorted list of by-chromosome, compressed VCF file names.")
    parser.add_argument("--ref_pgx_vcf", required=True, type=str,
                        help="(Optional) a sorted VCF of PGx variants.")
    parser.add_argument("--ref_seq",
                        help="(Optional) the Human Reference Genome GRCh38/hg38 in the fasta format.")
    parser.add_argument("--sample_file",
                        help="(Optional) a file of samples to be prepared for the PharmCAT, one sample at a line.")
    parser.add_argument("--path_to_bcftools",
                        help="(Optional) an alternative path to the executable bcftools.")
    parser.add_argument("--path_to_tabix",
                        help="(Optional) an alternative path to the executable tabix.")
    parser.add_argument("--output_folder", type=str,
                        help="(Optional) directory for outputs, by default, directory of the first input VCF.")
    parser.add_argument("--output_prefix", default='pharmcat_ready_vcf', type=str,
                        help="(Optional) prefix of the output VCF, default = \'pharmcat_ready_vcf\'")
    parser.add_argument("--keep_intermediate_files", action='store_true',
                        help="(Optional) keep intermediate files, false by default.")

    # parse arguments
    args = parser.parse_args()

    # normalize variant representations and reconstruct multi-allelic variants in the input VCF
    run(args)
