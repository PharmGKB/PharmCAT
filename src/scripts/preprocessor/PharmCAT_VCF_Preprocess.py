#!/usr/bin/env python3

__author__ = 'BinglanLi'

import os
import subprocess
import sys
from pathlib import Path
from timeit import default_timer as timer

import vcf_preprocess_utilities as util


def run(args):
    ## normalize and prepare the input VCF for PharmCAT
    start = timer()

    # validate input
    if bool(args.input_vcf) == bool(args.input_list):
        print("Invalid input file(s). Provide either [--input_vcf] or [--input_list].")
        sys.exit(1)

    # validate bcftools
    bcftools_executable_path = args.path_to_bcftools if args.path_to_bcftools else 'bcftools'
    try:
        subprocess.run([bcftools_executable_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except:
        print('Error: %s not found or not executable' % bcftools_executable_path)
        sys.exit(1)
    # validate tabix
    tabix_executable_path = args.path_to_tabix if args.path_to_tabix else 'tabix'
    try:
        subprocess.run([tabix_executable_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except:
        print('Error: %s not found or not executable' % tabix_executable_path)
        sys.exit(1)


    # list of files to be deleted
    tmp_files_to_be_removed = []
    # define working directory, default the directory of the first input VCF
    if args.output_folder:
        output_dir = args.output_folder
    elif args.input_vcf:
        output_dir = os.path.split(args.input_vcf)[0]
    else:
        with open(args.input_list, 'r') as file:
            for line in file:
                if os.path.isfile(line):
                    output_dir = os.path.split(line)[0]
                    break
        file.close()
    # create the output folder
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # download the human reference sequence if not provided
    path_to_ref_seq = util.download_grch38_ref_fasta_and_index(output_dir) if not args.ref_seq else args.ref_seq

    # index ref_pgx_vcf if not indexed
    if not os.path.exists(args.ref_pgx_vcf + '.tbi'):
        util.tabix_index_vcf(tabix_executable_path, args.ref_pgx_vcf)

    # read the sample list
    sample_list = []
    if args.sample_file:
        with open(args.sample_file, 'r') as file:
            for line in file:
                line = line.strip()
                sample_list.append(line)
        file.close()
    elif args.input_vcf:
        sample_list = util.obtain_vcf_sample_list(bcftools_executable_path, args.input_vcf)
    else:
        with open(args.input_list, 'r') as file:
            for line in file:
                if os.path.isfile(line):
                    sample_list = util.obtain_vcf_sample_list(bcftools_executable_path, line)
                    break
        file.close()

    # shrink input VCF down to PGx allele defining regions and selected samples
    # modify input VCF chromosomes naming format to <chr##>
    if args.input_list:
        vcf_pgx_regions = util.extract_regions_from_multiple_files(bcftools_executable_path, tabix_executable_path,
                                                                                args.input_list, path_to_ref_seq, args.ref_pgx_vcf, output_dir, args.output_prefix, sample_list)
    else:
        vcf_pgx_regions = util.extract_regions_from_single_file(bcftools_executable_path, tabix_executable_path,
                                                                             args.input_vcf, args.ref_pgx_vcf, output_dir, args.output_prefix, sample_list)
    tmp_files_to_be_removed.append(vcf_pgx_regions)

    # normalize the input VCF
    vcf_normalized = util.normalize_vcf(bcftools_executable_path, tabix_executable_path, vcf_pgx_regions, path_to_ref_seq)
    tmp_files_to_be_removed.append(vcf_normalized)

    # extract the specific PGx genetic variants in the reference PGx VCF
    # this step also generates a report of missing PGx positions in the input VCF
    vcf_normalized_pgxOnly = util.filter_pgx_variants(bcftools_executable_path, tabix_executable_path, vcf_normalized, path_to_ref_seq, args.ref_pgx_vcf, output_dir, args.output_prefix)
    tmp_files_to_be_removed.append(vcf_normalized_pgxOnly)

    # output PharmCAT-ready single-sample VCF
    # retain only the PharmCAT allele defining positions in the output VCF file
    util.output_pharmcat_ready_vcf(bcftools_executable_path, vcf_normalized_pgxOnly, output_dir, args.output_prefix, sample_list)

    # remove intermediate files
    if not args.keep_intermediate_files:
        for single_path in tmp_files_to_be_removed:
           util.remove_vcf_and_index(single_path)

    end = timer()
    print("Successfully preprocessed input VCF in %s seconds"%(str(end-start)))

if __name__ == "__main__":
    import argparse

    # describe the 
    parser = argparse.ArgumentParser(description='Prepare an input VCF for the PharmCAT')

    # list arguments
    parser.add_argument("--input_vcf", type = str, help="Load a compressed VCF file.")
    parser.add_argument("--input_list", type = str, help="A sorted list of by-chromosome, compressed VCF file names.")
    parser.add_argument("--ref_seq", help="Load the Human Reference Genome GRCh38/hg38 in the fasta format.")
    parser.add_argument("--ref_pgx_vcf", required=True, type = str, help="Load a VCF file of PGx variants. This file is available from the PharmCAT GitHub release.")
    parser.add_argument("--sample_file", help="A file of samples to be prepared for the PharmCAT, one sample at a line.")
    parser.add_argument("--path_to_bcftools", help="Load an alternative path to the executable bcftools.")
    parser.add_argument("--path_to_tabix", help="Load an alternative path to the executable tabix.")
    parser.add_argument("--output_folder", type = str, help="Directory of the output VCF, by default, the directory of the first input VCF.")
    parser.add_argument("--output_prefix", default = 'pharmcat_ready_vcf', type = str, help="Prefix of the output VCF")
    parser.add_argument("--keep_intermediate_files", action='store_true', help="Keep intermediate files, false by default.")

    # parse arguments
    args = parser.parse_args()

    # to be removed at the PharmCAT v1.0
    print("Warning: For PharmCAT v1.0.0, please use the \'pharmcat_positions_1.0.0_sorted.vcf.gz\' under "
          "\'src/scripts/test/pharmcat_pgx_vcf/\'")

    # normalize variant representations and reconstruct multi-allelic variants in the input VCF
    run(args)
