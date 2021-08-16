#!/usr/bin/env python3

__author__ = 'BinglanLi'

import os
import sys
import subprocess
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

    # organize args
    current_working_dir = os.getcwd()

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

    # create the output folder if not existing
    Path(args.output_folder).mkdir(parents=True, exist_ok=True)

    # list of files to be deleted
    tmp_files_to_be_removed = []

    # download the human reference sequence if not provided
    path_to_ref_seq = util.download_grch38_ref_fasta_and_index(current_working_dir) if not args.ref_seq else args.ref_seq

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
        first_valid_file = False
        with open(args.input_list, 'r') as file:
            for line in file:
                if not first_valid_file:
                    line = line.strip()
                    if os.path.isfile(line):
                        sample_list = util.obtain_vcf_sample_list(bcftools_executable_path, line)
                        first_valid_file = True
                else:
                    break
        file.close()

    # shrink input VCF down to PGx allele defining regions and selected samples
    # modify input VCF chromosomes naming format to <chr##>
    if args.input_list:
        intermediate_vcf_pgx_regions = util.extract_pharmcat_regions_from_multiple_files(bcftools_executable_path, tabix_executable_path,
                                                                                         args.input_list, path_to_ref_seq, args.ref_pgx_vcf, args.output_folder, args.output_prefix, sample_list)
    else:
        intermediate_vcf_pgx_regions = util.extract_pharmcat_regions_from_single_file(bcftools_executable_path, tabix_executable_path,
                                                                         args.input_vcf, args.ref_pgx_vcf, args.output_folder, args.output_prefix, sample_list)
    tmp_files_to_be_removed.append(intermediate_vcf_pgx_regions)

    # merge the input VCF with the PGx position file provided by '--ref_pgx_vcf'
    # run this step to ensure the output VCF will have THE SAME VARIANT REPRESENTATION as PharmCAT does
    intermediate_vcf_pgx_merged = util.merge_vcfs(bcftools_executable_path, tabix_executable_path, intermediate_vcf_pgx_regions, args.ref_pgx_vcf)
    tmp_files_to_be_removed.append(intermediate_vcf_pgx_merged)

    # normalize the input VCF
    intermediate_vcf_normalized = util.normalize_vcf(bcftools_executable_path, tabix_executable_path, intermediate_vcf_pgx_merged, path_to_ref_seq, args.ref_pgx_vcf)
    tmp_files_to_be_removed.append(intermediate_vcf_normalized)

    # generate a report of missing PGx positions in VCF format
    util.output_missing_pgx_positions(bcftools_executable_path, intermediate_vcf_normalized, args.ref_pgx_vcf, args.output_folder, args.output_prefix)

    # output PharmCAT-ready single-sample VCF
    # retain only the PharmCAT allele defining positions in the output VCF file
    util.output_pharmcat_ready_vcf(bcftools_executable_path, intermediate_vcf_normalized, args.output_folder, args.output_prefix, sample_list)

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
    parser.add_argument("--output_folder", default = os.getcwd(), type = str, help="Directory of the output VCF, by default, current working directory.")
    parser.add_argument("--output_prefix", default = 'pharmcat_ready_vcf', type = str, help="Prefix of the output VCF")
    parser.add_argument("--keep_intermediate_files", action='store_true', help="Keep intermediate files, false by default.")

    # parse arguments
    args = parser.parse_args()

    # to be removed at the PharmCAT v1.0
    print("Warning: For PharmCAT v0.8.0, please use the \'pharmcat_positions_0.8.0_updated_06222021.vcf.gz\' under "
          "\'src/scripts/test/pharmcat_pgx_vcf/\'")

    # normalize variant representations and reconstruct multi-allelic variants in the input VCF
    run(args)
