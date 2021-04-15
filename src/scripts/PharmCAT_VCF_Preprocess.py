#! /usr/bin/env python

__author__ = 'BinglanLi'

import os
from pathlib import Path
from timeit import default_timer as timer

import vcf_preprocess_utilities as util

def run(args):
    ## normalize and prepare the input VCF for PharmCAT
    start = timer()

    # organize args
    current_working_dir = os.getcwd()
    tabix_executable_path = args.path_to_tabix if args.path_to_tabix else 'tabix'
    bcftools_executable_path = args.path_to_bcftools if args.path_to_bcftools else 'bcftools'
    # read the sample list
    if args.sample_file:
        sample_list = []
        with open(args.sample_file, 'r') as file:
            for line in file:
                line = line.strip()
                sample_list.append(line)
        file.close()

    # create the output folder if not existing
    Path(args.output_folder).mkdir(parents=True, exist_ok=True)

    tmp_files_to_be_removed = []

    # download the human reference sequence if not provided
    path_to_ref_seq = util.download_grch38_ref_fasta_and_index(current_working_dir) if not args.ref_seq else args.ref_seq

    # if the input VCF file is not indexed (.tbi doesn't exist), create an index file in the input folder using tabix
    if not os.path.exists(args.input_vcf + '.tbi'):
        util.tabix_index_vcf(tabix_executable_path, args.input_vcf)

    # index ref_pgx_vcf if not indexed
    if not os.path.exists(args.ref_pgx_vcf + '.tbi'):
        util.tabix_index_vcf(tabix_executable_path, args.ref_pgx_vcf)

    # shrink input VCF down to PGx allele defining positions to speed up
    # modify input VCF chromosomes naming format to <chr##>
    intermediate_vcf_pgx_regions = util.extract_pharmcat_pgx_regions(tabix_executable_path, args.input_vcf, args.output_folder, args.ref_pgx_vcf, sample_list)
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
    util.output_pharmcat_ready_vcf(intermediate_vcf_normalized, args.output_folder, args.output_prefix)

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
    parser.add_argument("--input_vcf", required=True, type = str, help="Load a compressed VCF file.")
    parser.add_argument("--ref_seq", help="Load the Human Reference Genome GRCh38/hg38 in the fasta format.")
    parser.add_argument("--ref_pgx_vcf", required=True, type = str, help="Load a VCF file of PGx variants. This file is available from the PharmCAT GitHub release.")
    parser.add_argument("--sample_file", help="A file of samples to be prepared for the PharmCAT.")
    parser.add_argument("--path_to_bcftools", help="Load an alternative path to the executable bcftools.")
    parser.add_argument("--path_to_tabix", help="Load an alternative path to the executable tabix.")
    parser.add_argument("--output_folder", default = os.getcwd(), type = str, help="Directory of the output VCF, by default, current working directory.")
    parser.add_argument("--output_prefix", default = 'pharmcat_ready_vcf', type = str, help="Prefix of the output VCF")
    parser.add_argument("--keep_intermediate_files", action='store_true', help="Keep intermediate files, false by default.")

    # parse arguments
    args = parser.parse_args()

    # normalize variant representations and reconstruct multi-allelic variants in the input VCF
    run(args)
