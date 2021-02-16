#! /usr/bin/env python

__author__ = 'BinglanLi'

import copy
import os
import vcf_preprocess_utilities as Utilities

import convert_to_normalized_vcf
import convert_to_pharmcat_vcf

def run(args):
    ## this section first normalizes and prepares the input VCF file for the PharmCAT

    # normalize the input VCF file
    convert_to_normalized_vcf_args = copy.copy(args)
    convert_to_normalized_vcf_args.output_folder = os.getcwd() if not args.output_folder else args.output_folder
    convert_to_normalized_vcf_args.output_prefix = Utilities.obtain_vcf_file_prefix(args.input_vcf) + ".temp_normalized"
    convert_to_normalized_vcf.run(convert_to_normalized_vcf_args)

    # convert the normalized VCF to PharmCAT-ready format
    convert_to_pharmcat_vcf_args = copy.copy(args)
    convert_to_pharmcat_vcf_args.input_vcf = os.path.join(convert_to_normalized_vcf_args.output_folder, convert_to_normalized_vcf_args.output_prefix + ".vcf.gz")
    convert_to_pharmcat_vcf.run(convert_to_pharmcat_vcf_args)

    # remove temporary intermediate files
    try:
        os.remove(convert_to_pharmcat_vcf_args.input_vcf)
        os.remove(convert_to_pharmcat_vcf_args.input_vcf + ".tbi")
    except OSError as error_remove_tmp:
        print("Error: %s : %s" % (convert_to_pharmcat_vcf_args.input_vcf, error_remove_tmp.strerror))


if __name__ == "__main__":
    import argparse

    # describe the 
    parser = argparse.ArgumentParser(description='Prepare an input VCF for the PharmCAT')

    # list arguments
    parser.add_argument("--input_vcf", required=True, type = str, help="Load a compressed VCF file.")
    parser.add_argument("--ref_seq", help="Load the Human Reference Genome GRCh38/hg38 in the fasta format.")
    parser.add_argument("--rename_chrs", help="Load a chromosome rename map file. This is a must for VCF normalization if the chromosomes are not named 'chr##'.")
    parser.add_argument("--ref_pgx_vcf", required=True, type = str, help="Load a VCF file of PGx variants. This file is available from the PharmCAT GitHub release.")
    parser.add_argument("--path_to_bcftools", help="Load an alternative path to the executable bcftools.")
    parser.add_argument("--path_to_tabix", help="Load an alternative path to the executable tabix.")
    parser.add_argument("--output_folder", help="Directory of the output VCF, by default, current working directory.")
    parser.add_argument("--output_prefix", required=True, type = str, help="Prefix of the output VCF")

    # parse arguments
    args = parser.parse_args()

    # normalize variant representations and reconstruct multi-allelic variants in the input VCF
    run(args)
