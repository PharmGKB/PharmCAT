#! /usr/bin/env python

__author__ = 'BinglanLi'

import os
import subprocess
import vcf_preprocess_utilities as Utilities

def run(args):
    ## this section invokes external tabix and bcftools to normalize an input VCF file
    ## tabix commands are exclusively "tabix -p vcf <input_vcf>", which generates an index file (.tbi) for an input file (<input_file>) whose file type is specified by "-p vcf". The .tabi is, by default, output to the current working directory. 
    ## "bcftools <common_options> <input_vcf>". "-O" (capitalized letter O) specifies the output type as compressed VCF (z). "-o" writes to a file rather than to default standard output.
    ## "bcftools annotate <options> <input_vcf>" annotates and converts the VCF column info. "--rename-chrs <chr_name_mapping_file>" convert chromosome names from the old ones to the new ones.
    ## "bcftools norm <options> <input_vcf>" normalizes the input VCF file. "-m+" join biallelic sites into multiallelic records. "-f" reference genome sequence in fasta format. "--check-ref" when incorrect or missing REF allele is encountered: exit (e), warn (w), exclude (x), or set/fix (s) bad sites.

    # assign the current working directory as the output path if not specified
    output_folder = os.getcwd() if not args.output_folder else args.output_folder
    # organize args
    tabix_executable_path = args.path_to_tabix if args.path_to_tabix else "tabix"
    bcftools_executable_path = args.path_to_bcftools if args.path_to_bcftools else "bcftools"
    input_folder = os.path.split(args.input_vcf)[0]

    # if the input VCF file is not index (.tbi doesn't exist), index the file using tabix
    if not os.path.exists(args.input_vcf + ".tbi"):
        subprocess.run([tabix_executable_path, "-p", "vcf", args.input_vcf], cwd = input_folder)

    # obtain list of samples
    input_vcf_sample_list = subprocess.check_output([bcftools_executable_path, "query", "-l", args.input_vcf], universal_newlines=True, cwd = input_folder).split('\n')
    input_vcf_sample_list.pop()

    # generate PharmCAT-ready single-sample VCF file(s)
    for input_vcf_single_sample in input_vcf_sample_list:
        output_full_name = os.path.join(output_folder, args.output_prefix + "." + input_vcf_single_sample + ".vcf.gz")
        bcftools_command = [bcftools_executable_path, "view", "-Oz", "-o", output_full_name]
        bcftools_command.extend(["-s", input_vcf_single_sample, args.input_vcf]) if args.ref_pgx_vcf else bcftools_command.append(args.input_vcf) 

        subprocess.run(bcftools_command, cwd = args.output_folder)

    # output a VCF file of missing PGx positions
    if args.ref_pgx_vcf:
        output_missing_pos_full_name = os.path.join(output_folder, args.output_prefix + ".missing_pgx_var.vcf.gz")
        # explain bcftools flags
        bcftools_command = [bcftools_executable_path, "isec", "-c", "indels", "-w1", "-Oz", "-o", output_missing_pos_full_name, "-C", args.ref_pgx_vcf, args.input_vcf]
        subprocess.run(bcftools_command, cwd = args.output_folder)


if __name__ == "__main__":
    import argparse

    # describe the 
    parser = argparse.ArgumentParser(description='convert_to_normalized_vcf.py :  Will normalize the input VCF file, as well as reconstruct multi-allelic variants')

    # list arguments
    parser.add_argument("--input_vcf", required=True, type = str, help="Load a compressed, normalized VCF file.")
    parser.add_argument("--ref_pgx_vcf", required=True, type = str, help="Load a VCF file of PGx variants. This file is available from the PharmCAT GitHub release.")
    parser.add_argument("--path_to_bcftools", help="Load an alternative path to the executable bcftools.")
    parser.add_argument("--path_to_tabix", help="Load an alternative path to the executable tabix.")
    parser.add_argument("--output_folder", help="Directory of the output VCF, by default, current working directory.")
    parser.add_argument("--output_prefix", required=True, type = str, help="Prefix of the output VCF")

    # parse arguments
    args = parser.parse_args()

    # normalize variant representations and reconstruct multi-allelic variants in the input VCF
    run(args)
