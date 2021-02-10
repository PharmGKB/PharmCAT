#! /usr/bin/env python

__author__ = 'BinglanLi'

import os
import subprocess

def run(args):
    ## this section invokes external tabix and bcftools to normalize an input VCF file
    ## tabix commands are exclusively "tabix -p vcf <input_vcf>", which generates an index file (.tbi) for an input file (<input_file>) whose file type is specified by "-p vcf". The .tabi is, by default, output to the current working directory. 
    ## "bcftools <common_options> <input_vcf>". "-O" (capitalized letter O) specifies the output type as compressed VCF (z). "-o" writes to a file rather than to default standard output.
    ## "bcftools annotate <options> <input_vcf>" annotates and converts the VCF column info. "--rename-chrs <chr_name_mapping_file>" convert chromosome names from the old ones to the new ones.
    ## "bcftools norm <options> <input_vcf>" normalizes the input VCF file. "-m+" join biallelic sites into multiallelic records. "-f" reference genome sequence in fasta format.

    # (function pending) validate the input arguments
    #validate(args)

    # assign the current working directory as the output path if not specified
    output_folder = os.getcwd() if not args.output_folder else args.output_folder
    # organize args
    output_full_name = os.path.join(output_folder, args.output_prefix + ".vcf.gz")
    tabix_executable_path = args.path_to_tabix if args.path_to_tabix else "tabix"
    bcftools_executable_path = args.path_to_bcftools if args.path_to_bcftools else "bcftools"

    # if the input VCF file is not indexed (.tbi doesn't exist), index the file using tabix
    if not os.path.exists(args.input_vcf + ".tbi"):
        input_folder = os.path.split(args.input_vcf)[0]
        subprocess.run([tabix_executable_path, "-p", "vcf", args.input_vcf], cwd = input_folder)

    # rename chromosomes
    input_renamed_chr = os.path.join(output_folder, args.output_prefix + ".temp_chr_renamed.vcf.gz")
    # explain bcftools flags (find a place in the script to explain it)
    subprocess.run([bcftools_executable_path, "annotate", "--rename-chrs", args.rename_chrs, "-Oz", "-o", input_renamed_chr, args.input_vcf])
    subprocess.run([tabix_executable_path, "-p", "vcf", input_renamed_chr], cwd = output_folder)

    # normalize the input VCF
    # modify this part to comply to the PharmCAT VCF requirements and PharmCAT only
    # explain bcftools flags
    subprocess.run([bcftools_executable_path, "norm", "-m+", "-Oz", "-o", output_full_name, "-f", args.ref_seq, input_renamed_chr])

    # index the output VCF
    subprocess.run([tabix_executable_path, "-p", "vcf", output_full_name], cwd = args.output_folder)

    # remove temporary files
    try:
        os.remove(input_renamed_chr)
        os.remove(input_renamed_chr + ".tbi")
    except OSError as error_remove_tmp:
        print("Error: %s : %s" % (input_renamed_chr, error_remove_tmp.strerror))
    
    

if __name__ == "__main__":
    import argparse

    # describe the 
    parser = argparse.ArgumentParser(description='convert_to_normalized_vcf.py :  Will normalize the input VCF file, as well as reconstruct multi-allelic variants')

    # list arguments
    parser.add_argument("--input_vcf", help="Load a compressed VCF file.")
    parser.add_argument("--ref_seq", help="Load the Human Reference Genome GRCh38/hg38 in the fasta format.")
    parser.add_argument("--rename_chrs", help="Load a chromosome rename map file")
    parser.add_argument("--path_to_bcftools", help="Load an alternative path to the executable bcftools.")
    parser.add_argument("--path_to_tabix", help="Load an alternative path to the executable tabix.")
    parser.add_argument("--output_folder", help="Directory of the output VCF, by default, current working directory.")
    parser.add_argument("--output_prefix", help="Prefix of the output VCF")

    # parse arguments
    args = parser.parse_args()

    # normalize variant representations and reconstruct multi-allelic variants in the input VCF
    run(args)
