#! /usr/bin/env python

__author__ = 'BinglanLi'

import os
import subprocess

def run(args):
    # assign the current working directory as the output path if not specified
    output_folder = os.getcwd() if not args.output_folder else args.output_folder
    # organize args
    output_full_name = os.path.join(output_folder, args.output_prefix + ".vcf.gz")
    tabix_executable_path = args.path_to_tabix if args.path_to_tabix else "tabix"
    bcftools_executable_path = args.path_to_bcftools if args.path_to_bcftools else "bcftools"

    # validate the input arguments, write another function
    #validate(args)

    # if the input VCF file is not index (.tbi doesn't exist), index the file using tabix
    if not os.path.exists(args.input_vcf + ".tbi"):
        saved_wd = os.getcwd()
        input_folder = os.path.split(args.input_vcf)[0]
        subprocess.run([tabix_executable_path, "-p", "vcf", args.input_vcf], cwd = input_folder)
        os.chdir(saved_wd)

    # rename chromosomes
    output_temp = os.path.join(output_folder, args.output_prefix + ".temp_chr_renamed.vcf.gz")
    print(" ".join([bcftools_executable_path, "annotate", "--rename-chrs", args.rename_chrs, "-Oz", "-o", output_temp, args.input_vcf]))
    subprocess.run([bcftools_executable_path, "annotate", "--rename-chrs", args.rename_chrs, "-Oz", "-o", output_temp, args.input_vcf])
    #saved_wd = os.getcwd()
    subprocess.run([tabix_executable_path, "-p", "vcf", output_temp], cwd = output_folder)
    #os.chdir(saved_wd)

    # normalize the input VCF
    if args.extract_regions:
        subprocess.run([bcftools_executable_path, "norm", "-m+", "-Oz", "-o", output_full_name, "-R", args.extract_regions, "-f", args.ref_seq, output_temp])
    else:
        subprocess.run([bcftools_executable_path, "norm", "-m+", "-Oz", "-o", output_full_name, "-f", args.ref_seq, output_temp])

    # index the output VCF
    subprocess.run([tabix_executable_path, "-p", "vcf", output_full_name], cwd = args.output_folder)

    # remove temporary files
    try:
        os.remove(output_temp)
        os.remove(output_temp + ".tbi")
    except OSError as error_remove_tmp:
        print("Error: %s : %s" % (output_temp, error_remove_tmp.strerror))
    
    

if __name__ == "__main__":
    import argparse

    # describe the 
    parser = argparse.ArgumentParser(description='convert_to_normalized_vcf.py :  Will normalize the input VCF file, as well as reconstruct multi-allelic variants')

    # list arguments
    parser.add_argument("--input_vcf", help="Load a compressed VCF file.")
    parser.add_argument("--extract_regions", help="Extract the regions from the VCF. The region file has to be tab-delimited with CHROM, POS, and optionally END-POS.")
    parser.add_argument("--ref_seq", help="Load the Human Reference Genome GRCh38/hg38 in the fasta format.")
    parser.add_argument("--rename_chrs", help="Load a chromosome rename map file")
    parser.add_argument("--path_to_bcftools", help="Load an alternative path to the executable bcftools.")
    parser.add_argument("--path_to_tabix", help="Load an alternative path to the executable tabix.")
    parser.add_argument("--output_folder", help="Directory of the output VCF, by default, current working directory.")
    parser.add_argument("--output_prefix", help="Name of the output VCF")

    # parse arguments
    args = parser.parse_args()

    # normalize variant representations and reconstruct multi-allelic variants in the input VCF
    run(args)