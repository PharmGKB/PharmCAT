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
    ## "bcftools norm <options> <input_vcf>" normalizes the input VCF file. "-m+" join biallelic sites into multiallelic records. "-f" reference genome sequence in fasta format.

    # set up variables
    url_grch38_fasta = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz'
    url_grch38_fasta_index = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai'
    ## organize args
    current_working_dir = os.getcwd()
    output_folder = os.getcwd() if not args.output_folder else args.output_folder
    output_full_name = os.path.join(output_folder, args.output_prefix + ".vcf.gz")
    path_to_ref_seq = args.ref_seq
    tabix_executable_path = args.path_to_tabix if args.path_to_tabix else "tabix"
    bcftools_executable_path = args.path_to_bcftools if args.path_to_bcftools else "bcftools"

    # download the human reference sequence if not provided
    if not path_to_ref_seq:
        path_to_ref_seq = Utilities.download_from_url(url_grch38_fasta, current_working_dir)
        subprocess.run(['gunzip', path_to_ref_seq], cwd = current_working_dir)
        path_to_ref_seq = '.'.join(path_to_ref_seq.split('.')[:-1])
        Utilities.download_from_url(url_grch38_fasta_index, current_working_dir)

    # if the input VCF file is not indexed (.tbi doesn't exist), index the file using tabix
    if not os.path.exists(args.input_vcf + ".tbi"):
        input_folder = os.path.split(args.input_vcf)[0]
        subprocess.run([tabix_executable_path, "-p", "vcf", args.input_vcf], cwd = input_folder)

    # rename chromosomes
    input_renamed_chr = os.path.join(output_folder, args.output_prefix + ".temp_chr_renamed.vcf.gz")
    subprocess.run([bcftools_executable_path, "annotate", '--no-version', "--rename-chrs", args.rename_chrs, "-Oz", "-o", input_renamed_chr, args.input_vcf])
    subprocess.run([tabix_executable_path, "-p", "vcf", input_renamed_chr], cwd = output_folder)

    # normalize the input VCF
    # modify this part to comply to the PharmCAT VCF requirements and PharmCAT only
    subprocess.run([bcftools_executable_path, "norm", "-m+", "-Oz", "-o", output_full_name, "-f", path_to_ref_seq, input_renamed_chr])

    # index the output VCF
    subprocess.run([tabix_executable_path, "-p", "vcf", output_full_name], cwd = output_folder)

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
    parser.add_argument("--input_vcf", required=True, type = str, help="Load a compressed VCF file.")
    parser.add_argument("--ref_seq", help="Load the Human Reference Genome GRCh38/hg38 in the fasta format.")
    parser.add_argument("--rename_chrs", help="Load a chromosome rename map file")
    parser.add_argument("--path_to_bcftools", help="Load an alternative path to the executable bcftools.")
    parser.add_argument("--path_to_tabix", help="Load an alternative path to the executable tabix.")
    parser.add_argument("--output_folder", help="Directory of the output VCF, by default, current working directory.")
    parser.add_argument("--output_prefix", required=True, type = str, help="Prefix of the output VCF")

    # parse arguments
    args = parser.parse_args()

    # normalize variant representations and reconstruct multi-allelic variants in the input VCF
    run(args)
