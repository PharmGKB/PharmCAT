#!/usr/bin/env python3

__author__ = 'BinglanLi'

import os
import re
import sys
from pathlib import Path
from timeit import default_timer as timer

import vcf_preprocess_utilities as util


# expected tool versions
MIN_BCFTOOLS_VERSION = '1.16'
MIN_BGZIP_VERSION = '1.16'
# default filenames
PHARMCAT_POSITIONS = 'pharmcat_positions.vcf.bgz'
REFERENCE_FASTA = 'reference.fna.bgz'


def preprocess(bcftools_path, bgzip_path, pharmcat_positions_vcf, reference_genome, input_vcf, input_list,
               sample_file, output_dir, base_filename, keep_intermediate_files=False, missing_to_ref=False,
               concurrent_mode=False, max_processes=1, verbose=False):
    """ normalize and prepare the input VCF for PharmCAT """
    if not input_vcf and not input_list:
        print("Missing VCF input")
        sys.exit(1)
    start = timer()

    print("Saving output to", output_dir)
    if input_vcf:
        # compress if the input vcf is not bgzipped
        input_vcf = util.bgzip_vcf(bgzip_path, input_vcf, verbose)

    # prep pharmcat_positions helper files
    util.prep_pharmcat_positions(pharmcat_positions_vcf, bcftools_path, reference_genome, verbose)

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
    # get input basename
    if input_list:
        input_basename = os.path.basename(os.path.splitext(input_list)[0])
    else:
        input_basename = util.get_vcf_prefix(input_vcf)

    """
    normalize and prepare vcf for PharmCAT
    """
    # shrink input VCF down to PGx allele defining regions and selected samples
    # modify input VCF chromosomes naming format to <chr##>
    if input_list:
        vcf_pgx_regions = util.extract_regions_from_multiple_files(bcftools_path, bgzip_path, input_list,
                                                                   pharmcat_positions_vcf, output_dir, input_basename,
                                                                   sample_list, verbose)
    else:
        vcf_pgx_regions = util.extract_regions_from_single_file(bcftools_path, input_vcf, pharmcat_positions_vcf,
                                                                output_dir, input_basename, sample_list, verbose)
    tmp_files_to_be_removed.append(vcf_pgx_regions)

    # normalize the input VCF
    vcf_normalized = util.normalize_vcf(bcftools_path, vcf_pgx_regions, reference_genome, output_dir, verbose)
    tmp_files_to_be_removed.append(vcf_normalized)

    # extract the specific PGx genetic variants in the reference PGx VCF
    # this step also generates a report of missing PGx positions in the input VCF
    vcf_normalized_pgx_only = util.filter_pgx_variants(bcftools_path, bgzip_path, vcf_normalized, reference_genome,
                                                       pharmcat_positions_vcf, missing_to_ref, output_dir,
                                                       input_basename, verbose)
    tmp_files_to_be_removed.append(vcf_normalized_pgx_only)

    # output PharmCAT-ready single-sample VCF
    # retain only the PharmCAT allele defining positions in the output VCF file
    print()
    util.output_pharmcat_ready_vcf(bcftools_path, vcf_normalized_pgx_only, output_dir, base_filename, sample_list,
                                   concurrent_mode=concurrent_mode, max_processes=max_processes)

    # remove intermediate files
    if not keep_intermediate_files:
        print()
        print("Removing intermediate files...")
        for single_path in tmp_files_to_be_removed:
            util.remove_vcf_and_index(single_path, verbose=verbose)

    end = timer()
    print()
    print("Done.")
    print("Preprocessed input VCF in %.2f seconds" % (end - start))


if __name__ == "__main__":
    import argparse

    # describe the tool
    parser = argparse.ArgumentParser(description='Prepare an input VCF for the PharmCAT')

    # list arguments
    parser.add_argument("-vcf", "--vcf", type=str, required=True,
                        help="Path to a VCF file or a file of paths to VCF files (one file per line), "
                             "sorted by chromosome position.")
    parser.add_argument("-refVcf", "--reference-pgx-vcf", type=str,
                        default=os.path.join(os.getcwd(), PHARMCAT_POSITIONS),
                        help="(Optional) A sorted VCF of PharmCAT PGx variants, gzipped with preprocessor scripts. "
                             "Default = \'" + PHARMCAT_POSITIONS + "\' in the current working directory.")
    parser.add_argument("-refFna", "--reference-genome",
                        help="(Optional) the Human Reference Genome GRCh38/hg38 in the fasta format.")
    parser.add_argument("-S", "--sample-file",
                        help="(Optional) a file of samples to be prepared for the PharmCAT, one sample at a line.")
    parser.add_argument("-bcftools", "--path-to-bcftools",
                        help="(Optional) an alternative path to the executable bcftools.")
    parser.add_argument("-bgzip", "--path-to-bgzip",
                        help="(Optional) an alternative path to the executable bgzip.")
    parser.add_argument("-o", "--output-dir", type=str,
                        help="(Optional) directory for outputs, by default, directory of the first input VCF.")
    parser.add_argument("-bf", "--base-filename", type=str,
                        help="(Optional) output prefix (without file extensions), "
                             "by default the same base name as the input.")
    parser.add_argument("-k", "--keep-intermediate-files", action='store_true',
                        help="(Optional) keep intermediate files, false by default.")
    parser.add_argument("-0", "--missing-to-ref", action='store_true',
                        help="(Optional) assume genotypes at missing PGx sites are 0/0.  DANGEROUS!.")

    parser.add_argument("-c", "--concurrent-mode", action="store_true",
                        help="(Optional) use multiple processes - maximum number of processes spawned will default to "
                             "to one less than the number of cpu cores.")
    parser.add_argument("-cp", "--max-processes", type=int, default=None,
                        help="(Optional) maximum number of processes to use when multiprocessing.")

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="(Optional) print more verbose messages")

    # parse arguments
    args = parser.parse_args()

    m_bcftools_path = util.validate_bcftools(args.path_to_bcftools, MIN_BCFTOOLS_VERSION)
    m_bgzip_path = util.validate_bgzip(args.path_to_bgzip, MIN_BGZIP_VERSION)

    # print warnings here
    # alternatively, could use the "warnings" module
    if args.missing_to_ref:
        print("""
        =============================================================
        Warning: Argument "-0"/"--missing-to-ref" supplied
              
        THIS SHOULD ONLY BE USED IF: you sure your data is reference
        at the missing positions instead of unreadable/uncallable at
        those positions.
        
        Running PharmCAT with positions as missing vs reference can
        lead to different results.
        =============================================================

        """)

    # validate input vcf or file list
    m_input_vcf = None
    m_input_list = None
    # check whether input is a vcf or a list
    if re.search('[.]vcf([.]b?gz)?$', os.path.basename(args.vcf)):
        m_input_vcf = args.vcf
        if not os.path.exists(m_input_vcf):
            print("Cannot find VCF", m_input_vcf)
            sys.exit(1)
    else:
        m_input_list = args.vcf
        if not os.path.exists(m_input_list):
            print("Cannot find", m_input_list)
            sys.exit(1)

    # validate sample file
    m_sample_file = None
    if args.sample_file:
        m_sample_file = args.sample_file
        if not os.path.isfile(m_sample_file):
            print("Cannot find", m_sample_file)
            sys.exit(1)

    # define output base name, default to empty string
    m_base_filename = args.base_filename if args.base_filename else ''
    # define working directory, defaulting to the directory of the input VCF
    if args.output_dir:
        m_output_dir = args.output_dir
        if os.path.exists(m_output_dir):
            if not os.path.isdir(m_output_dir):
                print(m_output_dir, "is not a directory")
                sys.exit(1)
        else:
            # create the output folder
            Path(m_output_dir).mkdir(parents=True, exist_ok=True)
    elif m_input_vcf:
        m_output_dir = os.path.dirname(os.path.realpath(m_input_vcf))
    else:
        m_output_dir = os.path.dirname(os.path.realpath(m_input_list))

    # make sure we have pharmcat_positions.vcf.bgz
    if os.path.exists(args.reference_pgx_vcf):
        m_pharmcat_positions_vcf = args.reference_pgx_vcf
    else:
        # check script dir
        m_pharmcat_positions_vcf = os.path.join(os.path.dirname(os.path.realpath(__file__)), PHARMCAT_POSITIONS)
        if not os.path.exists(m_pharmcat_positions_vcf):
            print('Error: Cannot find', PHARMCAT_POSITIONS, 'in current working directory or script directory')
            sys.exit(1)

    # make sure we have reference FASTA
    if args.reference_genome:
        m_reference_genome = args.reference_genome
        if not os.path.isfile(m_reference_genome):
            print('Error: Cannot find reference genome at', m_reference_genome)
            sys.exit(1)
    else:
        ref_dirs = list({
            os.path.dirname(m_pharmcat_positions_vcf),
            os.path.dirname(os.path.realpath(__file__)),
            os.getcwd(),
        })
        m_reference_genome = None
        for ref_dir in ref_dirs:
            ref_file = os.path.join(str(ref_dir), REFERENCE_FASTA)
            if os.path.isfile(ref_file):
                m_reference_genome = ref_file
                if args.verbose:
                    print("Using reference FASTA at", m_reference_genome)
                break
        if m_reference_genome is None:
            m_reference_genome = util.get_default_grch38_ref_fasta_and_index(ref_dirs[0])
            print("Downloading reference FASTA to %s" % m_reference_genome)

    m_max_processes = 1
    if args.concurrent_mode:
        print("Using concurrent mode...")
        if args.max_processes is not None:
            if args.max_processes < 2:
                print("Cannot specify less than 2 processes in concurrent mode")
                sys.exit(1)
            if args.max_processes > os.cpu_count():
                print("Warning:", args.max_processes, "processes requested, but only", os.cpu_count(),
                      "CPUs available")
        else:
            m_max_processes = max(1, os.cpu_count() - 1)
    elif args.max_processes is not None:
        print("-cp/--max_processes will be ignored (not running in multiprocess mode)")

    # normalize variant representations and reconstruct multi-allelic variants in the input VCF
    preprocess(bcftools_path=m_bcftools_path,
               bgzip_path=m_bgzip_path,
               pharmcat_positions_vcf=m_pharmcat_positions_vcf,
               reference_genome=m_reference_genome,
               input_vcf=m_input_vcf,
               input_list=m_input_list,
               sample_file=m_sample_file,
               output_dir=m_output_dir,
               base_filename=m_base_filename,
               keep_intermediate_files=args.keep_intermediate_files,
               missing_to_ref=args.missing_to_ref,
               concurrent_mode=args.concurrent_mode,
               max_processes=m_max_processes,
               verbose=args.verbose,
               )
