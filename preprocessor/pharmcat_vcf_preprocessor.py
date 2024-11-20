#!/usr/bin/env python3

__author__ = 'BinglanLi'

import sys
from pathlib import Path
from timeit import default_timer as timer
from typing import List, Optional

import preprocessor


if __name__ == "__main__":
    import argparse

    # describe the tool
    parser = argparse.ArgumentParser(description='Prepares a VCF file for use by PharmCAT.')

    # list arguments
    parser.add_argument("-vcf", "--vcf", type=str, required=True, metavar='<file>',
                        help="Path to a VCF file or a file of paths to VCF files (one file per line), "
                             "sorted by chromosome position.")
    sample_group = parser.add_mutually_exclusive_group()
    sample_group.add_argument('-s', '--samples', type=str, metavar='<samples>',
                              help='A comma-separated list of sample IDs.  '
                                   'Only applicable if you have multiple samples and only want to work on specific ones.')
    sample_group.add_argument('-S', '--sample-file', type=str, metavar='<txt_file>',
                              help='A file containing a list of sample IDs, one sample at a line.  '
                                   'Only applicable if you have multiple samples and only want to work on specific ones.')

    parser.add_argument("-0", "--missing-to-ref", action='store_true',
                        help="Assume genotypes at missing PGx sites are 0/0.  DANGEROUS!")
    parser.add_argument("--absent-to-ref", action='store_true',
                        help="Assume genotypes at absent PGx sites are 0/0.  DANGEROUS!")
    parser.add_argument("--unspecified-to-ref", action='store_true',
                        help="Assume unspecified genotypes ./. as 0/0 when every sample is \'./.\'. DANGEROUS!")
    parser.add_argument("-G", "--no-gvcf-check", action="store_true",
                        help="Bypass check if VCF file is in gVCF format.")
    # output args
    output_group = parser.add_argument_group('Output arguments')
    output_group.add_argument("-o", "--output-dir", type=str, metavar='<dir>',
                              help="Directory for outputs.  Defaults to the directory of the input VCF.")
    output_group.add_argument("-bf", "--base-filename", type=str, metavar='<name>',
                              help="Prefix for output file names.  Defaults to the same base name as the input file.")
    output_group.add_argument("-ss", "--single-samples", action="store_true",
                              help="Generate 1 VCF file per sample.")
    output_group.add_argument("-k", "--keep-intermediate-files", action='store_true',
                              help="Keep intermediate files, false by default.")
    # concurrency args
    concurrency_group = parser.add_argument_group('Concurrency arguments')
    concurrency_group.add_argument("-c", "--concurrent-mode", action="store_true",
                                   help="Use multiple processes - maximum number of processes spawned will "
                                        "default to to two less than the number of cpu cores.")
    concurrency_group.add_argument("-cp", "--max-concurrent-processes", type=int, metavar='<num processes>',
                                   default=None,
                                   help='The maximum number of processes to use when concurrent mode ' +
                                        'is enabled.')
    # advanced args
    advanced_group = parser.add_argument_group('Advanced arguments')
    advanced_group.add_argument("-refVcf", "--reference-pgx-vcf", type=str, metavar='<vcf_file>',
                                help='A sorted, compressed VCF of PharmCAT PGx variants.  Defaults to "' +
                                     preprocessor.PHARMCAT_POSITIONS_FILENAME + '" in the current working directory ' +
                                     'or the directory the preprocessor is in.')
    advanced_group.add_argument("-refFna", "--reference-genome", type=str, metavar='<fna_file>',
                                help="The Human Reference Genome GRCh38/hg38 in the fasta format.")
    advanced_group.add_argument("-R", "--retain-specific-regions",
                                type=str, metavar='<bed_file>',
                                help='A sorted .bed file indicating regions to retain in VCF. '
                                     'For research use only. '
                                     'Additional variants are not used by PharmCAT and will slow PharmCAT down.')
    advanced_group.add_argument("-bcftools", "--path-to-bcftools", type=str, metavar='</path/to/bcftools>',
                                help="Path to the bcftools program.  Defaults to bcftools in PATH.")
    advanced_group.add_argument("-bgzip", "--path-to-bgzip", type=str, metavar='</path/to/bgzip>',
                                help="Path to the bgzip program.  Defaults to bgzip in PATH.")

    parser.add_argument("-v", "--verbose", action="count", default=0,
                        help="print more verbose messages")
    parser.add_argument('-V', '--version', action='version',
                        version='PharmCAT VCF Preprocessor v%s' % preprocessor.PHARMCAT_VERSION)

    # parse arguments
    args = parser.parse_args()

    # print the version number
    print('PharmCAT VCF Preprocessor version: %s' % preprocessor.PHARMCAT_VERSION)
    # print warnings here
    if args.missing_to_ref:
        print("""
        =============================================================
        Warning: Argument "-0"/"--missing-to-ref" supplied
              
        THIS SHOULD ONLY BE USED IF: you sure your data is reference
        at the missing positions instead of unreadable/uncallable at
        those positions.
        
        Running PharmCAT with positions as missing vs reference can
        lead to different results.
        
        "--missing-to-ref" is equivalent to providing both 
        "--absent-to-ref" and "--unspecified-to-ref".
        =============================================================

        """)
    if args.absent_to_ref:
        print("""
        =============================================================
        Warning: Argument "--absent-to-ref" supplied
              
        THIS SHOULD ONLY BE USED IF: you are sure your data is reference
        at the absent positions, instead of unreadable/uncallable.
        
        Running PharmCAT with positions as absent vs reference can 
        lead to different results.
        =============================================================

        """)
    if args.unspecified_to_ref:
        print("""
        =============================================================
        Warning: Argument "--unspecified-to-ref" supplied
        
        THIS SHOULD ONLY BE USED IF: you are confident in assuming
        the unspecified genotypes ./. in your data as reference.
        
        Running PharmCAT with positions as unspecified vs reference can
        lead to different results.
        =============================================================

        """)

    try:
        # make sure we have required tools
        m_bcftools_path = preprocessor.validate_bcftools(args.path_to_bcftools)
        m_bgzip_path = preprocessor.validate_bgzip(args.path_to_bgzip)

        script_dir: Path = Path(globals().get("__file__", "./_")).absolute().parent

        # make sure we have pharmcat_positions.vcf.bgz
        m_pharmcat_positions_vcf: Path
        if args.reference_pgx_vcf:
            m_pharmcat_positions_vcf = preprocessor.validate_file(args.reference_pgx_vcf)
        else:
            m_pharmcat_positions_vcf = preprocessor.find_file(preprocessor.PHARMCAT_POSITIONS_FILENAME,
                                                              [Path.cwd(), script_dir])
            if m_pharmcat_positions_vcf is None:
                print('Downloading pharmcat_positions.vcf...')
                m_pharmcat_positions_vcf = preprocessor.download_pharmcat_accessory_files(script_dir,
                                                                                          verbose=args.verbose)

        # make sure we have reference FASTA
        m_reference_genome: Path
        if args.reference_genome:
            m_reference_genome = preprocessor.validate_file(args.reference_genome)
        else:
            m_reference_genome = preprocessor.find_file(preprocessor.REFERENCE_FASTA_FILENAME,
                                                        [m_pharmcat_positions_vcf.parent, Path.cwd(), script_dir])
            if m_reference_genome is None:
                print('Downloading reference FASTA.  This may take a while...')
                m_reference_genome = preprocessor.download_reference_fasta_and_index(m_pharmcat_positions_vcf.parent,
                                                                                     verbose=args.verbose)

        # validate .bed file if retaining specific regions
        m_regions_bed: Optional[Path] = None
        if args.retain_specific_regions:
            m_regions_bed = preprocessor.validate_file(args.retain_specific_regions)

        # prep pharmcat_positions helper files
        preprocessor.prep_pharmcat_positions(m_pharmcat_positions_vcf, m_reference_genome, verbose=args.verbose)

        # validate input vcf or file list
        vcf_path: Path = Path(args.vcf)
        if vcf_path.is_dir():
            print('%s is a directory not a file' % vcf_path)
            sys.exit(1)
        elif not vcf_path.is_file():
            print("Error: no VCF input")
            sys.exit(1)

        m_vcf_files: List[Path] = []
        m_input_basename: str
        if preprocessor.is_vcf_file(vcf_path):
            m_vcf_files.append(vcf_path)
            m_input_basename = preprocessor.get_vcf_basename(vcf_path)
        else:
            if args.verbose:
                print("Looking up VCF files listed in", vcf_path)
            with open(vcf_path, 'r') as in_f:
                for line in in_f:
                    line = line.strip()
                    m_vcf_files.append(preprocessor.validate_file(line))
            m_input_basename = vcf_path.stem

        if len(m_vcf_files) == 0:
            print("Error: no VCF input")
            sys.exit(1)

        # check whether input is a gVCF, which currently does not support yet
        if args.no_gvcf_check:
            print('\nBypass the gVCF check.\n')
        else:
            for file in m_vcf_files:
                if preprocessor.is_gvcf_file(file):
                    print('%s is a gVCF file, which is not currently supported.\n'
                          'See https://github.com/PharmGKB/PharmCAT/issues/79 for details.\n'
                          'If this is not a gVCF file, use -G to bypass.' % str(file))
                    sys.exit(1)

        # check whether a sample file is provided, otherwise proceed with all samples in the input VCF
        m_samples: List[str] = []
        if args.sample_file:
            # validate sample file
            m_samples = preprocessor.read_sample_file(preprocessor.validate_file(args.sample_file),
                                                      verbose=args.verbose)
        elif args.samples:
            m_samples = preprocessor.parse_samples(args.samples)
        if len(m_samples) == 0:
            m_samples = preprocessor.read_vcf_samples(m_vcf_files[0], verbose=args.verbose)

        # make sure xx-to-ref commands are provided correctly
        m_absent_to_ref: bool
        m_unspecified_to_ref: bool
        if args.missing_to_ref:
            m_absent_to_ref = True
            m_unspecified_to_ref = True
        else:
            m_absent_to_ref = args.absent_to_ref
            m_unspecified_to_ref = args.unspecified_to_ref

        # define output base name, default to empty string
        m_output_basename: str = args.base_filename if args.base_filename else ''
        # define working directory, defaulting to the directory of the input VCF
        m_output_dir: Path
        if args.output_dir:
            m_output_dir = preprocessor.validate_dir(args.output_dir, create_if_not_exist=True)
        else:
            m_output_dir = vcf_path.parent

        m_max_processes: int = 1
        if args.concurrent_mode:
            print('Concurrent mode enabled...')
            m_max_processes = preprocessor.check_max_processes(args.max_concurrent_processes)
        elif args.max_concurrent_processes is not None:
            print("-cp/--max_processes will be ignored (not running in multiprocess mode)")

        # get pgx_regions we care about
        custom_regions: bool = bool(m_regions_bed)
        pgx_regions = preprocessor.get_pgx_regions(m_regions_bed or m_pharmcat_positions_vcf)

        start = timer()

        # normalize variant representations and reconstruct multi-allelic variants in the input VCF
        if args.verbose:
            print("Using reference FASTA at", m_reference_genome)
        print("Saving output to", m_output_dir.absolute())
        print()
        results = preprocessor.preprocess(pharmcat_positions_vcf=m_pharmcat_positions_vcf,
                                          reference_genome=m_reference_genome,
                                          regions_to_retain=pgx_regions,
                                          custom_regions=custom_regions,
                                          vcf_files=m_vcf_files,
                                          samples=m_samples,
                                          input_basename=m_input_basename,
                                          output_dir=m_output_dir,
                                          output_basename=m_output_basename,
                                          split_samples=args.single_samples,
                                          keep_intermediate_files=args.keep_intermediate_files,
                                          absent_to_ref=m_absent_to_ref,
                                          unspecified_to_ref=m_unspecified_to_ref,
                                          concurrent_mode=args.concurrent_mode,
                                          max_processes=m_max_processes,
                                          verbose=args.verbose,
                                          )
        end = timer()
        if len(results) == 1:
            print()
            print('Generated PharmCAT-ready VCF:', str(results[0].absolute()))
        elif args.verbose:
            print()
            print('Generated PharmCAT-ready VCFs:')
            for rez in results:
                print('* %s' % str(rez))
        print()
        print("Done.")
        print("Preprocessed input VCF in %.2f seconds" % (end - start))

    except preprocessor.ReportableException as e:
        print(e)
        sys.exit(1)
