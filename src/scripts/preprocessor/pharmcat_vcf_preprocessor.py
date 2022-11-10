#!/usr/bin/env python3

__author__ = 'BinglanLi'

import os
import sys
from pathlib import Path
from typing import Optional, List

import preprocessor


if __name__ == "__main__":
    import argparse

    # describe the tool
    parser = argparse.ArgumentParser(description='Prepares a VCF file for use by PharmCAT.')

    # list arguments
    parser.add_argument("-vcf", "--vcf", type=str, required=True, metavar='<vcf file>',
                        help="path to a VCF file or a file of paths to VCF files (one file per line), "
                             "sorted by chromosome position.")
    parser.add_argument("-refVcf", "--reference-pgx-vcf", type=str, metavar='<vcf_file>',
                        help='(Optional) a sorted, compressed VCF of PharmCAT PGx variants.  Defaults to "' +
                             preprocessor.PHARMCAT_POSITIONS_FILENAME + '" in the current working directory or the ' +
                             'directory the preprocessor is in.')
    parser.add_argument("-refFna", "--reference-genome", type=str, metavar='<fna_file>',
                        help="(Optional) the Human Reference Genome GRCh38/hg38 in the fasta format.")
    parser.add_argument("-S", "--sample-file", type=str, metavar='<txt_file>',
                        help="(Optional) a file of samples to be prepared for the PharmCAT, one sample at a line.")
    parser.add_argument("-bcftools", "--path-to-bcftools", type=str, metavar='</path/to/bcftools>',
                        help="(Optional) a path to the bcftools program.  Defaults to bcftools in PATH.")
    parser.add_argument("-bgzip", "--path-to-bgzip", type=str, metavar='</path/to/bgzip>',
                        help="(Optional) a path to the bgzip program.  Defaults to bgzip in PATH.")
    parser.add_argument("-o", "--output-dir", type=str, metavar='<dir>',
                        help="(Optional) directory for outputs.  Defaults to the directory of the input VCF.")
    parser.add_argument("-bf", "--base-filename", type=str, metavar='<name>',
                        help="(Optional) output prefix (without file extensions), "
                             "by default the same base name as the input.")
    parser.add_argument("-k", "--keep-intermediate-files", action='store_true',
                        help="(Optional) keep intermediate files, false by default.")
    parser.add_argument("-0", "--missing-to-ref", action='store_true',
                        help="(Optional) assume genotypes at missing PGx sites are 0/0.  DANGEROUS!.")

    parser.add_argument("-c", "--concurrent-mode", action="store_true",
                        help="(Optional) use multiple processes - maximum number of processes spawned will default to "
                             "to one less than the number of cpu cores.")
    parser.add_argument("-cp", "--max-concurrent-processes", type=int, metavar='<num processes>', default=None,
                        help="(Optional) the maximum number of processes to use when concurrent mode is enabled.")

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="(Optional) print more verbose messages")
    parser.add_argument('-V', '--version', action='version',
                        version='PharmCAT VCF Preprocessor v%s' % preprocessor.PHARMCAT_VERSION)

    # parse arguments
    args = parser.parse_args()

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
                                                              list({Path.cwd(), script_dir}))
            if m_pharmcat_positions_vcf is None:
                print('Error: Cannot find %s in current working directory or script directory' %
                      preprocessor.PHARMCAT_POSITIONS_FILENAME)
                sys.exit(1)
        # make sure we have reference FASTA
        m_reference_genome: Path
        if args.reference_genome:
            m_reference_genome = preprocessor.validate_file(args.reference_genome)
        else:
            m_reference_genome = preprocessor.find_file(preprocessor.REFERENCE_FASTA_FILENAME,
                                                        list({m_pharmcat_positions_vcf.parent, Path.cwd(), script_dir}))
            if m_reference_genome is None:
                m_reference_genome = preprocessor.download_reference_fasta_and_index(m_pharmcat_positions_vcf.parent,
                                                                                     verbose=args.verbose)
                print('Downloading reference FASTA.  This may take a while...')

        # validate input vcf or file list
        m_vcf_files: List[Path] = []
        vcf_path: Path = Path(args.vcf)
        m_input_basename: str
        if vcf_path.is_dir():
            m_vcf_files = preprocessor.find_vcf_files(vcf_path, args.verbose)
            m_input_basename = vcf_path.name
        elif vcf_path.is_file():
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
        else:
            print("Error:", vcf_path, "does not exist")
            sys.exit(1)
        if len(m_vcf_files) == 0:
            print("Error: no VCF input")
            sys.exit(1)

        m_samples: List[str] = []
        if args.sample_file:
            # validate sample file
            sample_file: Optional[Path] = preprocessor.validate_file(args.sample_file)
            print("Reading samples from", sample_file, "...")
            with open(sample_file, 'r') as file:
                for line in file:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        m_samples.append(line)
            if len(m_samples) == 0:
                print("Warning: No samples found. Will use all samples listed in VCF.")
            else:
                preprocessor.validate_samples(m_samples)
        if len(m_samples) == 0:
            m_samples = preprocessor.read_vcf_samples(m_vcf_files[0])

        # define output base name, default to empty string
        m_output_basename: str = args.base_filename if args.base_filename else ''
        # define working directory, defaulting to the directory of the input VCF
        m_output_dir: Path
        if args.output_dir:
            m_output_dir = Path(args.output_dir)
            if m_output_dir.exists():
                if not m_output_dir.is_dir():
                    print(m_output_dir, "is not a directory")
                    sys.exit(1)
            else:
                # create the output folder
                m_output_dir.mkdir(parents=True, exist_ok=True)
        else:
            m_output_dir = vcf_path.parent

        m_max_processes: int = 1
        if args.concurrent_mode:
            if args.max_concurrent_processes is not None:
                if args.max_concurrent_processes < 2:
                    print("Error: Cannot specify less than 2 processes in concurrent mode")
                    sys.exit(1)
                if args.max_concurrent_processes > os.cpu_count():
                    print("Warning:", args.max_concurrent_processes, "processes requested, but only", os.cpu_count(),
                          "CPUs available")
                m_max_processes = args.max_concurrent_processes
            else:
                m_max_processes = max(1, os.cpu_count() - 1)
        elif args.max_concurrent_processes is not None:
            print("-cp/--max_processes will be ignored (not running in multiprocess mode)")

        # normalize variant representations and reconstruct multi-allelic variants in the input VCF
        preprocessor.preprocess(pharmcat_positions_vcf=m_pharmcat_positions_vcf,
                                reference_genome=m_reference_genome,
                                vcf_files=m_vcf_files,
                                samples=m_samples,
                                input_basename=m_input_basename,
                                output_dir=m_output_dir,
                                output_basename=m_output_basename,

                                keep_intermediate_files=args.keep_intermediate_files,
                                missing_to_ref=args.missing_to_ref,
                                concurrent_mode=args.concurrent_mode,
                                max_processes=m_max_processes,
                                verbose=args.verbose,
                                )
    except preprocessor.ReportableException as e:
        print(e)
        sys.exit(1)
