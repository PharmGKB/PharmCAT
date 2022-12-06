import shutil
from pathlib import Path
from timeit import default_timer as timer
from typing import Optional, List

from . import utilities as util
from .exceptions import ReportableException


def preprocess(pharmcat_positions_vcf: Path, reference_genome: Path,
               vcf_files: List[Path], samples: Optional[List[str]], input_basename: str,
               output_dir: Path, output_basename: Optional[str] = '', split_samples=False,
               keep_intermediate_files=False, missing_to_ref=False,
               concurrent_mode=False, max_processes=1, verbose=False) -> None:
    """
    Normalize and prepare the input VCF for PharmCAT.
    """
    if len(vcf_files) == 0:
        raise ReportableException('Missing VCF input')

    start = timer()
    if verbose:
        print("Using reference FASTA at", reference_genome)
    print("Saving output to", output_dir.absolute())
    if concurrent_mode:
        print("Using concurrent mode with max of %d processes..." % max_processes)
    print()

    # prep pharmcat_positions helper files
    util.prep_pharmcat_positions(pharmcat_positions_vcf, reference_genome, verbose=verbose)

    # make sure we have samples
    if samples is None or len(samples) == 0:
        samples = util.read_vcf_samples(vcf_files[0])

    # list of files to be deleted
    tmp_files_to_be_removed: List[Path] = []

    # shrink input VCF down to PGx allele defining regions and selected samples
    # modify input VCF chromosomes naming format to <chr##>
    pgx_region_vcf: Path = util.extract_pgx_regions(pharmcat_positions_vcf, vcf_files, samples, output_dir,
                                                    input_basename,
                                                    concurrent_mode=concurrent_mode, max_processes=max_processes,
                                                    verbose=verbose)
    tmp_files_to_be_removed.append(pgx_region_vcf)

    # normalize the input VCF
    normalized_vcf = util.normalize_vcf(reference_genome, pgx_region_vcf, output_dir, input_basename, verbose=verbose)
    tmp_files_to_be_removed.append(normalized_vcf)

    # extract the specific PGx genetic variants in the reference PGx VCF
    # this step also generates a report of missing PGx positions in the input VCF
    pgx_variants_vcf: Path = util.extract_pgx_variants(pharmcat_positions_vcf, reference_genome, normalized_vcf,
                                                       output_dir, input_basename, missing_to_ref=missing_to_ref,
                                                       verbose=verbose)
    print()
    if split_samples:
        util.index_vcf(pgx_variants_vcf, verbose)
        tmp_files_to_be_removed.append(pgx_variants_vcf)

        # output PharmCAT-ready single-sample VCF
        # retain only the PharmCAT allele defining positions in the output VCF file
        util.output_pharmcat_ready_vcf(pgx_variants_vcf, samples, output_dir, output_basename or input_basename,
                                       concurrent_mode=concurrent_mode, max_processes=max_processes)
    else:
        final_vcf: Path = output_dir / ((output_basename or input_basename) + '.preprocessed.vcf.bgz')
        shutil.move(pgx_variants_vcf, final_vcf)
        print('Generated PharmCAT-ready VCF:', final_vcf)

    # remove intermediate files
    if not keep_intermediate_files:
        if verbose:
            print()
            print("Removing intermediate files...")
        for single_path in tmp_files_to_be_removed:
            util.delete_vcf_and_index(single_path, verbose=verbose)

    end = timer()
    print()
    print("Done.")
    print("Preprocessed input VCF in %.2f seconds" % (end - start))

