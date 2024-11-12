import shutil
from pathlib import Path
from typing import Optional, List

from . import utilities as util
from .exceptions import ReportableException


def preprocess(pharmcat_positions_vcf: Path, reference_genome: Path,
               vcf_files: List[Path], samples: Optional[List[str]], input_basename: str,
               output_dir: Path, output_basename: Optional[str] = '', split_samples: bool = False,
               keep_intermediate_files: bool = False,
               absent_to_ref: bool = False, unspecified_to_ref: bool = False,
               reference_regions_to_retain: Path = None,
               concurrent_mode: bool = False, max_processes: int = 1, verbose: int = 0) -> List[Path]:
    """
    Normalize and prepare the input VCF for PharmCAT.
    """
    if len(vcf_files) == 0:
        raise ReportableException('Missing VCF input')

    # make sure we have samples
    if samples is None or len(samples) == 0:
        samples = util.read_vcf_samples(vcf_files[0], verbose=verbose)
    else:
        # make sure samples are in vcf file
        vcf_samples = util.read_vcf_samples(vcf_files[0], verbose=verbose)
        for sample in samples:
            if sample not in vcf_samples:
                raise ReportableException('Sample "%s" not in VCF' % sample)

    basename = output_basename or input_basename
    multisample_vcf = _preprocess(pharmcat_positions_vcf, reference_genome,
                                  vcf_files, samples,
                                  output_dir, basename,
                                  keep_intermediate_files, absent_to_ref, unspecified_to_ref,
                                  reference_regions_to_retain,
                                  concurrent_mode, max_processes, verbose)
    if split_samples and len(samples) > 1:
        util.index_vcf(multisample_vcf, verbose)
        # output PharmCAT-ready single-sample VCF
        # retain only the PharmCAT allele defining positions in the output VCF file
        results = util.export_single_sample_vcf(multisample_vcf, samples, output_dir, basename,
                                                concurrent_mode=concurrent_mode, max_processes=max_processes)
        if not keep_intermediate_files:
            util.delete_vcf_and_index(multisample_vcf, verbose=verbose)
        return results

    else:
        final_vcf: Path = finalize_multisample_vcf(multisample_vcf, output_dir, basename)
        return [final_vcf]


def preprocess_multiple_files(pharmcat_positions_vcf: Path, reference_genome: Path,
                              vcf_files: List[Path], samples: Optional[List[str]],
                              output_dir: Path, output_basename: Optional[str] = '',
                              keep_intermediate_files: bool = False,
                              absent_to_ref: bool = False, unspecified_to_ref: bool = False,
                              reference_regions_to_retain: Path = None,
                              concurrent_mode: bool = False, max_processes: int = 1, verbose: int = 0) -> List[Path]:
    """
    Normalize and prepare the input VCF for PharmCAT.
    """
    if len(vcf_files) == 0:
        return []

    results: List[Path] = []
    for file in vcf_files:
        # make sure we have samples
        file_samples: List[str] = []
        if samples is None or len(samples) == 0:
            file_samples = util.read_vcf_samples(file, verbose=verbose)
        else:
            # make sure samples are in the VCF file
            vcf_samples = util.read_vcf_samples(vcf_files[0], verbose=verbose)
            for sample in samples:
                if sample in vcf_samples:
                    file_samples.append(sample)
        if len(file_samples) == 0:
            continue

        basename = output_basename or util.get_vcf_basename(file)
        multisample_vcf = _preprocess(pharmcat_positions_vcf, reference_genome,
                                      [file], file_samples,
                                      output_dir, basename,
                                      keep_intermediate_files,
                                      absent_to_ref, unspecified_to_ref,
                                      reference_regions_to_retain,
                                      concurrent_mode, max_processes, verbose)
        results.append(finalize_multisample_vcf(multisample_vcf, output_dir, basename))

    return results


def _preprocess(pharmcat_positions_vcf: Path, reference_genome: Path,
                vcf_files: List[Path], samples: Optional[List[str]],
                output_dir: Path, output_basename: Optional[str] = '',
                keep_intermediate_files: bool = False,
                absent_to_ref: bool = False, unspecified_to_ref: bool = False,
                reference_regions_to_retain: Path = None,
                concurrent_mode=False, max_processes=1, verbose: int = 0) -> Path:

    # shrink input VCF down to PGx allele defining regions and selected samples
    # standardize chromosome names to <chr##>
    pgx_region_vcf: Path = util.extract_pgx_regions(pharmcat_positions_vcf, vcf_files, samples,
                                                    output_dir, output_basename,
                                                    reference_regions_to_retain,
                                                    concurrent_mode=concurrent_mode, max_processes=max_processes,
                                                    verbose=verbose)
    # normalize the input VCF
    normalized_vcf = util.normalize_vcf(reference_genome, pgx_region_vcf, output_dir, output_basename, verbose=verbose)

    # extract the specific PGx genetic variants in the reference PGx VCF
    # this step also generates a report of missing PGx positions in the input VCF
    pgx_variants_vcf: Path = util.extract_pgx_variants(pharmcat_positions_vcf, reference_genome, normalized_vcf,
                                                       output_dir, output_basename,
                                                       absent_to_ref=absent_to_ref,
                                                       unspecified_to_ref=unspecified_to_ref,
                                                       retain_specific_regions=bool(reference_regions_to_retain),
                                                       verbose=verbose)

    if not keep_intermediate_files:
        util.delete_vcf_and_index(pgx_region_vcf, verbose=verbose)
        util.delete_vcf_and_index(normalized_vcf, verbose=verbose)

    return pgx_variants_vcf


def finalize_multisample_vcf(file: Path, output_dir: Path, basename: str) -> Path:
    final_vcf: Path = output_dir / (basename + '.preprocessed.vcf.bgz')
    shutil.move(file, final_vcf)
    return final_vcf
