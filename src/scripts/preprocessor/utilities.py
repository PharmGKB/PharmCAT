#! /usr/bin/env python
__author__ = 'BinglanLi'

import concurrent.futures
import copy
import gzip
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import textwrap
import urllib.parse
import urllib.request
from concurrent.futures import ALL_COMPLETED
from pathlib import Path
from typing import Optional, Union, List

import allel
import pandas as pd
from packaging import version

from exceptions import ReportableException, InappropriateVCFSuffix, InvalidURL


# version is dynamically updated - DO NOT MODIFY MANUALLY
PHARMCAT_VERSION = '2.1.1'

# expected tool versions
MIN_BCFTOOLS_VERSION = '1.16'
MIN_BGZIP_VERSION = '1.16'

UNIALLELIC_VCF_SUFFIX = '.uniallelic.vcf.bgz'
UNIALLELIC_VCF_FILENAME = 'pharmcat_positions' + UNIALLELIC_VCF_SUFFIX

this = sys.modules[__name__]
this.script_dir = Path(globals().get("__file__", "./_")).absolute().parent
this.chr_rename_file = this.script_dir / 'chr_rename_map.tsv'
this.bcftools_path = 'bcftools'
this.bgzip_path = 'bgzip'

# chromosome names
_chr_invalid = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
                "18", "19", "20", "21", "22", "X", "Y", "M", "MT", "chrMT"]
_chr_valid = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
              "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
              "chr22", "chrX", "chrY", "chrM", "chrM", "chrM"]
_chr_valid_sorter = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
                     "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
                     "chr22", "chrX", "chrY", "chrM"]


def find_uniallelic_file(pharmcat_positions: Path, must_exist: bool = True) -> Path:
    uniallelic_positions_vcf: Path = pharmcat_positions.parent / UNIALLELIC_VCF_FILENAME
    if must_exist and not uniallelic_positions_vcf.is_file():
        raise ReportableException('Cannot find %s' % UNIALLELIC_VCF_FILENAME)
    return uniallelic_positions_vcf


def run(command: List[str]):
    try:
        subprocess.run(command, check=True, stderr=subprocess.PIPE, universal_newlines=True)
    except FileNotFoundError:
        raise ReportableException('Error: %s not found' % command[0])
    except subprocess.TimeoutExpired:
        raise ReportableException('Error: %s took too long' % command[0])
    except subprocess.CalledProcessError as e:
        if e.stderr:
            raise ReportableException(e.stderr)
        else:
            raise ReportableException('Error: Failed to run %s' % ' '.join(command))


def validate_tool(tool_name: str, tool_path: str, min_version: Optional[str] = None):
    """
    Validates that tool is available and meets minimum version requirement.
    Version checking only works for samtools tools, which has standardized version info in its help text.

    :param tool_name: name of tool (used for error messages)
    :param tool_path: path to tool
    :param min_version: minimum version of tool
    :raises ReportableException if tool cannot be found or does not meet version requirement
    """
    try:
        help_message = subprocess.run([tool_path, '-h'], stdout=subprocess.PIPE, check=True, stderr=subprocess.PIPE,
                                      universal_newlines=True).stdout
    except FileNotFoundError:
        raise ReportableException('Error: %s not found' % tool_path)
    except subprocess.TimeoutExpired:
        raise ReportableException('Error: %s took too long' % tool_path)
    except subprocess.CalledProcessError as e:
        raise ReportableException(e.stderr if e.stderr else 'Error: Failed to run %s' % tool_path)

    if min_version is not None:
        # check that the minimum version requirement is met
        rez = re.search(r'Version: (\d+(\.\d+)*)', str(help_message), re.MULTILINE)
        if rez is not None:
            tool_version = rez.group(1)
            if version.parse(tool_version) < version.parse(min_version):
                raise ReportableException("Error: Please use %s %s or higher." % (tool_name, min_version))
        else:
            raise ReportableException(textwrap.dedent("""
            Error: Could not find the version information for %s.
            Please use %s %s or higher.
            """ % (tool_name, tool_name, min_version)))


def validate_bcftools(tool_path: Optional[str], min_version: Optional[str] = None):
    """
    Validates that bcftools is available and meets minimum version requirement.

    :param tool_path: path to bcftools
    :param min_version: minimum required bcftools version
    :raises ReportableException if bcftools cannot be found or does not meet version requirement
    """
    tool_path = tool_path if tool_path else 'bcftools'
    min_version = min_version if min_version else MIN_BCFTOOLS_VERSION
    validate_tool('bcftools', tool_path, min_version)
    this.bcftools_path = tool_path


def validate_bgzip(tool_path: Optional[str], min_version: Optional[str] = None):
    """
    Validates that bgzip is available and meets minimum version requirement.

    :param tool_path: path to bgzip
    :param min_version: minimum required bgzip version
    :raises ReportableException: if bgzip cannot be found or does not meet version requirement
    """
    bgzip_path = tool_path if tool_path else 'bgzip'
    bgzip_version = min_version if min_version else MIN_BGZIP_VERSION
    validate_tool('bgzip', bgzip_path, bgzip_version)
    this.bgzip_path = bgzip_path


def validate_file(file: Union[Path, str]) -> Path:
    """
    Checks that the specified file exists.

    :raises ReportableException: if it does not
    """
    if isinstance(file, str):
        file = Path(file)
    if not file.is_file():
        if file.exists():
            raise ReportableException("Error: %s is not a file" % file)
        else:
            raise ReportableException("Error: %s does not exist" % file)
    return file


def find_vcf_files(vcf_dir: Path, verbose: bool = False) -> List[Path]:
    """
    Finds all VCF files in the specified directory.
    VCF file can be compressed (either .bgz or .gz extension).
    If file exists in both compressed and uncompressed form, pick the compressed form.

    :raises ReportableException: if no VCF files can be found
    """
    if not vcf_dir.is_dir():
        raise ReportableException('%s is not a directory' % vcf_dir)
    if verbose:
        print('Looking for VCF files in', vcf_dir)
    vcf_dict: dict[str, List[str]] = {}
    for f in vcf_dir.iterdir():
        if f.is_file():
            if is_vcf_file(f.name):
                vcf_basename = get_vcf_basename(f.name)
                if vcf_basename in vcf_dict:
                    vcf_dict[vcf_basename].append(f.name)
                else:
                    vcf_dict[vcf_basename] = [f.name]
    vcf_files: List[Path] = []
    for basename in vcf_dict.keys():
        if basename.startswith('pharmcat_positions'):
            if verbose:
                print('  - ignoring:', ', '.join(map(str, vcf_dict[basename])))
            continue
        vcfs: List[str] = vcf_dict[basename]
        if len(vcfs) == 1:
            vcf_files.append(vcf_dir / vcfs[0])
        elif (basename + '.vcf.bgz') in vcfs:
            vcf_files.append(vcf_dir / (basename + '.vcf.bgz'))
        elif (basename + '.vcf.gz') in vcfs:
            vcf_files.append(vcf_dir / (basename + '.vcf.gz'))
        elif (basename + '.vcf') in vcfs:
            vcf_files.append(vcf_dir / (basename + '.vcf'))
    if len(vcf_files) == 0:
        raise ReportableException('Error: no VCF files found')
    return vcf_files


def find_file(filename: str, dirs: List[Path]) -> Optional[Path]:
    """
    Looks for filename in specified directories.
    """
    file: Optional[Path] = None
    for d in dirs:
        f = d / filename
        if f.is_file():
            file = f
            break
    return file


def is_vcf_file(file: Union[Path, str]):
    return re.search('\\.vcf(\\.b?gz)?$', str(file)) is not None


def get_vcf_basename(path: Union[Path, str]) -> str:
    """
    Gets the base filename of a VCF file (can be compressed).

    :raises InappropriateVCFSuffix: if file does not have expected file extension
    """
    if not isinstance(path, Path):
        path = Path(path)
    match = re.search('(.+?)((\\.pgx_regions)|(\\.normalized))*\\.vcf(\\.b?gz)?$', path.name)
    if not match:
        raise InappropriateVCFSuffix(path)
    return match.group(1)


def validate_samples(samples: List[str]):
    if any(',' in sample_name for sample_name in samples):
        raise ReportableException("Error: Please remove comma ',' from sample names "
                                  "(commas violate bcftools sample name convention).")


def read_sample_file(sample_file: Path) -> List[str]:
    print('Reading samples from', sample_file, '...')
    samples: List[str] = []
    with open(sample_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line and not line.startswith("#"):
                samples.append(line)
    if len(samples) == 0:
        print('  * WARNING: No samples found. Will use all samples listed in VCF.')
    else:
        validate_samples(samples)
    return samples


def read_vcf_samples(vcf_file: Path) -> List[str]:
    """
    Obtain a list of samples from the input VCF.

    "bcftools query <options> <input_vcf>". For bcftools common options, see running_bcftools().
    "-l" list sample names and exit.
    Samples are delimited by '\\\\n' and the last line ends as 'last_sample_ID\\\\n\\\\n'.
    """

    print('Reading samples from', vcf_file, '...')
    output = subprocess.check_output([this.bcftools_path, 'query', '-l', str(vcf_file)], universal_newlines=True)
    vcf_sample_list: List[str] = output.split('\n')[:-1]  # remove the black line at the end
    if len(vcf_sample_list) == 0:
        raise ReportableException('Error: No samples found in VCF.')
    validate_samples(vcf_sample_list)
    return vcf_sample_list


def is_gz_file(file: Path):
    """
    Checks whether a file is bgzip compressed.
    """
    with open(file, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


def bgzip_file(file: Path, verbose: bool = False):
    """
    bgzip the specified file.
    Will overwrite existing .gz/.bgz file.
    """
    if verbose:
        print('  * Bgzipping', file)
    run([this.bgzip_path, '-f', str(file)])
    # make sure it worked
    gz_path: Path = Path(str(file) + '.gz')
    if not gz_path.exists():
        raise ReportableException('Cannot find (b)gzipped file %s' % gz_path)
    # rename the file to use .bgz extension
    bgz_path: Path = gz_path.with_suffix('.bgz')
    # must delete existing because won't overwrite on Windows
    bgz_path.unlink(missing_ok=True)
    # use shutil.move instead of rename to deal with cross-device issues
    shutil.move(gz_path, bgz_path)
    return bgz_path


def bgzip_vcf(file: Path, verbose: bool = False) -> Path:
    """
    Make sure file is bgzipped.
    Will overwrite existing .gz/.bgz file.
    Will delete pre-existing .bgz indices.
    """
    if is_gz_file(file):
        return file
    file = bgzip_file(file, verbose)
    delete_index(file, '.tbi', verbose=verbose)
    delete_index(file, '.csi', verbose=verbose)
    return file


def delete_index(vcf_file: Path, suffix: str, verbose: bool = False):
    index = Path(str(vcf_file) + suffix)
    if index.exists():
        if index.is_file():
            if verbose:
                print('  * Removing pre-existing %s index' % suffix)
            index.unlink()
        else:
            raise ReportableException('Error: Cannot delete obsolete index.  %s is not a file' % index)


def find_index_file(vcf_file: Path) -> Optional[Path]:
    csi_file = Path(str(vcf_file) + '.csi')
    if csi_file.is_file():
        return csi_file
    tbi_file = Path(str(vcf_file) + '.tbi')
    if tbi_file.is_file():
        return tbi_file
    return None


def index_vcf(vcf_file: Path, verbose: bool = False) -> Path:
    """
    Index the input vcf using bcftools, and the output index file will be written to the working directory.
    """
    if verbose:
        print('  * Generating index for %s' % vcf_file)
    run([this.bcftools_path, 'index', str(vcf_file)])
    csi_file = Path(str(vcf_file) + '.csi')
    if not csi_file.exists():
        raise ReportableException('Cannot find indexed .csi file %s' % csi_file)
    return csi_file


def delete_vcf_and_index(vcf_file: Path, verbose: bool = False):
    """Delete compressed vcf as well as the index file."""
    if vcf_file.is_file():
        vcf_file.unlink()
    delete_index(vcf_file, '.csi', verbose)


def download_from_url(url: str, download_dir: Path, force_update: bool = False,
                      verbose: bool = False):
    """Download from a URL."""

    remote_basename = os.path.basename(urllib.parse.urlparse(url).path)
    if remote_basename:
        dl_file = download_dir / remote_basename
        if dl_file.exists():
            if force_update:
                dl_file.unlink()
            else:
                return dl_file

        if verbose:
            print('  Downloading from "%s"\n\tto "%s"' % (url, dl_file))
        with urllib.request.urlopen(url) as response:
            with open(dl_file, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
        return dl_file
    else:
        raise InvalidURL(url)


def download_reference_fasta_and_index(download_dir: Union[Path, str], force_update: bool = False,
                                       verbose: bool = False):
    """Download the human reference genome sequence GRCh38/hg38."""
    if not isinstance(download_dir, Path):
        download_dir = Path(download_dir)
    ref_file = download_dir / 'reference.fna.bgz'
    if ref_file.exists() and not force_update:
        return ref_file

    tar_file = download_from_url('https://zenodo.org/record/7288118/files/GRCh38_reference_fasta.tar',
                                 download_dir, force_update, verbose)
    with tarfile.open(tar_file, 'r') as tar:
        tar.extractall(path=download_dir)
    os.remove(tar_file)

    return ref_file


def prep_pharmcat_positions(pharmcat_positions_vcf: Optional[Path] = None,
                            reference_genome_fasta: Optional[Path] = None,
                            update_chr_rename_file: bool = False, verbose: bool = False):
    if pharmcat_positions_vcf is None:
        # assume it's in current working directory
        pharmcat_positions_vcf = Path('pharmcat_positions.vcf.bgz')
    if not pharmcat_positions_vcf.is_file():
        raise ReportableException('Cannot find %s' % pharmcat_positions_vcf)

    if reference_genome_fasta is None:
        # assume it's in current working directory
        reference_genome_fasta = Path('reference.fna.bgz')
    if not reference_genome_fasta.is_file():
        download_reference_fasta_and_index(pharmcat_positions_vcf.parent, verbose)

    csi_file = Path(str(pharmcat_positions_vcf) + '.csi')
    if not csi_file.is_file():
        index_vcf(pharmcat_positions_vcf, verbose)

    uniallelic_positions_vcf: Path = find_uniallelic_file(pharmcat_positions_vcf, must_exist=False)
    if not uniallelic_positions_vcf.is_file():
        if verbose:
            print('* Preparing uniallelic PharmCAT positions VCF')
        # convert reference PGx variants to the uniallelic format needed for extracting exact PGx positions
        # and generating an accurate missing report
        bcftools_command = [this.bcftools_path, 'norm', '--no-version', '-m-', '-c', 'ws', '-f',
                            str(reference_genome_fasta), '-Oz', '-o', str(uniallelic_positions_vcf),
                            str(pharmcat_positions_vcf)]
        run(bcftools_command)
        index_vcf(uniallelic_positions_vcf, verbose)

    # create chromosome mapping file
    if not this.chr_rename_file.is_file() or update_chr_rename_file:
        if verbose:
            print('* Preparing chromosome rename file')
        with open(this.chr_rename_file, 'w+') as f:
            for i in range(len(_chr_invalid)):
                f.write(_chr_invalid[i] + "\t" + _chr_valid[i] + "\n")


def _get_vcf_pos_min_max(positions, flanking_bp=100):
    """ given input positions, return "<min_pos>-<max_pos>"  """
    return '-'.join([str(min(positions) - flanking_bp), str(max(positions) + flanking_bp)])


def _is_valid_chr(vcf_file: Path) -> bool:
    with gzip.open(vcf_file, mode='rt', encoding='utf-8') as in_f:
        for line in in_f:
            if line[0] != '#':
                fields = line.rstrip().split()
                if fields[0] in _chr_valid:
                    return True
                elif fields[0] in _chr_invalid:
                    return False
                else:
                    break
    raise ReportableException('The CHROM column does not conform with either "chr##" or "##" format.')


def extract_pgx_regions(pharmcat_positions: Path, vcf_files: List[Path], samples: List[str],
                        output_dir: Path, output_basename: str,
                        concurrent_mode: bool = False, max_processes: int = 1,
                        verbose: bool = False) -> Path:
    """
    Extracts PGx regions from input VCF file(s) into a single VCF file and rename chromosomes to match PharmCAT
    expectations.

    :return: compressed VCF file that only contains positions of interest
    """

    with tempfile.TemporaryDirectory() as td:
        tmp_dir: Path = Path(td)
        # generate sample list
        tmp_sample_file: Path = tmp_dir / 'samples.txt'
        with open(tmp_sample_file, 'w+') as w:
            for sample in samples:
                w.write('%s\n' % sample)

        pgx_region_vcf_file: Path = output_dir / (output_basename + '.pgx_regions.vcf.bgz')
        if len(vcf_files) == 1:
            # this should create pgx_region_vcf_file
            _extract_pgx_regions(pharmcat_positions, vcf_files[0], tmp_sample_file, output_dir, output_basename,
                                 verbose)
        else:
            # generate files to be concatenated
            tmp_files: List[Path] = []
            if concurrent_mode:
                with concurrent.futures.ProcessPoolExecutor(max_workers=check_max_processes(max_processes)) as e:
                    futures = []
                    for vcf_file in vcf_files:
                        futures.append(e.submit(_extract_pgx_regions, pharmcat_positions, vcf_file, tmp_sample_file,
                                                output_dir, get_vcf_basename(vcf_file), verbose))
                    concurrent.futures.wait(futures, return_when=ALL_COMPLETED)
                    for future in futures:
                        tmp_files.append(future.result())
            else:
                for vcf_file in vcf_files:
                    tmp_files.append(_extract_pgx_regions(pharmcat_positions, vcf_file, tmp_sample_file, output_dir,
                                                          get_vcf_basename(vcf_file), verbose))
            # write file names to txt file for bcftools
            tmp_file_list = tmp_dir / 'regions.txt'
            with open(tmp_file_list, 'w+') as w:
                for tf in tmp_files:
                    w.write(str(tf) + "\n")

            # concatenate vcfs
            if verbose:
                print('Concatenating PGx VCFs')
            bcftools_command = [this.bcftools_path, 'concat', '--no-version', '-a', '-f', str(tmp_file_list), '-Oz',
                                '-o', str(pgx_region_vcf_file)]
            run(bcftools_command)
            # index the VCF file
            index_vcf(pgx_region_vcf_file, verbose)

        return pgx_region_vcf_file


def _extract_pgx_regions(pharmcat_positions: Path, vcf_file: Path, sample_file: Path, output_dir: Path,
                         output_basename: Optional[str], verbose: bool = False) -> Path:
    """
    Does the actual work to extract PGx regions from input VCF file(s) into a single VCF file and
    rename chromosomes to match PharmCAT expectations.

    "bcftools annotate <options> <vcf_file>".
    "--rename-chrs" renames chromosomes according to the map in chr_rename_map.tsv.
    """
    print('Processing', vcf_file, '...')
    # make sure vcf is bgzipped and indexed
    bgz_file = bgzip_vcf(vcf_file, verbose)
    idx_file = find_index_file(bgz_file)
    if idx_file is None:
        index_vcf(bgz_file, verbose)

    # obtain PGx regions to be extracted
    df_ref_pgx_pos = allel.vcf_to_dataframe(str(pharmcat_positions))
    df_ref_pgx_pos['CHROM'] = df_ref_pgx_pos['CHROM'].astype('category').cat.set_categories(_chr_valid_sorter)
    ref_pgx_regions = df_ref_pgx_pos.groupby(['CHROM'], observed=True)['POS'].agg(_get_vcf_pos_min_max).reset_index()
    # add a special case for 'chrMT'
    idx_chr_m = ref_pgx_regions.index[ref_pgx_regions['CHROM'] == 'chrM']
    ref_pgx_regions = pd.concat([ref_pgx_regions, ref_pgx_regions.loc[idx_chr_m].assign(**{'CHROM': 'chrMT'})])

    # validate chromosome formats
    if _is_valid_chr(bgz_file):
        # format the pgx regions to be extracted
        ref_pgx_regions = ",".join(ref_pgx_regions.apply(lambda row: ':'.join(row.values.astype(str)), axis=1))
    else:
        # format pgx regions
        ref_pgx_regions = ",".join(ref_pgx_regions.apply(lambda row: ':'.join(row.values.astype(str)), axis=1)
                                   .replace({'chr': ''}, regex=True))

    # extract pgx regions and modify chromosome names if necessary
    if output_basename is None:
        output_basename = get_vcf_basename(vcf_file)
    pgx_regions_vcf = output_dir / (output_basename + '.pgx_regions.vcf.bgz')
    bcftools_command = [this.bcftools_path, 'annotate', '--no-version', '-S', str(sample_file),
                        '--rename-chrs', str(this.chr_rename_file), '-r', ref_pgx_regions, '-i', 'ALT="."', '-k',
                        '-Oz', '-o', str(pgx_regions_vcf), str(bgz_file)]
    if verbose:
        print('  * Extracting PGx regions and normalizing chromosome names')
    run(bcftools_command)
    # index the PGx VCF file
    index_vcf(pgx_regions_vcf, verbose)
    return pgx_regions_vcf


def normalize_vcf(reference_genome: Path, vcf_file: Path, output_dir: Path, output_basename: Optional[str],
                  verbose: bool = False):
    """
    Normalize the input VCF against the human reference genome sequence GRCh38/hg38.

    "bcftools norm <options> <vcf_file>". For bcftools common options, see running_bcftools().
    "-m +|-" joins bi-allelic sites into multi-allelic records (+)
        and convert multi-allelic records into uniallelic format (-).
    "-f <reference_genome_fasta>" reference sequence. Supplying this option turns on left-alignment and normalization.
    "-c ws" when incorrect or missing REF allele is encountered, warn (w) and set/fix(s) bad sites.  's' will swap
    alleles and update GT and AC acounts. Importantly, s will NOT fix strand issues in a VCF.
    """
    if output_basename is None:
        output_basename = get_vcf_basename(vcf_file)
    normalized_vcf = output_dir / (output_basename + '.normalized.vcf.bgz')
    bcftools_command = [this.bcftools_path, 'norm', '--no-version', '-m-', '-c', 'ws', '-Oz', '-o', str(normalized_vcf),
                        '-f', str(reference_genome), str(vcf_file)]
    if verbose:
        print('Normalizing VCF')
    run(bcftools_command)
    index_vcf(normalized_vcf, verbose)
    return normalized_vcf


def _is_phased(gt_field) -> bool:
    """
    Determines the phasing status of a position.
    If the GT fields can be split by '/', this means at least one sample is unphased.
    """
    for x in gt_field:
        if len(x.split('/')) > 1:
            return False
        else:
            return True


def extract_pgx_variants(pharmcat_positions: Path, reference_fasta: Path, vcf_file: Path,
                         output_dir: Path, output_basename: str, missing_to_ref: bool = False,
                         verbose: bool = False) -> Path:
    """
    Extract specific pgx positions that are present in the reference PGx VCF.
    Generate a report of PGx positions that are missing in the input VCF.

    "bcftools isec <options> <vcf_file>". For bcftools common options, see running_bcftools().
    "-c none" only records with the same CHR, POS, REF and ALT are considered identical
    "-C" outputs positions present only in the first file but missing in the others.
    "-n=2" positions shared by both inputs
    "-i 'TYPE=snp' " to include only SNPs
    "-w" lists input files to output given as 1-based indices.
        "-w1" extracts and writes records only present in the first file (the reference PGx positions).
    """

    with tempfile.TemporaryDirectory() as td:
        tmp_dir = Path(td)
        uniallelic_positions_vcf: Path = find_uniallelic_file(pharmcat_positions)

        # extracting PGx positions from input using "bcftools view"
        pgx_pos_only_bgz: Path = tmp_dir / (output_basename + '.pgx_pos_only.vcf.bgz')
        if verbose:
            print('  * Retaining PGx positions, regardless of alleles')
        run([this.bcftools_path, 'view', '--no-version', '-U', '-T', str(uniallelic_positions_vcf), '-Oz',
             '-o', str(pgx_pos_only_bgz), str(vcf_file)])
        index_vcf(pgx_pos_only_bgz, verbose)

        # add in missing multi-allelic variants as '0|0'
        # bcftools will take care of the rest (sorting, normalization, etc.)
        #
        # first create a dictionary of PharmCAT reference PGx positions
        # the ref_pos_static (dict) will be used to static dictionary of reference PGx positions
        # the ref_pos_dynamic (dict) will be used to retain only PGx pos in the input
        ref_pos_dynamic = {}
        with gzip.open(uniallelic_positions_vcf, mode='rt', encoding='utf-8') as in_f:
            for line in in_f:
                # skip headers
                if line[0] == '#':
                    continue
                # read file
                line = line.rstrip('\n')
                fields = line.split('\t')
                # ref_pos_dynamic: a nested dictionary
                # ref_pos_dynamic[(chr,pos)][(ref, alt)] = all fields except GT
                # initiate the dictionary if the first key pair doesn't exist yet
                if (fields[0], fields[1]) not in ref_pos_dynamic.keys():
                    ref_pos_dynamic[(fields[0], fields[1])] = {}
                ref_pos_dynamic[(fields[0], fields[1])][(fields[3], fields[4])] = fields[0:9]
        # deep copy the ref_pos_dynamic to have a consistent reference of PharmCAT PGx positions
        ref_pos_static = copy.deepcopy(ref_pos_dynamic)

        # use input_pos to record PGx positions that are present in the input VCF
        input_pos = []
        input_pos_phased = {}
        # this dictionary saves genetic variants concurrent at PGx positions
        dict_non_pgx_records = {}
        print('Updating VCF header and PGx position annotations')
        updated_pgx_pos_vcf: Path = tmp_dir / (output_basename + '.updated_pgx_pos_only.vcf')
        with open(updated_pgx_pos_vcf, 'w') as out_f:
            # get header of samples from merged vcf, add in new contig info
            with gzip.open(pgx_pos_only_bgz, mode='rt', encoding='utf-8') as in_f:
                for line in in_f:
                    # process header lines, skip contig
                    if line[0:8] == '##contig':
                        continue  # skip contig info in the original vcf
                    # process header lines, except lines about contig
                    elif line[0:2] == '##':
                        out_f.write(line)
                    # process the sample line, add PGx-related lines
                    elif line[0:6] == '#CHROM':
                        for single_chr in _chr_valid_sorter:
                            out_f.write('##contig=<ID=' + single_chr + ',assembly=GRCh38.p13,species="Homo sapiens">\n')
                        out_f.write('##INFO=<ID=PX,Number=.,Type=String,Description="Gene">\n')
                        out_f.write('##INFO=<ID=POI,Number=0,Type=Flag,Description="Position of Interest but not'
                                    ' part of an allele definition">\n')
                        out_f.write('##FILTER=<ID=PCATxREF,Description="Reference allele does not match PharmCAT '
                                    'reference alleles">\n')
                        out_f.write('##FILTER=<ID=PCATxALT,Description="Alternate alleles do not match PharmCAT '
                                    'alternate alleles">\n')
                        out_f.write('##FILTER=<ID=PCATxINDEL,Description="Unexpected format for INDELs">\n')
                        out_f.write(line)
                        # get the number of samples
                        line = line.rstrip('\n')
                        fields: List[str] = line.split('\t')
                        n_sample = len(fields) - 9
                    # scan the genotype data
                    else:
                        line = line.rstrip('\n')
                        fields = line.split('\t')
                        input_chr_pos = (fields[0], fields[1])

                        # check whether input is a block gVCF, which will be supported in the future
                        if re.search('END', fields[7]):
                            print('=============================================================\n'
                                  'The PharmCAT VCF Preprocessor will support block gVCF in the future.\n'
                                  '=============================================================\n')
                            sys.exit(1)

                        # match chromosome positions
                        if input_chr_pos in ref_pos_static:
                            '''
                            match REF and ALT alleles
                            
                            SNPs:
                                1. Matching REF and ALT: retain and update the rsID and FORMAT/PGx gene name
                                2. Matching REF and ALT=. or <*>: homozygous reference SNPs, retain, 
                                    update rsID and FORMAT/PGx gene name
                                3. Matching REF and mismatching ALT
                                    3.1. There is another matching ALT: retain it in the same line
                                    3.2. There is no other matching ALT: label "FILTER/PCATxALT"
                            INDELs        
                                1. Matching REF and ALT: retain and update the rsID and FORMAT/PGx gene name
                                2. Matching REF and mismatching ALT: retain as non-PGx record
                                3. Matching REF and ALT=<*>: uncertain nucleotide changes, warn and ignore
                            '''

                            # check whether the position has unspecified alt '<*>' or is homozygous reference
                            if fields[4] in ['.', '<*>']:
                                is_nonspecific_alt: bool = True
                            else:
                                is_nonspecific_alt: bool = False

                            # list out REF alleles at a position
                            ref_alleles = [x[0] for x in ref_pos_static[input_chr_pos].keys()]
                            alt_alleles = [x[1] for x in ref_pos_static[input_chr_pos].keys()]
                            # check if the reference PGx variant is a SNP
                            len_ref: list[int] = [len(i) for i in ref_alleles]
                            len_alt: list[int] = [len(i) for i in alt_alleles]
                            is_snp: bool = (len_ref == len_alt)

                            # a tuple variable for ref and alt alleles in this VCF record
                            input_ref_alt = (fields[3], fields[4])

                            # update id when the id is not in the input
                            ref_id = list(set([x[2] for x in ref_pos_static[input_chr_pos].values()]))
                            input_id = fields[2].split(';')
                            updated_id = ';'.join(set(ref_id + input_id)) if fields[2] != '.' else ';'.join(ref_id)

                            # update info
                            ref_info = list(set([x[7] for x in ref_pos_static[input_chr_pos].values()]))
                            if fields[7] != '.':
                                ref_info.append(fields[7])
                            updated_info = ';'.join(ref_info)

                            # positions with matching REF and ALT
                            if input_ref_alt in ref_pos_static[input_chr_pos]:
                                # record input PGx positions in a list
                                # use this to fill up multi-allelic ALT later
                                if input_chr_pos not in input_pos:
                                    input_pos.append(input_chr_pos)
                                    # determine the phasing status of a position
                                    # the (position is unphased if any sample is unphased)
                                    input_pos_phased[input_chr_pos] = _is_phased(fields[9:])

                                # update id
                                fields[2] = updated_id
                                # update info
                                fields[7] = updated_info

                                # concat and write to output file
                                out_f.write('\t'.join(fields) + '\n')
                                # elimination: remove the dictionary item so that the variant won't be matched again
                                if input_chr_pos in ref_pos_dynamic:
                                    ref_pos_dynamic[input_chr_pos].pop(input_ref_alt)
                                    # remove a position if all of its alts are present in the input
                                    if ref_pos_dynamic[input_chr_pos] == {}:
                                        del ref_pos_dynamic[input_chr_pos]
                            # SNPs - homozygous reference
                            elif is_snp and is_nonspecific_alt and (fields[3] in ref_alleles):
                                # if reference ALT alleles are exhausted, next positions
                                if input_chr_pos not in ref_pos_dynamic:
                                    continue
                                else:
                                    # fill up all ALT alleles info
                                    alt_alleles = [x[1] for x in ref_pos_dynamic[input_chr_pos].keys()]
                                    for i in range(len(alt_alleles)):
                                        # update contents, concatenate cols from input and from the ref PGx VCF
                                        fields[2] = updated_id
                                        fields[3] = ref_alleles[i]
                                        fields[4] = alt_alleles[i]
                                        fields[7] = updated_info
                                        # concat and write to output file
                                        out_f.write('\t'.join(fields) + '\n')

                                        # for hom ref SNPs, remove the position from the dict for record
                                        if input_chr_pos in ref_pos_dynamic:
                                            ref_pos_dynamic[input_chr_pos].pop((ref_alleles[i], alt_alleles[i]))
                                            # remove a position if all of its alts are present in the input
                                            if ref_pos_dynamic[input_chr_pos] == {}:
                                                del ref_pos_dynamic[input_chr_pos]
                            # SNPs - mismatching ALT
                            elif is_snp and (fields[3] in ref_alleles):
                                # record input PGx positions in a list
                                # use this to fill up multi-allelic ALT later
                                if input_chr_pos not in input_pos:
                                    input_pos.append(input_chr_pos)
                                    input_pos_phased[input_chr_pos] = _is_phased(fields[9:])
                                # update info
                                fields[7] = updated_info
                                # write to output
                                out_f.write('\t'.join(fields) + '\n')
                            # INDELs with a unspecified allele, ALT="<*>"
                            elif not is_snp and is_nonspecific_alt:
                                print('  * WARNING: ignore \"%s:%s REF=%s ALT=%s\" which is not a valid GT format '
                                      'for INDELs'
                                      % (fields[0], fields[1], fields[3], fields[4]))
                                # update filter
                                fields[6] = 'PCATxINDEL'
                                # save the line in the dictionary for non-PGx variants
                                dict_non_pgx_records[input_chr_pos] = '\t'.join(fields)
                            # flag if the variant doesn't match PharmCAT ALT
                            elif fields[3] in ref_alleles:
                                for i in range(len(ref_alleles)):
                                    print('  * WARNING: \"%s:%s REF=%s ALT=%s\" does not match PharmCAT expectation '
                                          'of ALT at "%s:%s REF=%s ALT=%s"'
                                          % (fields[0], fields[1], fields[3], fields[4],
                                             fields[0], fields[1], ref_alleles[i], alt_alleles[i]))
                                # update filter
                                fields[6] = 'PCATxALT'
                                # save the line in the dictionary for non-PGx variants
                                dict_non_pgx_records[input_chr_pos] = '\t'.join(fields)
                            # flag if the variant doesn't match PharmCAT REF
                            else:
                                for i in range(len(ref_alleles)):
                                    print('  * WARNING: \"%s:%s REF=%s ALT=%s\" does not match PharmCAT expectation '
                                          'of REF at "%s:%s REF=%s ALT=%s"'
                                          % (fields[0], fields[1], fields[3], fields[4],
                                             fields[0], fields[1], ref_alleles[i], alt_alleles[i]))
                                # update filter
                                fields[6] = 'PCATxREF'
                                # save the line in the dictionary for non-PGx variants
                                dict_non_pgx_records[input_chr_pos] = '\t'.join(fields)
            # complete lines of multi-allelic loci or missing positions
            for key_chr_pos in ref_pos_dynamic:
                for val in ref_pos_dynamic[key_chr_pos].values():
                    # complete multi-allelic loci
                    if key_chr_pos in input_pos:
                        if input_pos_phased[key_chr_pos]:
                            line = '\t'.join(val + ['0|0'] * n_sample)
                        else:
                            line = '\t'.join(val + ['0/0'] * n_sample)
                        out_f.write(line + '\n')
                    # if missing_to_ref is true, output lines of missing positions as homozygous reference
                    elif missing_to_ref:
                        line = '\t'.join(val + ['0|0'] * n_sample)
                        out_f.write(line + '\n')
                    else:
                        continue

        # sort vcf
        print('Sorting by chromosomal location...')
        sorted_bgz: Path = tmp_dir / (output_basename + '.sorted.vcf.bgz')
        bcftools_command = [this.bcftools_path, 'sort', '-Oz', '-o', str(sorted_bgz), str(updated_pgx_pos_vcf)]
        run(bcftools_command)
        index_vcf(sorted_bgz, verbose)

        # make sure output complies with the multi-allelic format
        print('Enforcing multi-allelic variant representation...')
        normed_bgz: Path = tmp_dir / (output_basename + '.normed.vcf.bgz')
        run([this.bcftools_path, 'norm', '--no-version', '-m+', '-c', 'ws', '-f', str(reference_fasta), '-Oz',
             '-o', str(normed_bgz), str(sorted_bgz)])

        filtered_bgz: Path = output_dir / (output_basename + '.multiallelic.vcf.bgz')

        # sort non-PGx variants according to genomic positions
        pos_list = ref_pos_static.keys()
        non_pgx_records = [dict_non_pgx_records[key] for key in pos_list if key in dict_non_pgx_records]
        # if there are non-PGx variant at PGx positions, need to put back the non-PGx variants into VCF
        if len(non_pgx_records) >= 1:
            # insert lines of concurrent non-PGx variants after the PGx positions
            filtered_vcf: Path = output_dir / (output_basename + '.multiallelic.vcf')
            with open(filtered_vcf, 'w') as out_f:
                print('Adding back non-PGx variants at PGx positions...')
                with gzip.open(normed_bgz, mode='rt', encoding='utf-8') as in_f:
                    for line in in_f:
                        # print all header lines
                        if line[0] == '#':
                            out_f.write(line)
                        # print the rest if there is no more non-PGx variants
                        elif len(non_pgx_records) == 0:
                            out_f.write(line)
                        # scan the genotype data
                        else:
                            # put back non-PGx variants before the next PGx position
                            line = line.rstrip('\n')
                            fields = line.split('\t')
                            while len(non_pgx_records):
                                non_pgx_fields = non_pgx_records[0].split('\t')
                                if _chr_valid_sorter.index(fields[0]) > _chr_valid_sorter.index(non_pgx_fields[0]):
                                    out_f.write(non_pgx_records[0] + '\n')
                                    non_pgx_records.pop(0)
                                elif _chr_valid_sorter.index(fields[0]) == _chr_valid_sorter.index(non_pgx_fields[0]) \
                                        and int(fields[1]) > int(non_pgx_fields[1]):
                                    out_f.write(non_pgx_records[0] + '\n')
                                    non_pgx_records.pop(0)
                                else:
                                    break
                            # output original line
                            out_f.write(line + '\n')
                    # print the rest of the non-PGx variants
                    while len(non_pgx_records):
                        out_f.write(non_pgx_records[0] + '\n')
                        non_pgx_records.pop(0)
            # bgzip
            bgzip_vcf(filtered_vcf, verbose)
        else:
            # use shutil.move instead of rename to deal with cross-device issues
            shutil.move(normed_bgz, filtered_bgz)

        # report missing positions in the input VCF
        _print_missing_positions(pharmcat_positions, ref_pos_dynamic, output_dir, output_basename, verbose)

        index_vcf(filtered_bgz, verbose)
        return filtered_bgz


def _print_missing_positions(pharmcat_positions: Path, ref_pos_dynamic, output_dir: Path, output_basename: str,
                             verbose: bool = False):
    # report missing positions in the input VCF
    print('Generating report of missing PGx allele defining positions...')
    missing_pos_file: Path = output_dir / (output_basename + '.missing_pgx_var.vcf')
    with open(missing_pos_file, 'w') as out_f:
        # get VCF header from pharmcat_positions
        with gzip.open(pharmcat_positions, mode='rt', encoding='utf-8') as in_f:
            for line in in_f:
                # process header lines
                if line[0:2] == '##':
                    out_f.write(line)
                # process the sample line, add PGx-related lines
                elif line[0:6] == '#CHROM':
                    line = line.rstrip('\n')
                    fields: List[str] = line.split('\t')
                    fields[-1] = 'Missing'
                    out_f.write('\t'.join(fields) + '\n')
                else:
                    break
        # output positions that were not detected in the input VCF
        for input_ref_alt in ref_pos_dynamic.values():
            for val in input_ref_alt.values():
                line = '\t'.join(val + ['0|0'])
                out_f.write(line + '\n')
    # bgzip the missing report
    bgzip_vcf(missing_pos_file, verbose)


def _output_pharmcat_ready_vcf(vcf_file: Path, output_dir: Path, output_basename: str, single_sample: str):
    """
    Create final PharmCAT-ready VCF file.

    "bcftools view <options> <input_vcf>". For bcftools common options, see running_bcftools().
    "-U" exclude sites without any called genotype, i.e., all GT = './.'
    "--force-samples" only warn about unknown subset samples
    """
    if output_basename:
        output_file_name = output_dir / (output_basename + '.' + single_sample + '.preprocessed.vcf')
    else:
        output_file_name = output_dir / (single_sample + '.preprocessed.vcf')
    print('Generating PharmCAT-ready VCF for', single_sample)
    run([this.bcftools_path, 'view', '--no-version', '--force-samples', '-U', '-s', single_sample, '-Ov',
         '-o', str(output_file_name), str(vcf_file)])


def output_pharmcat_ready_vcf(vcf_file: Path, samples: List[str], output_dir: Path, output_basename: str,
                              concurrent_mode: bool = False, max_processes: int = 1):
    """
    Write final PharmCAT-ready VCF for each sample.
    """
    if concurrent_mode:
        with concurrent.futures.ProcessPoolExecutor(max_workers=check_max_processes(max_processes)) as e:
            futures = []
            for single_sample in samples:
                futures.append(e.submit(_output_pharmcat_ready_vcf, vcf_file, output_dir, output_basename,
                                        single_sample))
            concurrent.futures.wait(futures, return_when=ALL_COMPLETED)
    else:
        for single_sample in samples:
            _output_pharmcat_ready_vcf(vcf_file, output_dir, output_basename, single_sample)


def check_max_processes(max_processes: int) -> int:
    if max_processes is None:
        # number of cpus - 1, minimum of 1
        max_processes = max(1, os.cpu_count() - 1)
    elif max_processes < 1:
        max_processes = 1
    if os.name == 'nt' and max_processes > 61:
        # windows has max of 61 workers
        max_processes = 61
    return max_processes
