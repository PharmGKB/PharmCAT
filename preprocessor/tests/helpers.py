import gzip
import hashlib
import os
import shutil
import urllib.request
from pathlib import Path
from typing import Optional, List

import pcat
from pcat import download_reference_fasta_and_index


TEST_DOWNLOAD = False
NETWORK_AVAILABLE = True

# noinspection PyBroadException
try:
    with urllib.request.urlopen('https://google.com') as response:
        NETWORK_AVAILABLE = True
except:
    NETWORK_AVAILABLE = False

if os.environ.get('PHARMCAT_TEST_DOWNLOAD'):
    TEST_DOWNLOAD = True

test_dir: Path = Path(globals().get("__file__", "./_")).absolute().parent
src_dir: Path = test_dir / '../pcat'
pharmcat_positions_file: Path = test_dir / '../../pharmcat_positions.vcf.bgz'
uniallelic_pharmcat_positions_file: Path = test_dir / '../../pharmcat_positions.uniallelic.vcf.bgz'


def get_reference_fasta(pharmcat_positions: Path) -> Path:
    reference_fasta: Path = pharmcat_positions.parent / pcat.REFERENCE_FASTA_FILENAME
    if not reference_fasta.is_file():
        if TEST_DOWNLOAD:
            download_reference_fasta_and_index(pharmcat_positions.parent, True)
        else:
            raise RuntimeError("CANNOT TEST: no reference fasta and TEST_DOWNLOAD=False")
    # also make sure uniallelic positions vcf exists
    uniallelic_positions_vcf: Path = pcat.find_uniallelic_file(pharmcat_positions, must_exist=False)
    if not uniallelic_positions_vcf.is_file():
        pcat.create_uniallelic_vcf(uniallelic_positions_vcf, pharmcat_positions, reference_fasta)
    return reference_fasta


def touch(file: Path, text: Optional[str] = None):
    with open(file, 'a') as f:
        if text is not None:
            f.write(text)


def md5hash(file: Path):
    with open(file, 'rb') as f:
        file_hash = hashlib.md5()
        while chunk := f.read(8192):
            file_hash.update(chunk)
        return file_hash.hexdigest()


def read_vcf(file: Path, bgzipped: bool = False, skip_comments: bool = True):
    if bgzipped:
        with gzip.open(file, mode='rt', encoding='utf-8') as in_f:
            return _read_vcf(in_f, skip_comments=skip_comments)
    else:
        with open(file, mode='r', encoding='utf-8') as in_f:
            return _read_vcf(in_f, skip_comments=skip_comments)


def _read_vcf(in_f, skip_comments: bool = True):
    """Reads the VCF file and (1) strips trailing spaces, (2) removes empty lines and (3) normalizes line endings."""
    lines = []
    for line in in_f:
        line = line.rstrip()
        if line.startswith('##') and skip_comments:
            continue
        if line:
            lines.append(line)
    return '\n'.join(lines)


def compare_vcf_files(expected: Path, tmp_dir: Path, basename: str, sample: str = None, split_sample: bool = False,
                      copy_to_test_dir: bool = False, results: Optional[List[Path]] = None):
    key: str = basename
    orig_actual_vcf: Path
    if sample:
        key += '/%s' % sample
        actual: Path = tmp_dir / ('%s.%s.preprocessed.vcf' % (basename, sample))
        orig_actual_vcf = actual
    elif split_sample:
        actual: Path = tmp_dir / ('%s.preprocessed.vcf' % basename)
        orig_actual_vcf = actual
    else:
        actual_bgz: Path = tmp_dir / ('%s.preprocessed.vcf.bgz' % basename)
        assert actual_bgz.is_file(), '%s not found' % actual_bgz
        orig_actual_vcf = actual_bgz
        actual_gz: Path = tmp_dir / ('%s.preprocessed.vcf.gz' % basename)
        # use shutil.move instead of rename to deal with cross-device issues
        shutil.move(actual_bgz, actual_gz)
        pcat.run(['gunzip', str(actual_gz)])
        actual: Path = tmp_dir / ('%s.preprocessed.vcf' % basename)
    assert actual.is_file(), '%s not found' % actual
    if results is not None:
        assert orig_actual_vcf in results, ('%s not in results (%s)' % (actual, results))
    if copy_to_test_dir:
        shutil.copyfile(actual, actual.parent / actual.name)

    # compare vcfs line by line
    expected_lines = read_vcf(expected).split('\n')
    actual_lines = read_vcf(actual).split('\n')

    if len(expected_lines) != len(actual_lines):
        assert False, f'Different number of lines (expected {len(expected_lines)}, found {len(actual_lines)})'

    line_num = 0
    for expected_line, actual_line in zip(expected_lines, actual_lines):
        line_num += 1
        columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'Samples']
        expected_fields = expected_line.split('\t')
        actual_fields = actual_line.split('\t')

        if len(expected_fields) != len(actual_fields):
            assert False, f'Line {line_num}: different number of samples (expected {len(expected_fields)}, found {len(actual_fields)})'


        # print(f'EXPECTED: {expected_line}')
        # print(f'  ACTUAL: {actual_line}')
        for i, col in enumerate(columns):
            # compare the ALT alleles
            if i == 4 and set(actual_fields[3]) != set(expected_fields[3]):
                assert False, f'Line {line_num}: mismatched {col}\nexpected: {expected_fields[3]}\n  actual: {actual_fields[3]}'
            # compare genotypes
            if i == 9 and actual_fields[9:] != expected_fields[9:]:
                assert False, f'Line {line_num}: mismatched {col}\nexpected: {expected_fields[9]}\n  actual: {actual_fields[9]}'
            # compare the rest
            if actual_line[i] != expected_line[i]:
                assert False, f'Line {line_num}: mismatched {col}\nexpected: {expected_fields[i]}\n  actual: {actual_fields[i]}'
