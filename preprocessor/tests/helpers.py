import hashlib
import shutil
import urllib.request
from pathlib import Path
from typing import Optional, List

import preprocessor
from preprocessor import download_reference_fasta_and_index


TEST_DOWNLOAD = False
NETWORK_AVAILABLE = True

try:
    with urllib.request.urlopen('https://google.com') as response:
        NETWORK_AVAILABLE = True
except:
    NETWORK_AVAILABLE = False


test_dir: Path = Path(globals().get("__file__", "./_")).absolute().parent
src_dir: Path = test_dir / '../preprocessor'
pharmcat_positions_file: Path = test_dir / '../../pharmcat_positions.vcf.bgz'


def get_reference_fasta(pharmcat_positions: Path) -> Path:
    reference_fasta: Path = pharmcat_positions.parent / preprocessor.REFERENCE_FASTA_FILENAME
    if not reference_fasta.is_file():
        if TEST_DOWNLOAD:
            download_reference_fasta_and_index(pharmcat_positions.parent, True)
        else:
            raise RuntimeError("CANNOT TEST: no reference fasta and TEST_DOWNLOAD=False")
    # also make sure uniallelic positions vcf exists
    uniallelic_positions_vcf: Path = preprocessor.find_uniallelic_file(pharmcat_positions, must_exist=False)
    if not uniallelic_positions_vcf.is_file():
        preprocessor.create_uniallelic_vcf(uniallelic_positions_vcf, pharmcat_positions, reference_fasta)
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


def read_vcf(file: Path, skip_comments: bool = True):
    """Reads VCF file and (1) strips trailing spaces, (2) removes empty lines and (3) normalizes line endings."""
    with open(file, 'r') as f:
        lines = []
        for line in f:
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
        preprocessor.run(['gunzip', str(actual_gz)])
        actual: Path = tmp_dir / ('%s.preprocessed.vcf' % basename)
    assert actual.is_file(), '%s not found' % actual
    if results is not None:
        assert orig_actual_vcf in results, ('%s not in results (%s)' % (actual, results))
    if copy_to_test_dir:
        shutil.copyfile(actual, actual.parent / actual.name)
    print(read_vcf(actual))

    assert read_vcf(expected) == read_vcf(actual), '%s mismatch' % key
