import hashlib
import shutil
from pathlib import Path
from typing import Optional

import preprocessor
from preprocessor import download_reference_fasta_and_index


TEST_DOWNLOAD = False

test_dir: Path = Path(globals().get("__file__", "./_")).absolute().parent
src_dir: Path = test_dir / '../preprocessor'
pharmcat_positions_file: Path = src_dir / '../../../../pharmcat_positions.vcf.bgz'


def get_reference_fasta(pharmcat_positions: Path) -> Path:
    reference_fasta: Path = pharmcat_positions.parent / preprocessor.REFERENCE_FASTA_FILENAME
    if not reference_fasta.is_file():
        if TEST_DOWNLOAD:
            download_reference_fasta_and_index(pharmcat_positions.parent, True)
        else:
            raise RuntimeError("CANNOT TEST: no reference fasta and TEST_DOWNLOAD=False")
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


def read_vcf(file: Path):
    """Reads VCF file and (1) strips trailing spaces, (2) removes empty lines and (3) normalizes line endings."""
    with open(file, 'r') as f:
        lines = []
        for line in f:
            line = line.rstrip()
            if line:
                lines.append(line)
        return '\n'.join(lines)


def compare_vcf_files(expected: Path, tmp_dir: Path, basename: str, sample: str, copy_to_test_dir: bool = False):
    actual: Path = tmp_dir / ('%s.%s.preprocessed.vcf' % (basename, sample))
    if copy_to_test_dir:
        shutil.copyfile(actual, actual.parent / actual.name)
    assert actual.is_file(), '%s not found' % actual
    assert read_vcf(expected) == read_vcf(actual), '%s mismatch' % sample


