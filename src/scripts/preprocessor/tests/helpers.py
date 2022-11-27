import hashlib
import re
import shutil
import tempfile
from pathlib import Path
from typing import Optional, List

import preprocessor
from preprocessor import download_reference_fasta_and_index, ReportableException, bgzip_vcf


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


def compare_vcf_files(expected: Path, tmp_dir: Path, basename: str, sample: str = None, copy_to_test_dir: bool = False):
    if sample:
        actual: Path = tmp_dir / ('%s.%s.preprocessed.vcf' % (basename, sample))
    else:
        actual: Path = tmp_dir / ('%s.preprocessed.vcf' % basename)
    if copy_to_test_dir:
        shutil.copyfile(actual, actual.parent / actual.name)
    assert actual.is_file(), '%s not found' % actual
    assert read_vcf(expected) == read_vcf(actual), '%s mismatch' % sample


def split_test_vcf():
    test_vcf = test_dir / 'test.vcf.bgz'
    test1_vcf = test_dir / 'test1.vcf.bgz'
    test2_vcf = test_dir / 'test2.vcf.bgz'
    with tempfile.TemporaryDirectory() as td:
        tmp_dir: Path = Path(td)
        tmp_vcf_gz = tmp_dir / 'test.vcf.gz'
        shutil.copyfile(test_vcf, tmp_vcf_gz)

        preprocessor.run(['gunzip', str(tmp_vcf_gz)])
        tmp_vcf = tmp_dir / 'test.vcf'
        assert tmp_vcf.is_file()

        comments: List[str] = []
        data1: List[str] = []
        data2: List[str] = []
        pattern = re.compile(r'(?:chr)?([XYM]|\d+)\s.*')
        with open(tmp_vcf, 'r') as f:
            for line in f:
                line = line.rstrip()
                if line[0] == '#':
                    comments.append(line)
                else:
                    rez = pattern.match(line)
                    if not rez:
                        raise ReportableException('Unrecognized chromosome in: %s' % line)
                    if rez.group(1).isnumeric() and int(rez.group(1)) <= 10:
                        data1.append(line)
                    else:
                        data2.append(line)

        tmp_vcf1 = tmp_dir / 'test1.vcf'
        with open(tmp_vcf1, 'w') as f:
            f.write('\n'.join(comments))
            f.write('\n')
            f.write('\n'.join(data1))
            f.write('\n')
        tmp_bgz1 = bgzip_vcf(tmp_vcf1)

        tmp_vcf2 = tmp_dir / 'test2.vcf'
        with open(tmp_vcf2, 'w') as f:
            f.write('\n'.join(comments))
            f.write('\n')
            f.write('\n'.join(data2))
            f.write('\n')
        tmp_bgz2 = bgzip_vcf(tmp_vcf2)

        shutil.copyfile(tmp_bgz1, test1_vcf)
        shutil.copyfile(tmp_bgz2, test2_vcf)


if __name__ == "__main__":
    split_test_vcf()

