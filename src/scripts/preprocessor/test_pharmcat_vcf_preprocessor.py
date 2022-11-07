import shutil
import tempfile
from pathlib import Path
from unittest import TestCase

from pharmcat_vcf_preprocessor import preprocess
from test_utilities import get_reference_fasta, compare_vcf_files


script_dir: Path = Path(globals().get("__file__", "./_")).absolute().parent


class Test(TestCase):
    def test_preprocess(self):
        pharmcat_positions = Path(script_dir, '../../../pharmcat_positions.vcf.bgz')
        reference_fasta: Path = get_reference_fasta(pharmcat_positions)

        vcf_file = Path(script_dir, 'test/test.vcf.bgz')
        s1_file = script_dir / 'test/reference.Sample_1.preprocessed.vcf'
        s2_file = script_dir / 'test/reference.Sample_2.preprocessed.vcf'
        with tempfile.TemporaryDirectory() as td:
            tmp_dir = Path(td)
            tmp_vcf = tmp_dir / vcf_file.name
            shutil.copyfile(vcf_file, tmp_vcf)

            basename = 'preprocess'
            preprocess(pharmcat_positions, reference_fasta, [tmp_vcf], None, basename, tmp_dir, basename, verbose=True)

            compare_vcf_files(self, s1_file, tmp_dir, basename, 'Sample_1')
            compare_vcf_files(self, s2_file, tmp_dir, basename, 'Sample_2')

    def test_preprocess_multi_vcf(self):
        pharmcat_positions = Path(script_dir, '../../../pharmcat_positions.vcf.bgz')
        reference_fasta: Path = get_reference_fasta(pharmcat_positions)

        vcf1_file = Path(script_dir, 'test/test1.vcf.bgz')
        vcf2_file = Path(script_dir, 'test/test2.vcf.bgz')
        s1_file = script_dir / 'test/reference.Sample_1.preprocessed.vcf'
        s2_file = script_dir / 'test/reference.Sample_2.preprocessed.vcf'
        with tempfile.TemporaryDirectory() as td:
            tmp_dir = Path(td)
            tmp_vcf1 = tmp_dir / vcf1_file.name
            shutil.copyfile(vcf1_file, tmp_vcf1)
            tmp_vcf2 = tmp_dir / vcf2_file.name
            shutil.copyfile(vcf2_file, tmp_vcf2)

            basename = 'preprocess'
            preprocess(pharmcat_positions, reference_fasta, [tmp_vcf1, tmp_vcf2], None, basename, tmp_dir, basename)

            compare_vcf_files(self, s1_file, tmp_dir, basename, 'Sample_1')
            compare_vcf_files(self, s2_file, tmp_dir, basename, 'Sample_2')
