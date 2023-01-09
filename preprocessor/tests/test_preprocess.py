import os
import shutil
import tempfile
from pathlib import Path

import helpers
from preprocessor import preprocess


def test_preprocess():
    reference_fasta: Path = helpers.get_reference_fasta(helpers.pharmcat_positions_file)

    vcf_file = helpers.test_dir / 'raw.vcf.bgz'
    preprocessed_file = helpers.test_dir / 'raw.preprocessed.vcf'
    with tempfile.TemporaryDirectory() as td:
        tmp_dir = Path(td)
        tmp_vcf = tmp_dir / vcf_file.name
        shutil.copyfile(vcf_file, tmp_vcf)

        basename = 'preprocess'
        preprocess(helpers.pharmcat_positions_file, reference_fasta, [tmp_vcf], None, basename, tmp_dir, basename,
                   verbose=True)
        print(tmp_dir)
        files = os.listdir(tmp_dir)
        print(files)
        helpers.compare_vcf_files(preprocessed_file, tmp_dir, basename)


def test_preprocess_split_sample():
    reference_fasta: Path = helpers.get_reference_fasta(helpers.pharmcat_positions_file)

    vcf_file = helpers.test_dir / 'raw.vcf.bgz'
    s1_file = helpers.test_dir / 'raw.Sample_1.preprocessed.vcf'
    s2_file = helpers.test_dir / 'raw.Sample_2.preprocessed.vcf'
    with tempfile.TemporaryDirectory() as td:
        tmp_dir = Path(td)
        tmp_vcf = tmp_dir / vcf_file.name
        shutil.copyfile(vcf_file, tmp_vcf)

        basename = 'preprocess'
        preprocess(helpers.pharmcat_positions_file, reference_fasta, [tmp_vcf], None, basename, tmp_dir, basename,
                   split_samples=True, verbose=True)

        helpers.compare_vcf_files(s1_file, tmp_dir, basename, 'Sample_1')
        helpers.compare_vcf_files(s2_file, tmp_dir, basename, 'Sample_2')


def test_preprocess_concurrent():
    reference_fasta: Path = helpers.get_reference_fasta(helpers.pharmcat_positions_file)

    vcf_file = helpers.test_dir / 'raw.vcf.bgz'
    preprocessed_file = helpers.test_dir / 'raw.preprocessed.vcf'
    with tempfile.TemporaryDirectory() as td:
        tmp_dir = Path(td)
        tmp_vcf = tmp_dir / vcf_file.name
        shutil.copyfile(vcf_file, tmp_vcf)

        basename = 'preprocess'
        preprocess(helpers.pharmcat_positions_file, reference_fasta, [tmp_vcf], None, basename, tmp_dir, basename,
                   concurrent_mode=True, max_processes=2, verbose=True)

        helpers.compare_vcf_files(preprocessed_file, tmp_dir, basename)


def test_preprocess_multi_vcf():
    reference_fasta: Path = helpers.get_reference_fasta(helpers.pharmcat_positions_file)

    vcf1_file = helpers.test_dir / 'raw-p1.vcf.bgz'
    vcf2_file = helpers.test_dir / 'raw-p2.vcf.bgz'
    preprocessed_file = helpers.test_dir / 'raw.preprocessed.vcf'
    with tempfile.TemporaryDirectory() as td:
        tmp_dir = Path(td)
        tmp_vcf1 = tmp_dir / vcf1_file.name
        shutil.copyfile(vcf1_file, tmp_vcf1)
        tmp_vcf2 = tmp_dir / vcf2_file.name
        shutil.copyfile(vcf2_file, tmp_vcf2)

        basename = 'preprocess'
        preprocess(helpers.pharmcat_positions_file, reference_fasta, [tmp_vcf1, tmp_vcf2], None, basename, tmp_dir,
                   basename)

        helpers.compare_vcf_files(preprocessed_file, tmp_dir, basename)


def test_preprocess_multi_vcf_concurrent():
    reference_fasta: Path = helpers.get_reference_fasta(helpers.pharmcat_positions_file)

    vcf1_file = helpers.test_dir / 'raw-p1.vcf.bgz'
    vcf2_file = helpers.test_dir / 'raw-p2.vcf.bgz'
    preprocessed_file = helpers.test_dir / 'raw.preprocessed.vcf'
    with tempfile.TemporaryDirectory() as td:
        tmp_dir = Path(td)
        tmp_vcf1 = tmp_dir / vcf1_file.name
        shutil.copyfile(vcf1_file, tmp_vcf1)
        tmp_vcf2 = tmp_dir / vcf2_file.name
        shutil.copyfile(vcf2_file, tmp_vcf2)

        basename = 'preprocess'
        preprocess(helpers.pharmcat_positions_file, reference_fasta, [tmp_vcf1, tmp_vcf2], None, basename, tmp_dir,
                   basename, concurrent_mode=True, max_processes=2)

        helpers.compare_vcf_files(preprocessed_file, tmp_dir, basename)
