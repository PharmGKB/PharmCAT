import hashlib
import os
import shutil
import tempfile
from pathlib import Path
from timeit import default_timer as timer
from typing import Optional, List
from unittest import TestCase

import utilities as utils
from exceptions import ReportableException, InappropriateVCFSuffix


script_dir: Path = Path(globals().get("__file__", "./_")).absolute().parent
TEST_DOWNLOAD = False


def get_reference_fasta(pharmcat_positions: Path) -> Path:
    reference_fasta: Path = pharmcat_positions.parent / 'reference.fna.bgz'
    if not reference_fasta.is_file():
        if TEST_DOWNLOAD:
            utils.download_reference_fasta_and_index(pharmcat_positions.parent, True)
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


def compare_vcf_files(self, expected: Path, tmp_dir: Path, basename: str, sample: str, copy_to_test_dir: bool = False):
    actual: Path = tmp_dir / ('%s.%s.preprocessed.vcf' % (basename, sample))
    if copy_to_test_dir:
        shutil.copyfile(actual, actual.parent / actual.name)
    self.assertTrue(actual.is_file(), '%s not found' % actual)
    self.assertEqual(read_vcf(expected), read_vcf(actual), ('%s mismatch' % sample))


class Test(TestCase):
    def test_chr_arrays(self):
        # make sure chr arrays are of the same length
        self.assertTrue(len(utils._chr_invalid) == len(utils._chr_valid))

    def test_validate_tool(self):
        utils.validate_tool('bcftools', 'bcftools')

    def test_validate_tool_fail(self):
        with self.assertRaises(ReportableException) as context:
            utils.validate_tool('bcftools', 'bcftools', '99')
        print(context.exception)
        self.assertIn('use bcftools 99 or higher', context.exception.msg)

    def test_validate_bcftools(self):
        utils.validate_bcftools(None)
        self.assertEqual('bcftools', utils.bcftools_path)

        utils.validate_bcftools(None, '1.16')
        self.assertEqual('bcftools', utils.bcftools_path)

        with self.assertRaises(ReportableException) as context:
            utils.validate_bcftools('foo/bcftools', '1.16')
        print(context.exception)
        self.assertIn('not found', context.exception.msg)
        self.assertEqual('bcftools', utils.bcftools_path)

        with self.assertRaises(ReportableException) as context:
            utils.validate_bcftools(None, '99')
        print(context.exception)

    def test_validate_bgzip(self):
        utils.validate_bgzip(None)
        self.assertEqual('bgzip', utils.bgzip_path)

        utils.validate_bgzip(None, '1.16')
        self.assertEqual('bgzip', utils.bgzip_path)

        with self.assertRaises(ReportableException) as context:
            utils.validate_bgzip('foo/bgzip')
        print(context.exception)
        self.assertIn('not found', context.exception.msg)
        self.assertEqual('bgzip', utils.bgzip_path)

        with self.assertRaises(ReportableException) as context:
            utils.validate_bgzip(None, '99')
        print(context.exception)

    def test_find_vcf_files(self):
        vcf_files = utils.find_vcf_files(script_dir / 'test')
        self.assertTrue(len(vcf_files) > 1)

    def test_find_vcf_files_fail(self):
        with self.assertRaises(ReportableException) as context:
            utils.find_vcf_files(script_dir)
        self.assertIn('no VCF files found', context.exception.msg)

        with self.assertRaises(ReportableException) as context:
            utils.find_vcf_files(script_dir / 'bad-dir')
        print(context.exception)
        self.assertIn('not a directory', context.exception.msg)

    def test_is_vcf_file(self):
        valid_paths = [
            Path('/this/dir/file.vcf'),
            Path('/this/dir/file.vcf.bgz'),
            Path('/this/dir/file.vcf.gz')
        ]
        for p in valid_paths:
            self.assertTrue(utils.is_vcf_file(p), msg=str(p))
            self.assertTrue(utils.is_vcf_file(p.name), msg=str(p))

        invalid_paths = [
            Path('/this/dir/'),
            Path('/this/dir'),
            Path('/this/dir/file.txt'),
            Path('/this/dir/file.vcf.zip'),
        ]
        for p in invalid_paths:
            self.assertFalse(utils.is_vcf_file(p), msg=str(p))
            self.assertFalse(utils.is_vcf_file(p.name), msg=str(p))

    def test_get_vcf_basename(self):
        valid_paths = [
            Path('/this/dir/file.vcf'),
            Path('/this/dir/file.vcf.bgz'),
            Path('/this/dir/file.vcf.gz'),
            Path('/this/dir/file.pgx_regions.vcf.gz'),
            Path('/this/dir/file.normalized.vcf.gz'),
            Path('/this/dir/file.pgx_regions.normalized.vcf.gz')
        ]
        for p in valid_paths:
            self.assertEqual('file', utils.get_vcf_basename(p), msg=str(p))
            self.assertEqual('file', utils.get_vcf_basename(p.name), msg=str(p))

        invalid_path = Path('/this/dir/file.txt')
        with self.assertRaises(InappropriateVCFSuffix) as context:
            utils.get_vcf_basename(invalid_path)
        print(context.exception)

        with self.assertRaises(InappropriateVCFSuffix) as context:
            utils.get_vcf_basename(invalid_path.name)
        print(context.exception)

    def test_read_sample_file(self):
        samples = utils.read_sample_file(script_dir / 'test/test-empty-samples.txt')
        self.assertIsNotNone(samples)
        self.assertEqual(0, len(samples))

        samples = utils.read_sample_file(script_dir / 'test/test-samples.txt')
        self.assertTrue(samples is not None)
        self.assertEqual(2, len(samples))
        self.assertCountEqual(['SAMPLE_1', 'SAMPLE_2'], samples)

    def test_read_vcf_samples(self):
        samples = utils.read_vcf_samples(script_dir / 'test/reference.Sample_1.preprocessed.vcf')
        self.assertIsNotNone(samples)
        print(samples)
        self.assertEqual(1, len(samples))
        self.assertIn('Sample_1', samples)

        samples = utils.read_vcf_samples(script_dir / 'test/test.vcf.bgz')
        self.assertIsNotNone(samples)
        print(samples)
        self.assertEqual(2, len(samples))
        self.assertCountEqual(['Sample_1', 'Sample_2'], samples)

    def test_is_gz_file(self):
        self.assertTrue(utils.is_gz_file(script_dir / 'test/test.vcf.bgz'))
        self.assertFalse(utils.is_gz_file(script_dir / 'test/reference.Sample_1.preprocessed.vcf'))

    def test_bgzip_file(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_file = Path(tmp_dir, 'foo.txt')
            touch(tmp_file, 'hello, world')
            bgz_file = utils.bgzip_file(tmp_file, True)
            self.assertTrue(utils.is_gz_file(bgz_file))

    def test_bgzip_file_fail(self):
        with self.assertRaises(ReportableException) as context:
            utils.bgzip_file(Path('/no-such-file'))
        print(context.exception)
        self.assertIn('No such file', context.exception.msg)

    def test_bgzip_vcf(self):
        vcf_file = Path(script_dir / 'test/reference.Sample_1.preprocessed.vcf')
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_file = Path(tmp_dir, 's1.vcf')
            shutil.copyfile(vcf_file, tmp_file)
            csi_file = Path(tmp_dir, 's1.vcf.bgz.csi')
            touch(csi_file)
            tbi_file = Path(tmp_dir, 's1.vcf.bgz.tbi')
            touch(tbi_file)
            bgz_file = utils.bgzip_vcf(tmp_file, True)
            self.assertTrue(utils.is_gz_file(bgz_file))
            files = os.listdir(tmp_dir)
            self.assertNotIn('s1.vcf.bgz.csi', files)
            self.assertNotIn('s1.vcf.bgz.tbi', files)

    def test_index_vcf(self):
        vcf_file = Path(script_dir / 'test/reference.Sample_1.preprocessed.vcf')
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_file = Path(tmp_dir, 's1.vcf')
            shutil.copyfile(vcf_file, tmp_file)
            bgz_file = utils.bgzip_vcf(tmp_file, True)
            utils.index_vcf(bgz_file, True)

    def test_delete_vcf_and_index(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            vcf_file = Path(tmp_dir, 'foo.vcf')
            self.assertFalse(vcf_file.exists())
            utils.delete_vcf_and_index(vcf_file, True)

            csi_file = Path(str(vcf_file) + '.csi')
            touch(vcf_file)
            touch(csi_file)
            self.assertTrue(vcf_file.exists())
            self.assertTrue(csi_file.exists())
            utils.delete_vcf_and_index(vcf_file)
            self.assertFalse(vcf_file.exists())
            self.assertFalse(csi_file.exists())

    def test_download_reference_fasta_and_index(self):
        if TEST_DOWNLOAD:
            with tempfile.TemporaryDirectory() as tmp_dir:
                ref_file = utils.download_reference_fasta_and_index(tmp_dir, verbose=True)
                files = os.listdir(tmp_dir)
                self.assertIn(ref_file.name, files)

    def test_prep_pharmcat_positions(self):
        pharmcat_positions = Path(script_dir) / '../../../pharmcat_positions.vcf.bgz'
        with tempfile.TemporaryDirectory() as tmp_dir:
            os.chdir(tmp_dir)

            # no pharmcat_positions - should fail
            with self.assertRaises(ReportableException) as context:
                utils.prep_pharmcat_positions(None, None, verbose=True)
            # print(context.exception)
            self.assertIn('Cannot find pharmcat_positions.vcf.bgz', context.exception.msg)

            # add pharmcat positions to tmp_dir
            tmp_positions = Path(tmp_dir, pharmcat_positions.name)
            shutil.copyfile(pharmcat_positions, tmp_positions)
            tmp_uniallelic = Path(tmp_dir, 'pharmcat_positions.uniallelic.vcf.bgz')

            reference_fasta: Path
            if TEST_DOWNLOAD:
                # use pharmcat_positions from cwd and download reference fasta
                start = timer()
                utils.prep_pharmcat_positions(None, None, verbose=True)
                full_time = timer() - start
                print("time for full preparation:", full_time)
                reference_fasta = Path(tmp_dir, 'reference.fna.bgz')
                self.assertTrue(reference_fasta.is_file())

                self.assertTrue(tmp_uniallelic.is_file())
                tmp_uniallelic.unlink()
            else:
                # assumes that reference fasta is available next to pharmcat_positions
                reference_fasta = pharmcat_positions.parent / 'reference.fna.bgz'
                self.assertTrue(reference_fasta.is_file(), 'Cannot find reference FASTA for testing!')

            # already have reference fasta, so this should be faster, but still need to generate uniallelic positions
            start = timer()
            utils.prep_pharmcat_positions(tmp_positions, reference_fasta, verbose=True)
            build_time = timer() - start
            print("time to build uniallelic positions:", build_time)
            self.assertTrue(build_time < 1)

            tmp_uniallelic = Path(tmp_dir, 'pharmcat_positions.uniallelic.vcf.bgz')
            self.assertTrue(tmp_uniallelic.is_file())
            uniallelic_mtime = tmp_uniallelic.stat().st_mtime

            # already have reference fasta and uniallelic positions, so this should be very fast
            start = timer()
            # utils.prep_pharmcat_positions(None, tmp_reference, verbose=True)
            utils.prep_pharmcat_positions(tmp_positions, reference_fasta, verbose=True)
            check_time = timer() - start
            print("time to run check preparation:", check_time)
            self.assertTrue(check_time < build_time)
            # is an order of magnitude faster
            self.assertTrue(check_time < (build_time / 10))

            self.assertTrue(tmp_uniallelic.is_file())
            self.assertEqual(uniallelic_mtime, tmp_uniallelic.stat().st_mtime)

    def test_extract_pgx_regions(self):
        pharmcat_positions = Path(script_dir, '../../../pharmcat_positions.vcf.bgz')
        vcf_file = Path(script_dir, 'test/test.vcf.bgz')
        vcf_file1 = Path(script_dir, 'test/test1.vcf.bgz')
        vcf_file2 = Path(script_dir, 'test/test2.vcf.bgz')
        with tempfile.TemporaryDirectory() as td:
            tmp_dir = Path(td)
            tmp_positions = tmp_dir / pharmcat_positions.name
            shutil.copyfile(pharmcat_positions, tmp_positions)
            shutil.copyfile(utils.chr_rename_file, tmp_dir / utils.chr_rename_file.name)
            tmp_vcf = tmp_dir / vcf_file.name
            shutil.copyfile(vcf_file, tmp_vcf)
            tmp_vcf1 = tmp_dir / vcf_file1.name
            shutil.copyfile(vcf_file1, tmp_vcf1)
            tmp_vcf2 = tmp_dir / vcf_file2.name
            shutil.copyfile(vcf_file2, tmp_vcf2)

            vcf_files: List[Path] = [tmp_vcf1, tmp_vcf2]
            samples: List[str] = ['Sample_1', 'Sample_2']
            combo_pgx_vcf_file = utils.extract_pgx_regions(tmp_positions, vcf_files, samples, tmp_dir, 'combo_test',
                                                           verbose=True)
            self.assertTrue(combo_pgx_vcf_file.is_file())
            index_file = utils.find_index_file(combo_pgx_vcf_file)
            self.assertIsNotNone(index_file)
            combo_hash = md5hash(combo_pgx_vcf_file)
            print("combo:", combo_hash)

            pgx_vcf_file = utils.extract_pgx_regions(tmp_positions, [tmp_vcf], samples, tmp_dir, 'test', verbose=True)
            self.assertTrue(pgx_vcf_file.is_file())
            single_hash = md5hash(pgx_vcf_file)
            print("single:", single_hash)

            self.assertEqual(combo_hash, single_hash)

    def test_normalize_vcf(self):
        pharmcat_positions = Path(script_dir, '../../../pharmcat_positions.vcf.bgz')
        reference_fasta: Path = get_reference_fasta(pharmcat_positions)

        vcf_file = Path(script_dir, 'test/test.pgx_regions.vcf.bgz')
        with tempfile.TemporaryDirectory() as td:
            tmp_dir = Path(td)
            tmp_vcf = tmp_dir / vcf_file.name
            shutil.copyfile(vcf_file, tmp_vcf)

            normalized_vcf = utils.normalize_vcf(reference_fasta, tmp_vcf, tmp_dir, 'test', verbose=True)
            self.assertTrue(normalized_vcf.is_file())
            # shutil.copyfile(normalized_vcf, vcf_file.parent / 'test.normalized.vcf.bgz')

    def test_extract_pgx_variants(self):
        pharmcat_positions = Path(script_dir, '../../../pharmcat_positions.vcf.bgz')
        reference_fasta: Path = get_reference_fasta(pharmcat_positions)

        vcf_file = Path(script_dir, 'test/test.normalized.vcf.bgz')
        with tempfile.TemporaryDirectory() as td:
            tmp_dir = Path(td)
            tmp_vcf = tmp_dir / vcf_file.name
            shutil.copyfile(vcf_file, tmp_vcf)

            multiallelic_vcf = utils.extract_pgx_variants(pharmcat_positions, reference_fasta, tmp_vcf, tmp_dir, 'test',
                                                          verbose=True)
            self.assertTrue(multiallelic_vcf.is_file())
            # shutil.copyfile(multiallelic_vcf, vcf_file.parent / multiallelic_vcf.name)

    def test_output_pharmcat_ready_vcf(self):
        vcf_file = Path(script_dir, 'test/test.multiallelic.vcf.bgz')
        s1_file = script_dir / 'test/reference.Sample_1.preprocessed.vcf'
        s2_file = script_dir / 'test/reference.Sample_2.preprocessed.vcf'
        with tempfile.TemporaryDirectory() as td:
            tmp_dir = Path(td)
            tmp_vcf = tmp_dir / vcf_file.name
            shutil.copyfile(vcf_file, tmp_vcf)

            basename = 'test1'
            utils.output_pharmcat_ready_vcf(tmp_vcf, ['Sample_1', 'Sample_2'], tmp_dir, basename)

            compare_vcf_files(self, s1_file, tmp_dir, basename, 'Sample_1')
            compare_vcf_files(self, s2_file, tmp_dir, basename, 'Sample_2')

    def test_output_pharmcat_ready_vcf_partial(self):
        vcf_file = Path(script_dir, 'test/test.multiallelic.vcf.bgz')
        s1_file = script_dir / 'test/reference.Sample_1.preprocessed.vcf'
        with tempfile.TemporaryDirectory() as td:
            tmp_dir = Path(td)
            tmp_vcf = tmp_dir / vcf_file.name
            shutil.copyfile(vcf_file, tmp_vcf)

            basename = 'test_partial'
            utils.output_pharmcat_ready_vcf(tmp_vcf, ['Sample_1'], tmp_dir, basename)

            compare_vcf_files(self, s1_file, tmp_dir, basename, 'Sample_1')
            tmp_s2: Path = tmp_dir / ('%s.Sample_2.preprocessed.vcf' % basename)
            self.assertFalse(tmp_s2.is_file(), '%s found!' % tmp_s2)

    def test_output_pharmcat_ready_vcf_concurrent(self):
        vcf_file = Path(script_dir, 'test/test.multiallelic.vcf.bgz')
        s1_file = script_dir / 'test/reference.Sample_1.preprocessed.vcf'
        s2_file = script_dir / 'test/reference.Sample_2.preprocessed.vcf'
        with tempfile.TemporaryDirectory() as td:
            tmp_dir = Path(td)
            tmp_vcf = tmp_dir / vcf_file.name
            shutil.copyfile(vcf_file, tmp_vcf)

            basename = 'test2'
            utils.output_pharmcat_ready_vcf(tmp_vcf, ['Sample_1', 'Sample_2'], tmp_dir, basename, concurrent_mode=True,
                                            max_processes=2)

            compare_vcf_files(self, s1_file, tmp_dir, basename, 'Sample_1')
            compare_vcf_files(self, s2_file, tmp_dir, basename, 'Sample_2')
