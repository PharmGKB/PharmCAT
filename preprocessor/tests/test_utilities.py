import os
import shutil
import tempfile
from pathlib import Path
from timeit import default_timer as timer
from typing import List

import pytest

import helpers
import preprocessor
import preprocessor.utilities as utils
from preprocessor import ReportableException, InappropriateVCFSuffix, InvalidURL, common


def test_chr_arrays():
    # make sure chr arrays are of the same length
    assert len(utils._chr_invalid) == len(utils._chr_valid)


def test_validate_tool():
    utils.validate_tool('bcftools', 'bcftools')


def test_validate_tool_fail():
    with pytest.raises(ReportableException) as context:
        utils.validate_tool('bcftools', 'bcftools', '99')
    # print(context.value)
    assert 'use bcftools 99 or higher' in context.value.msg


def test_validate_bcftools():
    utils.validate_bcftools(None)
    assert 'bcftools' == preprocessor.BCFTOOLS_PATH

    utils.validate_bcftools(None, '1.18')
    assert 'bcftools' == preprocessor.BCFTOOLS_PATH

    with pytest.raises(ReportableException) as context:
        utils.validate_bcftools('foo/bcftools', '1.18')
    # print(context.value)
    assert 'not found' in context.value.msg
    assert 'bcftools' == preprocessor.BCFTOOLS_PATH

    with pytest.raises(ReportableException) as context:
        utils.validate_bcftools(None, '99')
    # print(context.value)
    assert '99 or higher' in context.value.msg


def test_validate_bgzip():
    utils.validate_bgzip(None)
    assert 'bgzip' == preprocessor.BGZIP_PATH

    utils.validate_bgzip(None, '1.18')
    assert 'bgzip' == preprocessor.BGZIP_PATH

    with pytest.raises(ReportableException) as context:
        utils.validate_bgzip('foo/bgzip')
    # print(context.value)
    assert 'not found' in context.value.msg
    assert 'bgzip' == preprocessor.BGZIP_PATH

    with pytest.raises(ReportableException) as context:
        utils.validate_bgzip(None, '99')
    # print(context.value)
    assert '99 or higher' in context.value.msg


def test_validate_java():
    utils.validate_java()
    assert 'java' == preprocessor.JAVA_PATH

    utils.validate_java('17')
    assert 'java' == preprocessor.JAVA_PATH

    with pytest.raises(ReportableException) as context:
        utils.validate_java('179')
    print(context.value)
    assert 'use Java 179 or higher' in context.value.msg


def test_validate_dir():
    with pytest.raises(ReportableException) as context:
        utils.validate_dir('/does/not/exist')
    assert 'does not exist' in context.value.msg

    with tempfile.TemporaryDirectory() as td:
        tmp_dir = Path(td)
        utils.validate_dir(td)
        utils.validate_dir(tmp_dir)

        with pytest.raises(ReportableException) as context:
            utils.validate_dir(td + '/does/not/exist')
        assert 'does not exist' in context.value.msg

        new_dir = utils.validate_dir(td + '/does/not/exist', create_if_not_exist=True)
        assert new_dir.exists()
        assert new_dir.is_dir()

    with tempfile.NamedTemporaryFile() as tf:
        with pytest.raises(ReportableException) as context:
            utils.validate_dir(tf.name)
        assert 'is not a directory' in context.value.msg


def test_validate_file():
    with pytest.raises(ReportableException) as context:
        utils.validate_file('/does/not/exist')
    assert 'does not exist' in context.value.msg

    with tempfile.NamedTemporaryFile() as tf:
        utils.validate_file(tf.name)

    with tempfile.TemporaryDirectory() as td:
        tmp_dir = Path(td)

        with pytest.raises(ReportableException) as context:
            utils.validate_file(tmp_dir)
        assert 'is not a file' in context.value.msg


def test_find_vcf_files():
    vcf_files = utils.find_vcf_files(helpers.test_dir)
    assert len(vcf_files) > 1


def test_find_vcf_files_fail():
    with pytest.raises(ReportableException) as context:
        utils.find_vcf_files(helpers.src_dir)
    assert 'no VCF files found' in context.value.msg

    with pytest.raises(ReportableException) as context:
        utils.find_vcf_files(helpers.src_dir / 'bad-dir')
    # print(context.value)
    assert 'not a directory' in context.value.msg


def test_is_vcf_file():
    valid_paths = [
        Path('/this/dir/file.vcf'),
        Path('/this/dir/file.vcf.bgz'),
        Path('/this/dir/file.vcf.gz')
    ]
    for p in valid_paths:
        assert utils.is_vcf_file(p), str(p)
        assert utils.is_vcf_file(p.name), str(p)

    invalid_paths = [
        Path('/this/dir/'),
        Path('/this/dir'),
        Path('/this/dir/file.txt'),
        Path('/this/dir/file.vcf.zip'),
    ]
    for p in invalid_paths:
        assert not utils.is_vcf_file(p), str(p)
        assert not utils.is_vcf_file(p.name), str(p)


def test_is_gvcf_file():
    test_vcf_file = helpers.test_dir / 'raw.vcf.bgz'
    test_gvcf_file = helpers.test_dir / 'test.g.vcf.bgz'
    with tempfile.TemporaryDirectory() as td:
        assert not utils.is_gvcf_file(test_vcf_file)

        tmp_dir = Path(td)
        f1 = tmp_dir / 'test.g.vcf.bgz'
        shutil.copyfile(test_vcf_file, f1)
        # is VCF file, but pass due to name
        assert utils.is_gvcf_file(f1)

        f2 = tmp_dir / 'test.genomic.vcf.gz'
        shutil.move(f1, f2)
        # is VCF file, but pass due to name
        assert utils.is_gvcf_file(f2)

        f1 = tmp_dir / 'test.genomic.vcf'
        shutil.move(f2, f1)
        # is VCF file, but pass due to name
        assert utils.is_gvcf_file(f1)

        f2 = tmp_dir / 'test.g.vcf'
        shutil.move(f1, f2)
        # is VCF file, but pass due to name
        assert utils.is_gvcf_file(f2)

        f1 = tmp_dir / 'test.vcf.bgz'
        shutil.copyfile(test_gvcf_file, f1)
        # pass due to compressed content check
        assert utils.is_gvcf_file(f1)

        f2 = tmp_dir / 'test.vcf.gz'
        shutil.move(f1, f2)
        utils.run(['gunzip', str(f2)])
        f = tmp_dir / 'test.vcf'
        assert utils.is_gvcf_file(f)


def test_get_vcf_basename():
    valid_paths = [
        Path('/this/dir/file.vcf'),
        Path('/this/dir/file.vcf.bgz'),
        Path('/this/dir/file.vcf.gz'),
        Path('/this/dir/file.pgx_regions.vcf.gz'),
        Path('/this/dir/file.normalized.vcf.gz'),
        Path('/this/dir/file.pgx_regions.normalized.vcf.gz')
    ]
    for p in valid_paths:
        assert 'file' == utils.get_vcf_basename(p), str(p)
        assert 'file' == utils.get_vcf_basename(p.name), str(p)

    invalid_path = Path('/this/dir/file.txt')
    with pytest.raises(InappropriateVCFSuffix) as context:
        utils.get_vcf_basename(invalid_path)
    # print(context.value)
    assert 'Inappropriate VCF suffix' in context.value.msg
    assert str(invalid_path) in context.value.msg

    with pytest.raises(InappropriateVCFSuffix) as context:
        utils.get_vcf_basename(invalid_path.name)
    # print(context.value)
    assert invalid_path.name in context.value.msg


def test_read_sample_file():
    samples = utils.read_sample_file(helpers.test_dir / 'test-empty-samples.txt')
    assert samples is not None
    assert 0 == len(samples)

    samples = utils.read_sample_file(helpers.test_dir / 'test-samples.txt')
    assert samples is not None
    assert 2 == len(samples)
    assert ['SAMPLE_1', 'SAMPLE_2'] == samples


def test_read_vcf_samples():
    samples = utils.read_vcf_samples(helpers.test_dir / 'raw.Sample_1.preprocessed.vcf')
    assert samples is not None
    # print(samples)
    assert 1 == len(samples)
    assert 'Sample_1' in samples

    samples = utils.read_vcf_samples(helpers.test_dir / 'raw.vcf.bgz')
    assert samples is not None
    # print(samples)
    assert 2 == len(samples)
    assert ['Sample_1', 'Sample_2'] == samples


def test_is_gz_file():
    assert utils.is_gz_file(helpers.test_dir / 'raw.vcf.bgz')
    assert not utils.is_gz_file(helpers.test_dir / 'raw.Sample_1.preprocessed.vcf')


def test_bgzip_file():
    with tempfile.TemporaryDirectory() as td:
        tmp_file = Path(td, 'foo.txt')
        helpers.touch(tmp_file, 'hello, world')
        bgz_file = utils.bgzip_file(tmp_file, True)
        assert utils.is_gz_file(bgz_file)


def test_bgzip_file_fail():
    with pytest.raises(ReportableException) as context:
        utils.bgzip_file(Path('/no-such-file'))
    # print(context.value)
    assert 'No such file' in context.value.msg


def test_bgzip_vcf():
    vcf_file = helpers.test_dir / 'raw.Sample_1.preprocessed.vcf'
    with tempfile.TemporaryDirectory() as td:
        tmp_dir: Path = Path(td)
        tmp_file = tmp_dir / 's1.vcf'
        shutil.copyfile(vcf_file, tmp_file)
        csi_file = tmp_dir / 's1.vcf.bgz.csi'
        helpers.touch(csi_file)
        tbi_file = tmp_dir / 's1.vcf.bgz.tbi'
        helpers.touch(tbi_file)
        bgz_file = utils.bgzip_vcf(tmp_file, True)
        assert utils.is_gz_file(bgz_file)
        files = os.listdir(tmp_dir)
        assert 's1.vcf.bgz.csi' not in files
        assert 's1.vcf.bgz.tbi' not in files


def test_index_vcf():
    vcf_file = Path(helpers.test_dir / 'raw.Sample_1.preprocessed.vcf')
    with tempfile.TemporaryDirectory() as td:
        tmp_dir: Path = Path(td)
        tmp_file = tmp_dir / 's1.vcf'
        shutil.copyfile(vcf_file, tmp_file)
        bgz_file = utils.bgzip_vcf(tmp_file, True)
        utils.index_vcf(bgz_file, True)
        files = os.listdir(tmp_dir)
        assert 's1.vcf.bgz.csi' in files
        assert 's1.vcf.bgz.tbi' not in files


def test_delete_vcf_and_index():
    with tempfile.TemporaryDirectory() as td:
        tmp_dir: Path = Path(td)
        vcf_file = tmp_dir / 'foo.vcf'
        assert not vcf_file.exists()
        utils.delete_vcf_and_index(vcf_file, True)

        csi_file = Path(str(vcf_file) + '.csi')
        helpers.touch(vcf_file)
        helpers.touch(csi_file)
        assert vcf_file.exists()
        assert csi_file.exists()
        utils.delete_vcf_and_index(vcf_file)
        assert not vcf_file.exists()
        assert not csi_file.exists()


@pytest.mark.skipif(not helpers.NETWORK_AVAILABLE, reason='No network')
def test_download_from_url_fail():
    with tempfile.TemporaryDirectory() as td:
        tmp_dir: Path = Path(td)
        with pytest.raises(InvalidURL) as context:
            preprocessor.download_from_url('https://no.such.org', tmp_dir)
        assert 'https://no.such.org' in context.value.msg


@pytest.mark.skipif(not helpers.TEST_DOWNLOAD, reason='Download disabled for tests')
def test_download_reference_fasta_and_index():
    with tempfile.TemporaryDirectory() as td:
        tmp_dir: Path = Path(td)
        ref_file = utils.download_reference_fasta_and_index(tmp_dir, verbose=1)
        files = os.listdir(tmp_dir)
        assert ref_file.name in files


@pytest.mark.skipif(not helpers.NETWORK_AVAILABLE, reason='No network')
def test_download_pharmcat_positions():
    orig_version = common.PHARMCAT_VERSION
    common.PHARMCAT_VERSION = 'nope'
    with pytest.raises(ReportableException) as context:
        with tempfile.TemporaryDirectory() as td:
            tmp_dir: Path = Path(td)
            utils.download_pharmcat_positions(tmp_dir, verbose=1)
    assert 'Cannot find pharmcat_positions file' in context.value.msg
    common.PHARMCAT_VERSION = orig_version

    with tempfile.TemporaryDirectory() as td:
        tmp_dir: Path = Path(td)
        utils.download_pharmcat_positions(tmp_dir, verbose=1)
        files = os.listdir(tmp_dir)
        assert common.PHARMCAT_POSITIONS_FILENAME in files
        assert ('%s.csi' % common.PHARMCAT_POSITIONS_FILENAME) in files


def test_prep_pharmcat_positions():
    with tempfile.TemporaryDirectory() as td:
        tmp_dir: Path = Path(td)
        os.chdir(tmp_dir)

        # no pharmcat_positions - should fail
        with pytest.raises(ReportableException) as context:
            utils.prep_pharmcat_positions(None, None, verbose=1)
        # print(context.value)
        assert 'Cannot find pharmcat_positions.vcf.bgz' in context.value.msg

        # add pharmcat positions to tmp_dir
        tmp_positions = tmp_dir / helpers.pharmcat_positions_file.name
        shutil.copyfile(helpers.pharmcat_positions_file, tmp_positions)
        tmp_uniallelic = tmp_dir / 'pharmcat_positions.uniallelic.vcf.bgz'

        reference_fasta: Path
        if helpers.TEST_DOWNLOAD:
            # use pharmcat_positions from cwd and download reference fasta
            start = timer()
            utils.prep_pharmcat_positions(None, None, verbose=1)
            full_time = timer() - start
            print("time for full preparation:", full_time)
            reference_fasta = tmp_dir / preprocessor.REFERENCE_FASTA_FILENAME
            assert reference_fasta.is_file()

            assert tmp_uniallelic.is_file()
            tmp_uniallelic.unlink()
            tmp_uniallelic_index = Path(str(tmp_uniallelic) + '.csi')
            assert tmp_uniallelic_index.is_file()
            tmp_uniallelic_index.unlink()
        else:
            # assumes that reference fasta is available next to pharmcat_positions
            reference_fasta = helpers.pharmcat_positions_file.parent / preprocessor.REFERENCE_FASTA_FILENAME
            assert reference_fasta.is_file(), 'Cannot find reference FASTA for testing!'

        # already have reference fasta, so this should be faster, but still need to generate uniallelic positions
        start = timer()
        utils.prep_pharmcat_positions(tmp_positions, reference_fasta, verbose=1)
        build_time = timer() - start
        print("time to build uniallelic positions:", build_time)
        assert build_time < 2

        tmp_uniallelic = tmp_dir / preprocessor.UNIALLELIC_VCF_FILENAME
        assert tmp_uniallelic.is_file()
        uniallelic_mtime = tmp_uniallelic.stat().st_mtime

        # already have reference fasta and uniallelic positions, so this should be very fast
        start = timer()
        # utils.prep_pharmcat_positions(None, tmp_reference, verbose=True)
        utils.prep_pharmcat_positions(tmp_positions, reference_fasta, verbose=1)
        check_time = timer() - start
        print("time to run check preparation:", check_time)
        assert check_time < build_time
        # is an order of magnitude faster
        assert check_time < (build_time / 10)

        assert tmp_uniallelic.is_file()
        assert uniallelic_mtime == tmp_uniallelic.stat().st_mtime


def test_extract_pgx_regions():
    vcf_file = helpers.test_dir / 'raw.vcf.bgz'
    vcf_file1 = helpers.test_dir / 'raw-p1.vcf.bgz'
    vcf_file2 = helpers.test_dir / 'raw-p2.vcf.bgz'
    with tempfile.TemporaryDirectory() as td:
        tmp_dir = Path(td)
        tmp_positions = tmp_dir / helpers.pharmcat_positions_file.name
        shutil.copyfile(helpers.pharmcat_positions_file, tmp_positions)
        shutil.copyfile(preprocessor.CHR_RENAME_FILE, tmp_dir / preprocessor.CHR_RENAME_MAP_FILENAME)
        tmp_vcf = tmp_dir / vcf_file.name
        shutil.copyfile(vcf_file, tmp_vcf)
        tmp_vcf1 = tmp_dir / vcf_file1.name
        shutil.copyfile(vcf_file1, tmp_vcf1)
        tmp_vcf2 = tmp_dir / vcf_file2.name
        shutil.copyfile(vcf_file2, tmp_vcf2)

        vcf_files: List[Path] = [tmp_vcf1, tmp_vcf2]
        samples: List[str] = ['Sample_1', 'Sample_2']
        combo_pgx_vcf_file = utils.extract_pgx_regions(tmp_positions, vcf_files, samples, tmp_dir, 'combo_test',
                                                       verbose=1)
        assert combo_pgx_vcf_file.is_file()
        index_file = utils.find_index_file(combo_pgx_vcf_file)
        assert index_file is not None
        combo_hash = helpers.md5hash(combo_pgx_vcf_file)
        print("combo:", combo_hash)

        pgx_vcf_file = utils.extract_pgx_regions(tmp_positions, [tmp_vcf], samples, tmp_dir, 'test', verbose=1)
        assert pgx_vcf_file.is_file()
        single_hash = helpers.md5hash(pgx_vcf_file)
        print("single:", single_hash)

        assert combo_hash == single_hash


def test_normalize_vcf():
    reference_fasta: Path = helpers.get_reference_fasta(helpers.pharmcat_positions_file)

    vcf_file = helpers.test_dir / 'test.pgx_regions.vcf.bgz'
    with tempfile.TemporaryDirectory() as td:
        tmp_dir = Path(td)
        tmp_vcf = tmp_dir / vcf_file.name
        shutil.copyfile(vcf_file, tmp_vcf)

        normalized_vcf = utils.normalize_vcf(reference_fasta, tmp_vcf, tmp_dir, 'test', verbose=1)
        assert normalized_vcf.is_file()
        # shutil.copyfile(normalized_vcf, vcf_file.parent / 'test.normalized.vcf.bgz')


def test_extract_pgx_variants():
    reference_fasta: Path = helpers.get_reference_fasta(helpers.pharmcat_positions_file)

    vcf_file = helpers.test_dir / 'test.normalized.vcf.bgz'
    with tempfile.TemporaryDirectory() as td:
        tmp_dir = Path(td)
        tmp_vcf = tmp_dir / vcf_file.name
        shutil.copyfile(vcf_file, tmp_vcf)

        multiallelic_vcf = utils.extract_pgx_variants(helpers.pharmcat_positions_file, reference_fasta, tmp_vcf,
                                                      tmp_dir, 'test', verbose=1)
        assert multiallelic_vcf.is_file()
        # shutil.copyfile(multiallelic_vcf, vcf_file.parent / multiallelic_vcf.name)


def test_export_single_sample():
    vcf_file = helpers.test_dir / 'raw.preprocessed.vcf'
    s1_file = helpers.test_dir / 'raw.Sample_1.preprocessed.vcf'
    s2_file = helpers.test_dir / 'raw.Sample_2.preprocessed.vcf'
    with tempfile.TemporaryDirectory() as td:
        tmp_dir = Path(td)
        tmp_vcf = tmp_dir / vcf_file.name
        shutil.copyfile(vcf_file, tmp_vcf)

        basename = 'test1'
        utils.export_single_sample_vcf(tmp_vcf, ['Sample_1', 'Sample_2'], tmp_dir, basename)

        helpers.compare_vcf_files(s1_file, tmp_dir, basename, 'Sample_1')
        helpers.compare_vcf_files(s2_file, tmp_dir, basename, 'Sample_2')


def test_export_single_sample_partial():
    vcf_file = helpers.test_dir / 'raw.preprocessed.vcf'
    s1_file = helpers.test_dir / 'raw.Sample_1.preprocessed.vcf'
    with tempfile.TemporaryDirectory() as td:
        tmp_dir = Path(td)
        tmp_vcf = tmp_dir / vcf_file.name
        shutil.copyfile(vcf_file, tmp_vcf)

        basename = 'test_partial'
        utils.export_single_sample_vcf(tmp_vcf, ['Sample_1'], tmp_dir, basename)

        for file in os.listdir(tmp_dir):
            print(file)

        helpers.compare_vcf_files(s1_file, tmp_dir, basename, split_sample=True)
        tmp_s2: Path = tmp_dir / ('%s.Sample_2.preprocessed.vcf' % basename)
        assert not tmp_s2.is_file(), '%s found!' % tmp_s2


def test_export_single_sample_concurrent():
    vcf_file = helpers.test_dir / 'raw.preprocessed.vcf'
    s1_file = helpers.test_dir / 'raw.Sample_1.preprocessed.vcf'
    s2_file = helpers.test_dir / 'raw.Sample_2.preprocessed.vcf'
    with tempfile.TemporaryDirectory() as td:
        tmp_dir = Path(td)
        tmp_vcf = tmp_dir / vcf_file.name
        shutil.copyfile(vcf_file, tmp_vcf)

        basename = 'test2'
        utils.export_single_sample_vcf(tmp_vcf, ['Sample_1', 'Sample_2'], tmp_dir, basename, concurrent_mode=True,
                                       max_processes=2)

        helpers.compare_vcf_files(s1_file, tmp_dir, basename, 'Sample_1')
        helpers.compare_vcf_files(s2_file, tmp_dir, basename, 'Sample_2')


def test_check_max_memory():
    assert utils.check_max_memory(None) is None
    assert utils.check_max_memory('') is None
    assert utils.check_max_memory('4M') == '4M'
    assert utils.check_max_memory('87G') == '87G'
    assert utils.check_max_memory('19m') == '19m'

    with pytest.raises(ReportableException):
        utils.check_max_memory('foo')
    with pytest.raises(ReportableException):
        utils.check_max_memory('87gb')
    with pytest.raises(ReportableException):
        utils.check_max_memory('9.7G')
