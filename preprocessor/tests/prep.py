import os
import re
import shutil
import sys
import tempfile
from pathlib import Path
from typing import List


# get the absolute path to the parent directory containing this file
module_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
# add the parent directory to the system path so we can import pcat
sys.path.insert(0, module_dir)

from pcat import bgzip_vcf, common, ReportableException, run
from tests.helpers import test_dir


if __name__ == "__main__":
    raw_vcf = test_dir / 'raw.vcf'
    raw_bgz = test_dir / 'raw.vcf.bgz'
    raw1_vcf = test_dir / 'raw-p1.vcf'
    raw1_bgz = test_dir / 'raw-p1.vcf.bgz'
    raw2_vcf = test_dir / 'raw-p2.vcf'
    raw2_bgz = test_dir / 'raw-p2.vcf.bgz'
    raw_bcf = test_dir / 'raw.bcf'
    raw_bcf_bgzf = test_dir / 'raw.bcf.bgzf'
    raw_chr_vcf = test_dir / 'raw_chr.vcf'
    raw_chr_bcf_bgzf = test_dir / 'raw_chr.bcf.bgzf'
    raw_bad_chr_vcf = test_dir / 'raw_bad_chr.vcf'
    raw_bad_chr_bcf_bgzf = test_dir / 'raw_bad_chr.bcf.bgzf'
    with tempfile.TemporaryDirectory() as td:
        tmp_dir: Path = Path(td)

        # bgzip raw file
        tmp_vcf = tmp_dir / 'raw.vcf'
        shutil.copyfile(raw_vcf, tmp_vcf)
        tmp_bgz = bgzip_vcf(tmp_vcf)
        shutil.copyfile(tmp_bgz, raw_bgz)

        # generate BCF file
        bcftools_command = [common.BCFTOOLS_PATH, 'view', '--no-version', '-Ob', '-o', str(raw_bcf_bgzf),
                            str(raw_vcf)]
        run(bcftools_command)
        bcftools_command = [common.BCFTOOLS_PATH, 'view', '--no-version', '-Ou', '-o', str(raw_bcf),
                            str(raw_vcf)]
        run(bcftools_command)
        bcftools_command = [common.BCFTOOLS_PATH, 'view', '--no-version', '-Ou', '-o', str(raw_chr_bcf_bgzf),
                            str(raw_chr_vcf)]
        run(bcftools_command)
        bcftools_command = [common.BCFTOOLS_PATH, 'view', '--no-version', '-Ou', '-o', str(raw_bad_chr_bcf_bgzf),
                            str(raw_bad_chr_vcf)]
        run(bcftools_command)

        # split raw file in half and bgzip them
        comments: List[str] = []
        data1: List[str] = []
        data2: List[str] = []
        pattern = re.compile(r'(?:chr)?([XYM]|\d+)\s.*')
        with open(raw_vcf, 'r') as f:
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

        with open(raw1_vcf, 'w') as f:
            f.write('\n'.join(comments))
            f.write('\n')
            f.write('\n'.join(data1))
            f.write('\n')
        bgzip_vcf(raw1_vcf)

        with open(raw2_vcf, 'w') as f:
            f.write('\n'.join(comments))
            f.write('\n')
            f.write('\n'.join(data2))
            f.write('\n')
        bgzip_vcf(raw2_vcf)
