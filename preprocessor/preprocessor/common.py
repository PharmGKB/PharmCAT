from pathlib import Path


# version is dynamically updated - DO NOT MODIFY MANUALLY
PHARMCAT_VERSION = '2.2.1'

# expected tool versions
MIN_BCFTOOLS_VERSION = '1.16'
MIN_BGZIP_VERSION = '1.16'
MIN_JAVA_VERSION = '17'

# default filenames
PHARMCAT_POSITIONS_FILENAME = 'pharmcat_positions.vcf.bgz'
REFERENCE_FASTA_FILENAME = 'reference.fna.bgz'
UNIALLELIC_VCF_SUFFIX = '.uniallelic.vcf.bgz'
UNIALLELIC_VCF_FILENAME = 'pharmcat_positions' + UNIALLELIC_VCF_SUFFIX
CHR_RENAME_MAP_FILENAME = 'chr_rename_map.tsv'

# paths
SCRIPT_DIR: Path = Path(globals().get("__file__", "./_")).absolute().parent
CHR_RENAME_FILE: Path = SCRIPT_DIR / CHR_RENAME_MAP_FILENAME
BCFTOOLS_PATH = 'bcftools'
BGZIP_PATH = 'bgzip'
JAVA_PATH = 'java'
