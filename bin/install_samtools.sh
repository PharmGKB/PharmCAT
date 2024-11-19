#!/bin/bash
#
#  Sets up PharmCAT in the current directory after checking for basic requirements.
#
#

set -e
set -u
set -o pipefail

BCFTOOLS_VERSION=1.21
HTSLIB_VERSION=1.21
SAMTOOLS_VERSION=1.21

# download the suite of tools
wget -nv https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2
wget -nv https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2
wget -nv https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2

# extract files
tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2
tar -xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2
tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2

# compile
cwd=$(pwd)
cd htslib-${HTSLIB_VERSION}/     && ./configure --prefix="$1" && make && make install
cd "$cwd"
cd bcftools-${BCFTOOLS_VERSION}/ && ./configure --prefix="$1" && make && make install
cd "$cwd"
cd samtools-${SAMTOOLS_VERSION}/ && ./configure --prefix="$1" && make && make install

# cleanup
cd "$cwd"
rm -r htslib-${HTSLIB_VERSION}* bcftools-${BCFTOOLS_VERSION}* samtools-${SAMTOOLS_VERSION}*
