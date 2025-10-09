#!/bin/bash
#
# Sets up samtools.
# Requires the following libraries on Ubuntu: libbz2-dev, libncurses-dev, liblzma-dev.
#
# Can be called with a single parameter, defining where to install (will be supplied as --prefix to configure).
# Defaults to current directory.
#

set -e
set -u
set -o pipefail

BCFTOOLS_VERSION=1.22
HTSLIB_VERSION=1.22
SAMTOOLS_VERSION=1.22

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
if [ $# -eq 0 ]; then
  install_dir=$cwd
else
  install_dir=$1
fi
cd htslib-${HTSLIB_VERSION}/     && ./configure --prefix="$install_dir" && make && make install
cd "$cwd"
cd bcftools-${BCFTOOLS_VERSION}/ && ./configure --prefix="$install_dir" && make && make install
cd "$cwd"
cd samtools-${SAMTOOLS_VERSION}/ && ./configure --prefix="$install_dir" && make && make install

# cleanup
cd "$cwd"
rm -r htslib-${HTSLIB_VERSION}* bcftools-${BCFTOOLS_VERSION}* samtools-${SAMTOOLS_VERSION}*
