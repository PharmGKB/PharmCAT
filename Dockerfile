# syntax=docker/dockerfile:1
#
# Base Dockerfile for PharmCAT
#
FROM python:3

# apt-utils line due to https://github.com/phusion/baseimage-docker/issues/319
RUN apt-get update && \
    apt-get install -y --no-install-recommends apt-utils && \
    apt-get -y upgrade && \
    apt-get -y install bzip2 build-essential wget


ENV BCFTOOLS_VERSION 1.13
ENV HTSLIB_VERSION 1.13

# download the suite of tools
WORKDIR /usr/local/bin/
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2

# extract files for the suite of tools
RUN tar -xjf /usr/local/bin/bcftools-${BCFTOOLS_VERSION}.tar.bz2 -C /usr/local/bin/
RUN tar -xjf /usr/local/bin/htslib-${HTSLIB_VERSION}.tar.bz2 -C /usr/local/bin/

# compile tools
RUN cd /usr/local/bin/htslib-${HTSLIB_VERSION}/ && ./configure
RUN cd /usr/local/bin/htslib-${HTSLIB_VERSION}/ && make && make install
RUN cd /usr/local/bin/bcftools-${BCFTOOLS_VERSION}/ && make && make install

# cleanup
RUN rm -f /usr/local/bin/bcftools-${BCFTOOLS_VERSION}.tar.bz2
RUN rm -rf /usr/local/bin/bcftools-${BCFTOOLS_VERSION}
RUN rm -f /usr/local/bin/htslib-${HTSLIB_VERSION}.tar.bz2
RUN rm -rf /usr/local/bin/htslib-${HTSLIB_VERSION}

# add pharmcat scritps
RUN mkdir /pharmcat
WORKDIR /pharmcat
COPY src/scripts/* ./
RUN pip3 install -r PharmCAT_VCF_Preprocess_py3_requirements.txt
RUN python -c 'import vcf_preprocess_utilities as utils; utils.download_grch38_ref_fasta_and_index("/pharmcat", "/pharmcat/reference.fasta")'
RUN rm *.fna.gz
