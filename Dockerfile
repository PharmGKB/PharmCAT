# syntax=docker/dockerfile:1
#
# Base Dockerfile for PharmCAT
#
FROM python:3.9

# apt-utils line due to https://github.com/phusion/baseimage-docker/issues/319
RUN apt-get update && \
    apt-get install -y --no-install-recommends apt-utils apt-transport-https gnupg && \
    apt-get -y upgrade && \
    apt-get -y install bzip2 build-essential wget

# install java (https://blog.adoptium.net/2021/12/eclipse-temurin-linux-installers-available/)
RUN wget -O - https://packages.adoptium.net/artifactory/api/gpg/key/public | apt-key add -
RUN echo "deb https://packages.adoptium.net/artifactory/deb $(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" \
    | tee /etc/apt/sources.list.d/adoptium.list
RUN apt-get update && \
    apt-get -y install --no-install-recommends temurin-17-jdk


RUN mkdir /pharmcat
WORKDIR /pharmcat
# download fasta files
RUN wget https://zenodo.org/record/7251599/files/GRCh38_reference_fasta.tar && \
    tar -xf GRCh38_reference_fasta.tar && \
    rm -f GRCh38_reference_fasta.tar


ENV BCFTOOLS_VERSION 1.16
ENV HTSLIB_VERSION 1.16
ENV SAMTOOLS_VERSION 1.16

# download the suite of tools
WORKDIR /usr/local/bin/
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2

# extract files for the suite of tools
RUN tar -xjf /usr/local/bin/htslib-${HTSLIB_VERSION}.tar.bz2 -C /usr/local/bin/
RUN tar -xjf /usr/local/bin/bcftools-${BCFTOOLS_VERSION}.tar.bz2 -C /usr/local/bin/
RUN tar -xjf /usr/local/bin/samtools-${SAMTOOLS_VERSION}.tar.bz2 -C /usr/local/bin/

# compile tools
RUN cd /usr/local/bin/htslib-${HTSLIB_VERSION}/ && ./configure
RUN cd /usr/local/bin/htslib-${HTSLIB_VERSION}/ && make && make install
RUN cd /usr/local/bin/bcftools-${BCFTOOLS_VERSION}/ && ./configure
RUN cd /usr/local/bin/bcftools-${BCFTOOLS_VERSION}/ && make && make install
RUN cd /usr/local/bin/samtools-${SAMTOOLS_VERSION}/ && ./configure
RUN cd /usr/local/bin/samtools-${SAMTOOLS_VERSION}/ && make && make install

# cleanup
RUN rm  -f /usr/local/bin/bcftools-${BCFTOOLS_VERSION}.tar.bz2
RUN rm -rf /usr/local/bin/bcftools-${BCFTOOLS_VERSION}
RUN rm  -f /usr/local/bin/htslib-${HTSLIB_VERSION}.tar.bz2
RUN rm -rf /usr/local/bin/htslib-${HTSLIB_VERSION}
RUN rm  -f /usr/local/bin/samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN rm -rf /usr/local/bin/samtools-${SAMTOOLS_VERSION}


WORKDIR /pharmcat
# setup python env
COPY src/scripts/preprocessor/PharmCAT_VCF_Preprocess_py3_requirements.txt ./
RUN pip3 install -r PharmCAT_VCF_Preprocess_py3_requirements.txt

# add pharmcat scripts
COPY src/scripts/preprocessor/*.py ./
COPY src/scripts/pharmcat ./
RUN chmod 755 *.py
RUN chmod 755 pharmcat
COPY pharmcat_positions.vcf* ./
COPY build/pharmcat.jar ./
