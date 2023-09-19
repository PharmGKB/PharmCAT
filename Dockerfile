# syntax=docker/dockerfile:1
#
# Base Dockerfile for PharmCAT
#
FROM python:3.10

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
RUN wget https://zenodo.org/record/7288118/files/GRCh38_reference_fasta.tar && \
    tar -xf GRCh38_reference_fasta.tar --no-same-owner && \
    rm -f GRCh38_reference_fasta.tar


ENV BCFTOOLS_VERSION 1.18
ENV HTSLIB_VERSION 1.18
ENV SAMTOOLS_VERSION 1.18

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

# setup python env
COPY preprocessor/requirements.txt ./
RUN pip3 install -r requirements.txt
RUN rm requirements.txt

# setup user
COPY src/main/config/bashrc /root/.bashrc

WORKDIR /pharmcat
# add pharmcat scripts
COPY preprocessor/pharmcat_vcf_preprocessor.py \
     preprocessor/pharmcat_pipeline \
     bin/pharmcat \
     build/pharmcat.jar \
     pharmcat_positions.vcf* \
     ./
RUN mkdir preprocessor
COPY preprocessor/preprocessor/*.py \
     preprocessor/preprocessor/*.tsv \
     preprocessor/
RUN mkdir data
RUN chmod 755 *.py pharmcat pharmcat_pipeline pharmcat_vcf_preprocessor.py preprocessor data
RUN python -c "import preprocessor; preprocessor.prep_pharmcat_positions()"
