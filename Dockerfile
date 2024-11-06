# syntax=docker/dockerfile:1
#
# Base Dockerfile for PharmCAT
#
FROM python:3.12

# apt-utils needed due to https://github.com/phusion/baseimage-docker/issues/319
# liblapack-dev and libatlas-base-dev for compiling plink
RUN apt-get update && \
    apt-get install -y --no-install-recommends apt-utils apt-transport-https gpg && \
    apt-get -y upgrade && \
    apt-get -y install bzip2 build-essential wget liblapack-dev libatlas-base-dev


# install java (https://blog.adoptium.net/2021/12/eclipse-temurin-linux-installers-available/)
RUN curl https://packages.adoptium.net/artifactory/api/gpg/key/public | gpg --dearmor -o /etc/apt/keyrings/adoptium.gpg && \
    echo "deb [signed-by=/etc/apt/keyrings/adoptium.gpg] https://packages.adoptium.net/artifactory/deb $(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" \
      | tee /etc/apt/sources.list.d/adoptium.list && \
    apt-get update && \
    apt-get -y install --no-install-recommends temurin-17-jdk

# install google cloud utils
RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /etc/apt/keyrings/cloud.google.gpg && \
    echo "deb [signed-by=/etc/apt/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main"  \
      | tee /etc/apt/sources.list.d/google-cloud-sdk.list && \
    apt-get update && \
    apt-get -y install google-cloud-cli

# install DNAnexus tools
RUN pip3 install dxpy


WORKDIR /usr/local/bin/

ENV BCFTOOLS_VERSION=1.21
ENV HTSLIB_VERSION=1.21
ENV SAMTOOLS_VERSION=1.21
# download the suite of tools
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
# extract files
RUN tar -xjf /usr/local/bin/htslib-${HTSLIB_VERSION}.tar.bz2 -C /usr/local/bin/ && \
    tar -xjf /usr/local/bin/bcftools-${BCFTOOLS_VERSION}.tar.bz2 -C /usr/local/bin/ && \
    tar -xjf /usr/local/bin/samtools-${SAMTOOLS_VERSION}.tar.bz2 -C /usr/local/bin/
# compile
RUN cd /usr/local/bin/htslib-${HTSLIB_VERSION}/ && ./configure && \
    cd /usr/local/bin/htslib-${HTSLIB_VERSION}/ && make && make install && \
    cd /usr/local/bin/bcftools-${BCFTOOLS_VERSION}/ && ./configure && \
    cd /usr/local/bin/bcftools-${BCFTOOLS_VERSION}/ && make && make install && \
    cd /usr/local/bin/samtools-${SAMTOOLS_VERSION}/ && ./configure && \
    cd /usr/local/bin/samtools-${SAMTOOLS_VERSION}/ && make && make install
# cleanup
RUN rm -rf /usr/local/bin/bcftools-${BCFTOOLS_VERSION} /usr/local/bin/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
           /usr/local/bin/htslib-${HTSLIB_VERSION}     /usr/local/bin/htslib-${HTSLIB_VERSION}.tar.bz2 \
           /usr/local/bin/samtools-${SAMTOOLS_VERSION} /usr/local/bin/samtools-${SAMTOOLS_VERSION}.tar.bz2

# install plink2
RUN git clone https://github.com/chrchang/plink-ng.git && \
    cd plink-ng/2.0 && ./build.sh && \
    mv bin/* /usr/local/bin && \
    rm -rf /usr/local/bin/plink-ng



RUN mkdir /pharmcat
WORKDIR /pharmcat
ENV PATH="$PATH:/pharmcat"

# download fasta files
RUN wget https://zenodo.org/record/7288118/files/GRCh38_reference_fasta.tar && \
    tar -xf GRCh38_reference_fasta.tar --no-same-owner && \
    rm -f GRCh38_reference_fasta.tar

# setup python env
COPY preprocessor/requirements.txt ./
RUN pip3 install -r requirements.txt && \
    rm requirements.txt

# setup user
COPY src/main/config/bashrc /root/.bashrc

# add pharmcat scripts
COPY preprocessor/pharmcat_vcf_preprocessor.py \
     preprocessor/pharmcat_pipeline \
     bin/pharmcat \
     build/pharmcat.jar \
     pharmcat_positions.vcf* \
     pharmcat_regions.bed \
     ./
RUN mkdir preprocessor
COPY preprocessor/preprocessor/*.py \
     preprocessor/preprocessor/*.tsv \
     preprocessor/
RUN mkdir data && \
    chmod 755 *.py pharmcat pharmcat_pipeline pharmcat_vcf_preprocessor.py preprocessor data && \
    python -c "import preprocessor; preprocessor.prep_pharmcat_positions()"

CMD ["/bin/bash"]
