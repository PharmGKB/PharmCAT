# syntax=docker/dockerfile:1
#
# Base Dockerfile for PharmCAT
#
FROM python:3.12

# apt-utils needed due to https://github.com/phusion/baseimage-docker/issues/319
RUN apt-get update && \
    apt-get install -y --no-install-recommends apt-utils apt-transport-https gpg vim && \
    apt-get -y upgrade && \
    apt-get -y install bzip2 build-essential wget


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

ENV BCFTOOLS_VERSION=1.22
ENV HTSLIB_VERSION=1.22
ENV SAMTOOLS_VERSION=1.22
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
#RUN apt-get -y install liblapack-dev libatlas-base-dev && \
#    git clone https://github.com/chrchang/plink-ng.git && \
#    cd plink-ng/2.0 && ./build.sh && \
#    mv bin/* /usr/local/bin && \
#    rm -rf /usr/local/bin/plink-ng

# install R
# hack to get debian codename from https://unix.stackexchange.com/a/253476/90413
#RUN gpg --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' && \
#    gpg --armor --export '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' | gpg --dearmor | tee /usr/share/keyrings/cran.gpg > /dev/null && \
#    echo "deb [signed-by=/usr/share/keyrings/cran.gpg] https://cloud.r-project.org/bin/linux/debian $(dpkg --status tzdata|grep Provides|cut -f2 -d'-')-cran40/" > /etc/apt/sources.list.d/cran.list


RUN mkdir /pharmcat
WORKDIR /pharmcat
ENV PATH="$PATH:/pharmcat"
ENV PCAT_PLATFORM=docker

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
COPY preprocessor/pharmcat_vcf_preprocessor \
     preprocessor/pharmcat_pipeline \
     bin/pharmcat \
     bin/calc_allele_freqs \
     build/pharmcat.jar \
     pharmcat_positions.vcf* \
     pharmcat_regions.bed \
     ./
RUN mkdir pcat
COPY preprocessor/pcat/*.py \
     preprocessor/pcat/*.tsv \
     pcat/
RUN mkdir data && \
    chmod 755 pharmcat pharmcat_pipeline pharmcat_vcf_preprocessor pcat data && \
    python -c "import pcat; pcat.prep_pharmcat_positions()"

CMD ["/bin/bash"]
