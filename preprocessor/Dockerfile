FROM continuumio/anaconda


MAINTAINER "Lester Carter" lester@pharmgkb.org

# update
RUN apt-get -y update
RUN apt-get -y install zlib1g-dev vim python-software-properties apt-file software-properties-common build-essential curl git python-setuptools ruby

# Install bioconda packages
RUN conda config --add channels bioconda && conda install -y -c daler tabix=0.2.6  && conda install -y bedops=2.4.20 && conda install -y -c jakirkham npm=1.4.12 && conda install -y -c bcbio crossmap=0.1.8 && conda install -y vcftools=0.1.14 && conda install -y bcftools=1.3.1 && conda install -y -c bioconda java-jdk=8.0.92 && conda install -y -c conda-cluster apache-maven=3.2.1

RUN echo "Get PharmCAT"
# Get PharmCAT
RUN git clone https://github.com/PharmGKB/PharmCAT.git
WORKDIR /PharmCAT
RUN ./gradlew shadowJar
WORKDIR /

# get pgx positions
RUN java -cp /PharmCAT/build/libs/pharmcat-*-all.jar org.pharmgkb.pharmcat.definition.ExtractPositions  -g hg38 -d /PharmCAT/src/main/resources/org/pharmgkb/pharmcat/definition/alleles/  -o /pgx.vcf

RUN bgzip -c pgx.vcf  > pgx.vcf.gz && tabix -p vcf pgx.vcf.gz

# Add prepare.sh
ADD prepare.sh prepare.sh

# Add parse.py
ADD parse.py parse.py


ENTRYPOINT ["/prepare.sh"]
