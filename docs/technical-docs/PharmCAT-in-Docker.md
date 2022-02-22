---
title: Docker
permalink: technical-docs/docker
parent: Using PharmCAT
nav_order: 3
---
# Docker

PharmCAT is available in a Docker container.

## Setup

If you are not familiar with Docker, this [overview](https://docs.docker.com/get-started/overview/) is a good starting point.

You must have Docker [installed](https://docs.docker.com/get-docker/) to use PharmCAT via Docker.

Then you can get PharmCAT from [Docker Hub](https://hub.docker.com/repository/docker/pgkb/pharmcat):

```commandline
# docker pull pgkb/pharmcat
```

## Usage

You will need to make your data accessible to the Docker container.  There are [many options](https://docs.docker.com/storage/), of which a volume mount is probably the best.

However, this tutorial will use bind mounts.

General usage:

```commandline
# docker run --rm -v /path/to/data:/pharmcat/data pgkb/pharmcat xxx
```

* `--rm` - cleans up the container automatically when you're done with it
* `-v` - bind mounts `/path/to/data` on your machine to `/pharmcat/data` in the Docker image.  This will make the data available under the `data` subdirectory.
* `xxx` - command to run

If you run `ls`, it will list the contents of the `/pharmcat` directory: 

```commandline
# docker run --rm -v /path/to/data:/pharmcat/data pgkb/pharmcat ls
GRCh38_reference_fasta.tar
PharmCAT_VCF_Preprocess.py
PharmCAT_VCF_Preprocess_py3_requirements.txt
data
pharmcat
pharmcat.jar
pharmcat_positions.vcf
pharmcat_positions.vcf.bgz
pharmcat_positions.vcf.bgz.tbi
reference.fasta.bgz
reference.fasta.bgz.fai
reference.fasta.bgz.gzi
vcf_preprocess_exceptions.py
vcf_preprocess_utilities.py
```

### Running the VCF preprocessor

Your VCF files needs to comply with [PharmCAT's requirements](VCF-Requirements).  PharmCAT has a [VCF preprocessor](Preprocessing-VCF-Files-for-PharmCAT) that will handle much of this for you.  Please consult the [preprocessor documention](Preprocessing-VCF-Files-for-PharmCAT) for details.

```commandline
# docker run --rm -v /path/to/data:/pharmcat/data pgkb/pharmcat ./PharmCAT_VCF_Preprocess.py
```

If you have a file `/path/to/data/sample.vcf`, you would use:

```commandline
# docker run --rm -v /path/to/data:/pharmcat/data pgkb/pharmcat ./PharmCAT_VCF_Preprocess.py --input_vcf data/sample.vcf
```

Note: the GRCh38 reference is included in the Docker image, so you do not need to provide it unless you have special reference requirements.


### Running PharmCAT

```commandline
# docker run --rm -v /path/to/data:/pharmcat/data pgkb/pharmcat ./pharmcat
```

After running the file `/path/to/data/sample.vcf` through the preprocessor, assuming the sample ID was "SAMPLE1", you would have gotten a file called `pharmcat_ready_vcf.SAMPLE1.vcf`.  You can then run this through PharmCAT with:

```commandline
# docker run --rm -v /path/to/data:/pharmcat/data pgkb/pharmcat ./pharmcat -vcf data/pharmcat_ready_vcf.SAMPLE1.vcf
```


> The Docker image includes the `pharmcat` script, which is just a wrapper around the call to Java.
> For details on using PharmCAT, please see the [Running PharmCAT](docs/technical-docs/Running-PharmCAT.md#running)
