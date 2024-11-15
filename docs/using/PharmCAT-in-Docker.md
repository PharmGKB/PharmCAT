---
parent: Using PharmCAT
title: PharmCAT in Docker
permalink: using/PharmCAT-in-Docker/
nav_order: 1
render_with_liquid: false
---
# PharmCAT in Docker

PharmCAT is available in a [Docker container](https://hub.docker.com/r/pgkb/pharmcat).
If you are not familiar with Docker, this [overview](https://docs.docker.com/get-started/overview/) is a good starting
point.

This page will cover the basics of getting started with Docker and PharmCAT.


## Setup

If you don't already have Docker installed,
[follow these instructions to install Docker](https://docs.docker.com/get-docker/).

Then you can get PharmCAT from [Docker Hub](https://hub.docker.com/r/pgkb/pharmcat):

```console
# docker pull pgkb/pharmcat
```

## Usage

You will need to make your data accessible to the Docker container. There are
[several options](https://docs.docker.com/storage/) to choose from, and you will have to decide what works best for you.
For example, a volume mount is the best for persisting data but will take some configuration.

This tutorial will use bind mounts because they are the easiest to use and require no prior configuration.

There are two ways to use the Docker image: on a per-command basis or interactively.
Choose the former if you only want to run PharmCAT as a one-off event.
If you intend to run multiple PharmCAT commands or are exploring the tool, then interactive mode would be a better 
option (and involve less typing).  


### Per-command usage

Use this to run PharmCAT with single commands.  

General usage:

```console
# docker run --rm -v /path/to/data:/pharmcat/data pgkb/pharmcat <xxx>
```

* __docker run__: The base Docker command
* __--rm__: Cleans up the container automatically when you're done with it
* __-v__: Bind mounts `/path/to/data` on your machine to `/pharmcat/data` in the Docker image.
This will make the data available under the `data` subdirectory.
* __pgkb/pharmcat__: The name of the PharmCAT image
* __&lt;xxx&gt;__: Command to run

##### Example 

Assuming you have a VCF file called `sample.vcf` in `/path/to/data`, you can run the
[PharmCAT Pipeline](/using/Running-PharmCAT-Pipeline) with:

```console
# docker run --rm -v /path/to/data:/pharmcat/data pgkb/pharmcat pharmcat_pipeline /pharmcat/data/sample.vcf
PharmCAT version: {{site.pharmcat_version}}

Processing /pharmcat/data/sample.vcf ...

Running PharmCAT...
Checking files...
* Found 1 VCF file
Saving named allele matcher JSON results to /pharmcat/data/sample.match.json
Saving phenotyper JSON results to /pharmcat/data/sample.phenotype.json
Saving reporter HTML results to /pharmcat/data/sample.report.html

Done.
#
```


### Interactive mode

Use this if you want to run PharmCAT multiple times or just explore the available options.

To start interactive mode:

```console
# docker run --rm -v /path/to/data:/pharmcat/data -it pgkb/pharmcat
```
* __docker run__: The base Docker command
* __--rm__: Cleans up the container automatically when you're done with it
* __-v__: Bind mounts `/path/to/data` on your machine to `/pharmcat/data` in the Docker image.
This will make the data available under the `data` subdirectory.
* __-it__: Starts interactive mode
* __pgkb/pharmcat__: The name of the PharmCAT image

##### Example

Once you are in interactive mode, you can proceed to use PharmCAT.
Assuming you have a VCF file called `sample.vcf` in `/path/to/data`, you can run the
[PharmCAT Pipeline](/using/Running-PharmCAT-Pipeline) with:

```console
# docker run --rm -v /path/to/data:/pharmcat/data -it pgkb/pharmcat
/pharmcat > pharmcat_pipeline /pharmcat/data/sample.vcf
PharmCAT version: {{site.pharmcat_version}}

Processing /pharmcat/data/sample.vcf ...

Running PharmCAT...
Checking files...
* Found 1 VCF file
Saving named allele matcher JSON results to /pharmcat/data/sample.match.json
Saving phenotyper JSON results to /pharmcat/data/sample.phenotype.json
Saving reporter HTML results to /pharmcat/data/sample.report.html

Done.
/pharmcat >
```


## Next Steps

Once you have Docker set up, learn more about how to use PharmCAT:

* [PharmCAT Pipeline](/using/Running-PharmCAT-Pipeline)
