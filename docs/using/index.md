---
title: Using PharmCAT
permalink: using/
nav_order: 3
has_children: true
has_toc: false
---
# Using PharmCAT
{: .no_toc }

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

## Quick Start

The quickest and easiest way to get started with PharmCAT is through
[Docker](https://docs.docker.com/get-started/docker-overview/).

Once you have Docker [installed](https://docs.docker.com/get-docker/), you can get PharmCAT from
[Docker Hub](https://hub.docker.com/r/pgkb/pharmcat):

```console
# docker pull pgkb/pharmcat
```

##### Run PharmCAT Pipeline
{: .no_toc }

You will need a VCF file to run PharmCAT.
Download our [sample VCF](/examples/pharmcat.example.vcf) and save it as `sample.vcf`.
Assuming you have `sample.vcf` in `/path/to/data`, you can run the PharmCAT pipeline with:

```console
# docker run --rm -v /path/to/data:/pharmcat/data pgkb/pharmcat pharmcat_pipeline data/sample.vcf
```

* `docker run`: The base Docker command
* `--rm`: Cleans up the container automatically when you're done with it
* `-v`: Bind mounts `/path/to/data` on your machine to `/pharmcat/data` in the Docker image.
  This will make the data available under the `data` subdirectory.
* `pgkb/pharmcat`: The name of the PharmCAT image
* `pharmcat_pipeline`: The PharmCAT pipeline script
* `data/sample.vcf`: The VCF file to process


And you should see the following output:  
```
PharmCAT version: {{site.pharmcat_version}}

Processing data/sample.vcf ...

Running PharmCAT...
Checking files...
* Found 1 VCF file
Saving named allele matcher JSON results to /pharmcat/data/sample.match.json
Saving phenotyper JSON results to /pharmcat/data/sample.phenotype.json
Saving reporter HTML results to /pharmcat/data/sample.report.html

Done.
```

## Inputs

Before continuing with PharmCAT, please take the time to understand PharmCAT's input requirements.
PharmCAT accepts two types of files: VCF files and "outside call" files.

VCF files contain genetic sequence variation and must meet __[PharmCAT's VCF Requirements](/using/VCF-Requirements)__.
This is used by PharmCAT to identify pharmacogenomic (PGx) genotypes and infer haplotypes, typically called star
alleles.

There are some genes, such as CYP2D6, HLA-A and HLA-B that PharmCAT cannot call. To get PharmCAT to generate guideline
recommendations for these genes, you will need to provide the diplotype calls directly to PharmCAT using an
__[Outside Call file](/using/Outside-Call-Format)__ (so named because the call was made outside PharmCAT).
You might also use this if you want to get guideline recommendations without using PharmCAT to call diplotypes from
VCF data.  Consult the [documentation](/using/Outside-Call-Format) for details.

* Details on [Calling CYP2D6](/using/Calling-CYP2D6)
* Details on [Calling HLA](/using/Calling-HLA)


## Docker

The quickest and easiest way to get started with PharmCAT is through
[Docker](https://docs.docker.com/get-started/docker-overview/).

See [PharmCAT in Docker](/using/PharmCAT-in-Docker) for a more complete walk-through on how to do so.


## PharmCAT Pipeline

The easiest way to run PharmCAT is to use the `pharmcat_pipeline` script.

[Running the PharmCAT Pipeline](/using/Running-PharmCAT-Pipeline) has all the details on using the script and plenty
of examples.

Note: the `pharmcat_pipeline` script caters for simplicity and ease of use.
If you really want to take advantage of PharmCAT, take a look at the [Advanced Usage](#advanced-usage) section below.


## Advanced Usage

The PharmCAT pipeline is composed of two components: the VCF Preprocessor and the core PharmCAT tool.
You can run either of these components separately.
In fact, the core PharmCAT tool is itself composed of multiple modules that can be run independently as well.

For details on using these components, see:

* [Running the VCF Preprocessor](/using/VCF-Preprocessor)
* [Running PharmCAT](/using/Running-PharmCAT)

If you do not wish to use Docker, or you need to set it up in your own environment, then take a look at
[Advanced PharmCAT Setup](/using/Setup-PharmCAT).


### Video Tutorial
{: .no_toc }

Prefer video tutorials?
We have a hands-on video tutorial on [how to run PharmCAT](https://youtu.be/d1IZPLOrPOE?si=LREY8RI-wz-5PoqN) on YouTube.
It will walk you through setting up and running PharmCAT.

Note that this tutorial does not cover the `pharmcat_pipeline` script, but using the VCF preprocessor and core PharmCAT
tool directly.

<iframe width="560" height="315" src="https://www.youtube-nocookie.com/embed/d1IZPLOrPOE?si=iq-8xFCpGik9C3w3&amp;controls=0" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>

### Going Further
{: .no_toc }

* [Multi-Sample Analysis](/using/Multi-Sample-Analysis)
* [Research Mode](/using/Research-Mode)

 
![haplocat](/images/haplocat.png)
