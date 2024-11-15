---
parent: Using PharmCAT
title: Advanced PharmCAT Setup
permalink: using/Setup-PharmCAT/
nav_order: 5
render_with_liquid: false
mermaid: true
---
# Advanced PharmCAT Setup

This document teaches you how to set up PharmCAT.


## Install Dependencies

Before installing PharmCAT, you will need to install the following dependencies:


1. [Java 17](https://adoptium.net/index.html?variant=openjdk17&jvmVariant=hotspot) or newer.
    _We currently recommend Java 17._
2. [Python 3.10.14](https://www.python.org/downloads/) or newer.
3. The following bioinformatic tools:
    * [bcftools >= v1.18](http://www.htslib.org/download/)
    * [htslib >= v1.18](http://www.htslib.org/download/) (for bgzip)
   You will need to have the binaries from these packages available in your PATH.

{: .info}
Mac users: bcftools and htslib are available via homebrew.


##### Environment variables

If you wish, you can customize which versions of PharmCAT's dependencies to use via environment variables:

* `JAVA_HOME` - the directory where the version of Java you want to use is installed
* `BCFTOOLS_PATH` - the full path to the `bcftools` program
* `BGZIP_PATH` - the full path to the `bgzip` program


## Install PharmCAT

If you are on a Mac or Unix system, you can run:

```console
curl -fsSL https://get.pharmcat.org | bash
```

This will check to make sure you have the required dependencies listed above and install the PharmCAT in the current
directory.

Alternatively, you can download the PharmCAT Pipeline tar file (`pharmcat-pipeline-xxx.tar.gz`) from our
[release page](https://github.com/PharmGKB/PharmCAT/releases/) and un-tar it.


Whichever way you choose, you will need to manually install the Python dependencies that the Preprocessor requires.
These dependencies are listed in a file called `requirements.txt`.
You can either use your preferred Python virtual environment to install these dependencies or run:
```console
pip3 install -r requirements.txt
```


## Done!

And that's it.  You should now have all the tools you need to run PharmCAT.

Next steps:

* Learn about [PharmCAT's VCF requirements](/using/VCF-Requirements)
* Learn how to [use the PharmCAT Pipeline](/using/Running-PharmCAT-Pipeline)
* Learn how to [use the PharmCAT VCF Preprocessor](/using/VCF-Preprocessor)
* Learn how to [use the core PharmCAT tool](/using/Running-PharmCAT)
