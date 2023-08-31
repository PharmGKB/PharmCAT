---
parent: Using PharmCAT
title: Calling CYP2D6
permalink: using/Calling-CYP2D6/
nav_order: 8
---
# Calling CYP2D6

PharmCAT supports CYP2D6 data with the following options, which are discussed in more detail below:
* use a different tool to determine CYP2D6 diplotypes and use PharmCAT to provide the translation into phenotypes and recommendations
* run PharmCAT in research mode to call CYP2D6 from VCF files

## CYP2D6 complications
### Problems with calling CYP2D6 from VCF

While PharmCAT supports CYP2D6, we **do NOT recommend calling CYP2D6 from VCF** due to the large influence of structural variation (SV) and copy number variation (CNV) on inferring CYP2D6 phenotype, which is beyond the scope of what can be called from SNPs or INDELs in a VCF file.

CYP2D6 phenotype prediction is sensitive to SV and CNV. Two or more copies of CYP2D6 on one chromosome have been reported for normal function (e.g. `*1`, `*2`), decreased function (e.g. `*10`, `*17`) and no function (e.g. `*4`, `*36`) alleles. The CYP2D6 ultrarapid metabolizer phenotype (UM, activity score >2.25) is predicted based on duplications of normal function alleles and the absence of alleles with an activity value of 0 or 0.25 (e.g. `CYP2D6*2x2/*1`, `CYP2D6*2x2/*17`) or multiplications of normal function alleles (e.g. `CYP2D6*1x3/*4`). As such, CYP2D6 UMs cannot be called using only SNPs and INDELs in a VCF file. Additionally, missing the copy number of a normal or decreased function allele can misrepresent a normal metabolizer (NM) (e.g. `*2x2/*4`, `*10/*17x2`) as an intermediate metabolizer (IM) (e.g. `*2/*4`, `*10/*17`). Moreover, SVs, such as CYP2D6 and CYP2D7 hybrid alleles (e.g. `*13`, `*68`) and the whole gene deletion (`*5`), are no function alleles. At the moment, it is still difficult to have accurate representation of these alleles defined by complex SV in a VCF. Omission or misrepresentation of these alleles, depending on other alleles, will lead to potential carrier being mistakenly reported as IM or NM.

In the specific case where a sample has the whole gene deletion (`*5`) on one CYP2D6 allele and presents variants on the other CYP2D6 allele, these hemizygous variants will be falsely presented as homozygous in the VCF, e.g. `*5/*29` will be detected as `*29/*29` due to this misrepresentation in a VCF file.


## Problems with whole-exome sequencing/low coverage whole-genome sequencing 

**We strongly discourage using whole-exome sequencing or low coverage whole-genome sequencing data to call CYP2D6** based on our comparison below.

For demonstration purposes, we are using [StellarPGx](https://github.com/SBIMB/StellarPGx), which supports calling CYP2D6 with structural variants and CNV from whole genome sequencing (WGS) CRAM/BAM files. StellarPGx calls have been compared against the CDC's [GeT-RM](https://www.cdc.gov/labquality/get-rm/inherited-genetic-diseases-pharmacogenetics/pharmacogenetics.html) dataset which provides consensus reference diplotypes for selected samples. According to the [StellarPGx paper](https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/cpt.2173), StellarPGx produces results with a higher concordance to GeT-RM than other available CYP2D6 callers. However, this benchmark was performed using GeT-RM's short-read 30x PCR-free WGS dataset. In addition to high-coverage WGS, GeT-RM also provides access to low-coverage WGS and whole exome sequencing (WES) data for some samples. We selected 6 random Get-RM samples and ran StellarPGx on the high-coverage WGS, low-coverage WGS, and WES and compared the results to the GeT-RM ground truth calls.

| Sample | WES | Low-Cov WGS | 30x WGS | GeT-RM |
| ------ | --- | ----------- | ------- | ------ |
| [HG00436](https://www.internationalgenome.org/data-portal/sample/HG00436) | No_call | No_call | *71/*2x2 | *2x2/*71 |
| [HG01086](https://www.internationalgenome.org/data-portal/sample/HG01086) | *1/*1x2 | No_call | *1/*31 | *1/*31 |
| [HG01190](https://www.internationalgenome.org/data-portal/sample/HG01190) | *1/*1 | *5/*1 | *5/*68+*4 | *68+*4/*5 |
| [NA07048](https://www.internationalgenome.org/data-portal/sample/NA07048) | No_call | *1/*1 | *139/*4 | *1/*4 |
| [NA18545](https://www.internationalgenome.org/data-portal/sample/NA18545) | *34/*34x7 | *1/*1x5 | *36+*10/*36+*10 | *5/*36x2+*10x2 |
| [NA21105](https://www.internationalgenome.org/data-portal/sample/NA21105) | No_call | *34/*34 | *111/*3 | *3/*111 |

There is drastic variation in the calls produced by the three sequencing technologies. In many cases, the WES or low-coverage WGS input resulted in a "no call" result. In some samples, improbably high copy numbers and calls that varied greatly from the GeT-RM calls were produced. The 30x WGS produced identical or equivalent results most of the time, and when it deviated the differences were relatively minor. In conclusion, we do not recommend the use WES or low-coverage WGS for calling CYP2D6 for both research and clinical purposes.


## Working with CYP2D6 in PharmCAT

### Using outside calls

PharmCAT supports pulling in results from other tools using what we call ["outside calls"](/using/Outside-Call-Format).  This allows you to select your preferred CYP2D6 caller and use PharmCAT for the phenotype translation and to include CYP2D6 results in a PharmCAT report.

To include outside calls you need to use the `-po` flag in PharmCAT to specify the outside calls file. For example:

```commandline
# java -jar pharmcat.jar -vcf patient_001.vcf -po patient_001_cyp2d6.txt
```

where `patient_001_cyp2d6.txt` takes the format of:

```text
CYP2D6	*1x2/[*2 + *6]
```

For more information: 

* see [Running PharmCAT](/using/Running-PharmCAT#phenotyper) for details on the `-po` flag
* see [Outside Call Format](/using/Outside-Call-Format) for details on the outside call file


#### Formatting StellarPGx output for PharmCAT

As mentioned above, if you have whole genome sequencing (WGS) CRAM/BAM files, your best option is a tool like StellarPGx.  Please note we are not affiliated with the StellarPGx team and offer no guarantees about its performance. Any questions or concerns on StellarPGx should be directed to the StellarPGx maintainer at [twesigomwedavid@gmail.com](mailto:twesigomwedavid@gmail.com). Usage instructions for StellarPGx can be found on their [GitHub repo](https://github.com/SBIMB/StellarPGx).

While this tutorial is StellarPGx specific, it should illustrate how to integrate genotype calls from other tools.

After running StellarPGx with CYP2D6 as the target gene, it should produce a `<run_name>_summary.txt` file which looks like:

```text
HG00436	*71/*2x2
HG01086	[*1/*31]	Possible novel allele or suballele present: interpret with caution; experimental validation and expert review through PharmVar is recommended
HG01190	*5/*68+*4
NA07048	*139/*4
NA18545	*36+*10/*36+*10
NA21105	*111/*3
```

Column 1 is the sample name and column 2 is the CYP2D6 call. In order to convert this into a PharmCAT-readable format you must split this into a separate file for each sample.

The following is a simple python script that parses StellarPGx's `summary.txt` file and outputs PharmCAT-ready outside calls files for use with the `-a` flag.

```python
import sys
import os

input_path = sys.argv[1] # summary.txt file from StellarPGx
output_dir = sys.argv[2] # output directory

with open(input_path, "r") as infile:
    entries = infile.readlines() # read all the StellarPGx file lines into a list

for entry in entries: # for each sample in the StellarPGx file
    split = entry.split("\t") # split the lines into columns
    sample = split[0] # extract the sample name
    call = split[1].strip("[] \n") # extract the sample call and strip excess spaces and brackets

    if not os.path.exists(output_dir): # create the output directory if it doesn't exist yet
        os.makedir(output_dir)

    with open("%s/%s_cyp2d6.txt" % (output_dir,sample), "w") as outfile:
        # write the sample to a text file in the specified output directory
        outfile.write("CYP2D6\t%s" % call)
```

Usage example:

```console
# python3 stellarPGx_to_PharmCAT.py summary.txt pharmcat_inputs
```

Where `summary.txt` is the StellarPGx output and `pharmcat_inputs` is the directory where the call files should be placed. After running this, calling `ls pharmcat_inputs` to list the directory contents should show:

```console
# ls pharmcat_inputs
HG00436_cyp2d6.txt  HG01086_cyp2d6.txt  HG01190_cyp2d6.txt  NA07048_cyp2d6.txt  NA18545_cyp2d6.txt  NA21105_cyp2d6.txt
```

and `HG00436_cyp2d6.txt` would look like

```text
CYP2D6  *71/*2x2
```

We urge you to manually verify that the calls in the PharmCAT-ready files correspond to the StellarPGx output.



### Calling CYP2D6 with PharmCAT's research mode

As stated above, we **do NOT recommend calling CYP2D6 from VCF** due to the large influence of SV and CNV on phenotype prediction.

PharmCAT can, however, call CYP2D6 star alleles that are defined based on SNPs and/or INDELs. For more information, read the [Research Mode](/using/Research-Mode) page.
