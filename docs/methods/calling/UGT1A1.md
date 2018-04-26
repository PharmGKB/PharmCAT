---
title: Calling UGT1A1
permalink: methods/calling/UGT1A1/
---

# Calling and Reporting on UGT1A1/atazanavir in PharmCAT

When getting returns for UGT1A1 from the caller, the reporter needs to
be able to determine if the return is from phased or unphased
information.

**Treat \*80+\*28 as \*28 and \*80+\*37 as \*37**


## If from phased

1. Then print the return as the diplotype
2. For the recommendation, metabolizer status:

* If there is only 1 \* allele per allele (1 \* allele on each side of the "/")
  * Then chose picker language for those 2 alleles  
  AND choose metabolizer status from that combo in the picker
* If there are \> 1 \* allele on either side of the "/"
  * Then use logic for unphased data below

## If from unphased

1. Then print a list of all \* alleles called with the copy number OR
print "heterozygous" for all alleles with 1 copy and "homozygous" for
all alleles with 2 copies in the diplotype field
2. For the recommendation:
  * If \[homozygous (2 copies) of \>=1 of \[\*28, \*6, \*37, \*27\]
  * OR If \[heterozygous (1 copy) for \>= 1 \[\*28\] AND \[heterozygous (1 copy) for \>=1 of \[\*6 , \*37, \*27\] \]
  * OR If \[heterozygous (1 copy) for \>=2 of \[\*6 , \*37, \*27\]
    * Then chose picker language for \*80/\*80  
    AND use "Poor metabolizer" as metabolizer status

  * If \[\[heterozygous (1 copy) of \>=1 of \[\*28\] OR =1 of \[\*6, \*37, \*27\]  
  AND No (0) copies for EVERY OTHER \* allele out of \[\*28, \*6, \*37, \*27\]
    * Then chose picker language for \*1/\*80  
    AND use "Intermediate metabolizer" as metabolizer status \]

  * If \[0 copies of \*6 AND \*37 AND \*27 AND \*28\]
    * Then chose picker language for \*1/\*1  
    AND use "Normal metabolizer" as metabolizer status \]

For all calls, phased or not, in Allele Function column of report, list
the functionality for every \* allele returned and go off of the
spreadsheet for that.
