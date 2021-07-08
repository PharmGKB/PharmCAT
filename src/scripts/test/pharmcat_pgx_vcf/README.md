## Overview

This folder contains the reference PGx allele defining positions used by PharmCAT in VCF format. As we are currently going through the prelease, there can be some variant representation errors in the PharmCAT VCF.

## Explanation of files

1. `pharmcat_positions_0.8.0.vcf`. This is the reference VCF that one can download from the PharmCAT GitHub release.
   
2. `pharmcat_positions_0.8.0_updated_06222021.vcf`. This is the reference VCF that we recommend PharmCAT users to use. This file has gone through bcftools normalization and manual inspection. This file is normalized by the bcftools command:
   
	```bcftools norm -m+ -c ws -Ov -o pharmcat_positions_0.8.0_updated_06222021.vcf -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna pharmcat_positions_0.8.0.vcf```

    Notes:
        - The bcftools norm did not correctly normalize rs746071566 from the `pharmcat_positions_0.8.0.vcf` due to the complex variant representation formats extracted from the source databases (PharmVar). The error has been manually corrected. 
        - Other positions passed the manual inspection.





