---
parent: Working with Large Datasets
title: Run Time and Cost
permalink: using/Runtime-Cost/
nav_order: 2
---
# PharmCAT Run Time and Cost
{: .no_toc }

This page records PharmCAT run times on different datasets and the estimates of computing resource cost
on different cloud-based research analytic platforms.

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

---


## UK Biobank 200K Integrated Call Set

These statistics were tested on the Stanford Sherlock high-performance computing center by Binglan Li from PharmGKB
at Stanford University.


### Case 1. The VCF Preprocessor

This was documented before the VCF Preprocessor supported multi-sample VCF by default.

#### Case 1.1. Using 31 processors

- Job: Running the VCF Preprocessor
  - Date: Dec 3, 2022
  - command: `python3 "$VCF_PREPROCESS_SCRIPT" -vcf "$INPUT_VCF" -refFna "$REF_SEQ" -refVcf "$REF_PGX_VCF" -o "$PREPROCESSED_VCF_OUTPUT_DIR"/ -c`
- Sample Size = 200,044
- Nodes: 1
- Utilized cores/processors: 31
  - A modern CPU is composed of many cores, typically 10 to 36.
- Elapsed time: 38 hr 12 min 19 sec
- Maximum memory used: 2.98 GB
- **Average speed** = 0.68 seconds/sample
- **Overall time** = 38 hours for 200K samples using 20 processors
- Cost estimate on the UK Biobank Research Analytic Platform using Swiss Army Knife
  - Instance requested: 70.3 GB total memory, 3600 GB total storage, 36 cores
  - Estimated Cost Per Hour: £0.4464
  - Suppose the whole job takes 38 hrs to finish
  - Total cost: ~ £17 or $18
- Cost estimate on the Google Cloud for All of Us
  - $166 per month for running 2 days per week using one instance of e2-highcpu-32 (vCPUs: 32, RAM: 32GB)
  - Conservatively speaking, for running one job, it might cost $50

#### Case 1.2. Using 132 processors

- Job: Running the VCF Preprocessor
    - Date: Dec 6, 2022
    - command: `python3 "$VCF_PREPROCESS_SCRIPT" -vcf "$INPUT_VCF" -refFna "$REF_SEQ" -refVcf "$REF_PGX_VCF" -o "$PREPROCESSED_VCF_OUTPUT_DIR"/ -c`
- Sample Size = 200,044
- Nodes: 1
- Utilized cores/processors: 132
- Elapsed time: 2 hr 10 min 04 sec
- Maximum memory used: 1.69 GB
- **Average speed** = 0.68 seconds/sample
- **Overall time** = 2 hours for 200K samples using 132 processors
- Cost estimate on the UK Biobank Research Analytic Platform using Swiss Army Knife
    - Instance requested: 187.5 GB total memory, 9600 GB total storage, 96 cores
    - Estimated Cost Per Hour: £1.1904
    - To make a conservative educated guess, suppose the whole job takes 5 hrs to finish
    - Total cost: ~ £6 or $6.36


### Case 2. the VCF Preprocessor

This was documented after the VCF Preprocessor supported multi-sample VCF by default.

- Job: Running the VCF Preprocessor
    - Date: Dec 6, 2022
    - command: `python3 "$VCF_PREPROCESS_SCRIPT" -vcf "$INPUT_VCF" -refFna "$REF_SEQ" -refVcf "$REF_PGX_VCF" -o "$PREPROCESSED_VCF_OUTPUT_DIR"/ -c`
- Sample Size = 200,044
- Nodes: 1
- Utilized cores/processors: 23
- Elapsed time: 2 hr 49 min 56 sec
- Maximum memory used: 1.56 GB
- **Average speed** = 0.05 seconds/sample
- **Overall time** = 3 hours for 200K samples using 23 processors (compare with case 1.1)
- Cost estimate on the UK Biobank Research Analytic Platform
  - Instance requested: 70.3 GB total memory, 3600 GB total storage, 36 cores
  - Estimated Cost Per Hour: £0.4464
  - To make a conservative educated guess, suppose the whole job takes 3 hrs to finish
  - Total cost: ~ £1.3


### Case 3. PharmCAT - 100 subsets

#### Case 3.1. Running 100 subsets in parallel

- Job: Running PharmCAT
- Sample Size = 200,044
- Nodes: 1
- Used builtin multiprocessing support: No
- Used a parallelization framework: Yes
- Cores/processors per node: 1
  - 200K samples were divided into 100 subsets.
  - Each of the 100 subsets was run using 1 processor in parallel on Stanford Sherlock HPC system.
  - Within each subset, samples were run sequentially.
- Elapsed time: 3 hr 17 min 00 sec
  - This was the time when the analyses on all subsets were finished
- Maximum memory used: 195.85 MB
- **Average speed** = 5.9 seconds/sample
- **Overall time** = 3 hours for 200K samples by running 100 parallel subsets
- Cost estimate on the UK Biobank Research Analytic Platform
    - Instance requested: 3.9 GB total memory, 200 GB total storage, 2 cores
    - Estimated Cost Per Hour: £0.0248
    - To make a conservative educated guess, suppose each subset job takes 4 hrs to finish
    - Total cost: ~ £10 or $10.60

#### Case 3.2. Running 996 subsets in parallel

- Job: Running PharmCAT
- Sample Size = 200,044
- Nodes: 1
- Used builtin multiprocessing support: No
- Used a parallelization framework: Yes
- Cores/processors: 1
    - 200K samples were divided into 996 subsets.
    - Each of the 996 subsets was run using 1 processor in parallel on Stanford Sherlock HPC system.
    - Within each subset, samples were run sequentially.
- Elapsed time: 4 min 21 sec
    - This was the time when the analyses on all subsets were finished
- Maximum memory used: 180.44 MB
- **Average speed** = 1.19 seconds/sample
- **Overall time** = 4 mins for 200K samples by running 996 parallel subsets
- Cost estimate on the UK Biobank Research Analytic Platform
    - Instance requested: 3.9 GB total memory, 200 GB total storage, 2 cores
    - Estimated Cost Per Hour: £0.0248
    - It's unclear if the UK Biobank Research Analytic Platform is charged based on hours or minutes
    - To make a conservative educated guess, suppose the platform is charged based on hours
    - Total cost: ~ £24.8 or $26

### Case 4. Extracting PharmCAT JSON results to a TSV file

This was documented before the VCF Preprocessor supported multi-sample VCF by default.

- Job: Extracting JSON content to a TSV file
- Sample Size = 200,044
- Nodes: 1
- Used builtin multiprocessing support: Yes
- Cores/processors: 40
- Elapsed time: 10 hr 9 min 14 sec
- Maximum memory used: 11.60 MB
- **Overall time** = 10hr for 200K samples
- Cost estimate on the UK Biobank Research Analytic Platform
  - Instance requested: 70.3 GB total memory, 3600 GB total storage, 36 cores
  - Estimated Cost Per Hour: £0.4464
  - To make a conservative educated guess, suppose the whole job takes 15 hrs to finish
  - Total cost: ~ £6.70 or $7.10


## Penn Medicine Biobank

The statistics were kindly provided by Karl Keat from Dr. Marylyn Ritchie's group at the University of Pennsylvania.
We thank Karl Keat and Dr. Marylyn Ritchie for their collaboration and contribution to PharmCAT.

### Case 1. the VCF Preprocessor on 43K samples

- Job: Running the VCF Preprocessor
- Sample Size = 43K
- Nodes: 1
- Cores/processors per node: 1
- Used builtin multiprocessing support `-c`: No
- Used a parallelization framework: No
- Elapsed time: 212196.12 sec
- Maximum memory used: ~2.1 GB
- **Average speed** = ~4.9 seconds/sample
- **Overall time** = 59 hours for 43K samples using 1 processors

### Case 2: PharmCAT on 43K samples

- Job: Running PharmCAT
- Sample Size = 43K
- Nodes: 1
- Cores/processors per node: 1
- Used builtin multiprocessing support `-c`: No
- Used a parallelization framework: Yes
  - 43K samples were divided into 504 subsets of 86 samples.
  - Each of the 504 subsets was run in parallel on the Penn HPC system.
  - Within each subset, 86 samples were run sequentially.
- Elapsed time: ~480 seconds
- Maximum memory used: ~200 MB
- **Average speed** = ~5.6 seconds/sample
- **Overall time** = 8 minutes for 43K samples by running 504 parallel subsets


## All of Us

The statistics were kindly provided by Andrew Haddad from Dr. Philip Empey's at the University of Pittsburgh.
We thank Andrew Haddad and Dr. Philip Empey for their collaboration and contribution to PharmCAT.

### Case 1. All on 100K WGS

- Job: Running all analyses on 100K WGS _All of Us_ data
  - This includes
  - (1) copying files from permanent storage on Google buckets to computing nodes
  - (2) running the VCF Preprocessor
  - (3) running PharmCAT
  - (4) organizing results from JSON files to a summary TSV file
- Sample Size = 98,622
- Nodes: 1
- CPUs per node: 96
- Used builtin multiprocessing support `-c`: No
- Used a parallelization framework: Yes
- Elapsed time: Approximately 15-16 hrs
- Maximum memory used: Unknown
- **Average speed** = 0.55-0.58 seconds/sample
- **Overall time** = 16 hours for 100K samples using 96 processors
- Cost estimate for All of Us
  - $100 per run of PharmCAT on 100K samples with a VM with 96 CPUs.

