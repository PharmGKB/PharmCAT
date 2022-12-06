---
parent: FAQs
title: PharmCAT Runtime
permalink: Runtime/
---

# PharmCAT Runtime

This page documents the records of PharmCAT runtime on different datasets.


### UK Biobank 200K Integrated Call Set

Case 1: the VCF Preprocessor using 20 processors
- Job: Running the VCF Preprocessor
- Sample Size = 200,044
- Nodes: 1
- Used builtin multiprocessing support `-c`: Yes
- Cores/processors per node: 20
  - A modern CPU is composed of numerous cores, typically 10 to 36.
  - The VCF Preprocessor detected 31 usable processors for concurrent processing.
- Elapsed time: 38 hr 12 min 19 sec
- Maximum memory utilized: 2.98 GB
- **Average speed** = 0.68 seconds/sample
- **Overall time** = 38 hours for 200K samples using 20 processors

Case 2: PharmCAT
- Job: Running PharmCAT
- Sample Size = 200,044
- Nodes: 1
- Used builtin multiprocessing support: No
- Used a parallelization framework: Yes
- Cores/processors per node: 100
  - 200K samples were divided into 100 subsets.
  - Each of the 100 subsets was run in parallel on Stanford Sherlock HPC system.
  - Within each subset, samples were run sequentially.
- Elapsed time: 3 hr 17 min 00 sec
  - This was the time when the analyses on all subsets were finished
- Maximum memory utilized: 195.85 MB
- **Average speed** = 5.9 seconds/sample
- **Overall time** = 3 hours for 200K samples by running 100 parallel subsets


### Penn Medicine Biobank

Case 1: the VCF Preprocessor on 43K samples
- Job: Running the VCF Preprocessor
- Sample Size = 43K
- Nodes: 1
- Cores/processors per node: 1
- Used builtin multiprocessing support `-c`: No
- Used a parallelization framework: No
- Elapsed time: 212196.12 sec
- Maximum memory utilized: ~2.1 GB
- **Average speed** = ~4.9 seconds/sample
- **Overall time** = 59 hours for 43K samples using 1 processors

Case 2: PharmCAT on 43K samples
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
- Maximum memory utilized: ~200 MB
- **Average speed** = ~5.6 seconds/sample
- **Overall time** = 8 minutes for 43K samples by running 504 parallel subsets


### All of Us

Case 1: All on 100K WGS
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
- Maximum memory utilized: Unknown
- **Average speed** = 0.55-0.58 seconds/sample
- **Overall time** = 16 hours for 100K samples using 96 processors



