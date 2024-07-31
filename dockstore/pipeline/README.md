# WDL to run PharmCAT_Pipeline

This WDL script executes the PharmCAT pipeline on a specified VCF file or a set of VCF files, processing genetic data to provide pharmacogenomic insights. The workflow automates the execution of the PharmCAT pipeline, streamlining the analysis of genetic variants to predict drug response and tailor medical treatment to individual patients' genetic profiles. By leveraging the Workflow Description Language (WDL), this script ensures reproducibility, scalability, and ease of use across various computational environments.


## Input Parameters

### Input Arguments
- `File vcf_file`: Path to a VCF file or a file of paths to VCF files (one file per line), sorted by chromosome position.
- `String sample_ids` (default: `""`): A comma-separated list of sample IDs.
- `File? sample_file` (default: `null`): A file containing a list of samples, one sample per line.

### Preprocessor Arguments
- `Boolean missing_to_ref` (default: `false`): Assume genotypes at missing PGx sites are 0/0. DANGEROUS!.
- `Boolean no_gvcf_check` (default: `false`): Bypass the gVCF check for the input VCF. DANGEROUS!.
- `Boolean retain_specific_regions` (default: `false`): Retain the genomic regions specified by `-refRegion`.
- `File? reference_regions` (default: `null`): A sorted bed file of specific PGx regions to retain. Must be used with the `-R` argument.

### Named Allele Matcher Arguments
- `Boolean run_matcher` (default: `false`): Run named allele matcher independently.
- `Boolean matcher_all_results` (default: `false`): Return all possible diplotypes, not just top hits.
- `Boolean matcher_save_html` (default: `false`): Save named allele matcher results as HTML.
- `String research_mode` (default: `""`): Comma-separated list of research features to enable: [cyp2d6, combinations].

### Phenotyper Arguments
- `Boolean run_phenotyper` (default: `false`): Run phenotyper independently.

### Reporter Arguments
- `Boolean run_reporter` (default: `false`): Run reporter independently.
- `String reporter_sources` (default: `""`): Comma-separated list of sources to limit report to: [CPIC, DPWG].
- `Boolean reporter_extended` (default: `false`): Output extended report.
- `Boolean reporter_save_json` (default: `false`): Save reporter results as JSON.

### Output Arguments
- `String base_filename` (default: `""`): Prefix for output files. Defaults to the same base name as the input.
- `Boolean delete_intermediate_files` (default: `false`): Delete intermediate PharmCAT files (saved by default).

### Concurrency/Memory Arguments
- `Int max_concurrent_processes` (default: `1`): The maximum number of processes to use when concurrent mode is enabled.
- `String max_memory` (default: `"4G"`): The maximum memory PharmCAT should use (e.g. "64G"). This is passed to Java using the -Xmx flag.

## Outputs
- `Array[File] results_all`: The results of the PharmCAT pipeline. These files are saved in the execution directory of the job.

## Running the PharmCAT Pipeline
For convenience, the `pharmcat_pipeline script simplifies the process of running the entire PharmCAT pipeline (the VCF Preprocessor and the core PharmCAT tool). The necessary dependencies are already included in the provided image.

### Prerequisites
The required dependencies (python3, java, bcftools, bgzip) are included in the provided image, so no additional installation is needed.

### Execution Platforms
This WDL script can be executed using [Cromwell](https://github.com/broadinstitute/cromwell) or on platforms such as [Terra](https://support.terra.bio/hc/en-us) and [AnVIL](https://anvil.terra.bio/#) that can be launched from [Dockstore](https://dockstore.org).

### Local Execution
To run the WDL with the PharmCAT-Pipeline locally, ensure that Docker and Cromwell are installed in your execution environment. Then, execute the following command in your bash:

```sh
$ java -jar {path}/cromwell-{version}.jar run {path}/pharmCAT_Pipeline.wdl -i {path}/inputs.json
```

Here is an example of how to provide the inputs in a JSON file:

```json
{
  "pharmcat_pipeline.vcf_file": "gs://your-bucket/path/to/your.vcf",
  "pharmcat_pipeline.sample_ids": "",
  "pharmcat_pipeline.sample_file": null,
  "pharmcat_pipeline.missing_to_ref": false,
  "pharmcat_pipeline.no_gvcf_check": false,
  "pharmcat_pipeline.retain_specific_regions": false,
  "pharmcat_pipeline.reference_regions": null,
  "pharmcat_pipeline.run_matcher": false,
  "pharmcat_pipeline.matcher_all_results": false,
  "pharmcat_pipeline.matcher_save_html": false,
  "pharmcat_pipeline.research_mode": "",
  "pharmcat_pipeline.run_phenotyper": false,
  "pharmcat_pipeline.run_reporter": false,
  "pharmcat_pipeline.reporter_sources": "",
  "pharmcat_pipeline.reporter_extended": false,
  "pharmcat_pipeline.reporter_save_json": false,
  "pharmcat_pipeline.base_filename": "",
  "pharmcat_pipeline.delete_intermediate_files": false,
  "pharmcat_pipeline.max_concurrent_processes": 1,
  "pharmcat_pipeline.max_memory": "4G"
}
```


## Software Version Notes

### WDL Version
This script is written in WDL (Workflow Description Language) 1.0 to ensure compatibility with Dockstore. For more information about WDL, visit the WDL Documentation.

### Cromwell Version
Successfully tested on v53 and v87

### PharmCAT and Imagem Version
v:2.13.0

### Important Notes
- Runtime parameters are optimized for Google Cloud Platform implementation.
- The provided JSON is a generic ready-to-use example template for the workflow. It is the user’s responsibility to correctly set the reference and resource variables for their own particular test case using the PharmCAT documentation.


## Documentation Links
- [PharmCAT Documentation](https://pharmcat.org/)
- [PharmCAT-Pipeline Documentation](https://pharmcat.org/using/Running-PharmCAT-Pipeline/)
- [PharmCAT Project](https://github.com/PharmGKB/PharmCAT)

## Contact
For technical questions or bug reports, [file an issue](https://github.com/PharmGKB/PharmCAT/issues).

For general questions about the PharmCAT project, contact [pharmcat@pharmgkb.org](mailto:pharmcat@pharmgkb.org).


## License
(C) 2024 Your Organization | BSD-3

This script is released under the [WDL open source code license](https://github.com/openwdl/wdl/blob/master/LICENSE) (BSD-3). Note, however, that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script.