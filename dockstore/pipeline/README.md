# WDL to run the PharmCAT pipeline

PharmCAT (Pharmacogenomics Clinical Annotation Tool) is a bioinformatics tool that analyzes genetic variants to predict
drug response and tailor medical treatment to an individual patient’s genetic profile. It does this in two phases:

1. Processes VCF files from next generation sequencing (NGS) or genotyping methods and identifies pharmacogenomic (PGx)
   genotypes and infers haplotypes, typically called star alleles.
2. Uses the pharmacogene diplotypes (combination of maternal and paternal star alleles) to predict PGx phenotypes and
   reports the corresponding drug-prescribing recommendations from [CPIC guidelines](https://cpicpgx.org/guidelines/),
   [PharmGKB-annotated DPWG guidelines](https://www.pharmgkb.org/page/dpwg) and
   [PharmGKB-annotated FDA-approved drug labels](https://www.pharmgkb.org/page/drugLabelLegend).

This WDL script executes the PharmCAT pipeline, taking VCF files as input and produces a recommendation report for each
sample. By leveraging the Workflow Description Language (WDL), this script ensures reproducibility, scalability, and
ease of use across various computational environments.

For details, see:

- [PharmCAT Pipeline Documentation](https://pharmcat.org/using/Running-PharmCAT-Pipeline/)
- [PharmCAT Documentation](https://pharmcat.org/)
- [PharmCAT Repository](https://github.com/PharmGKB/PharmCAT)



## Input Parameters

The only required input is a VCF file.
An example VCF file you can use to test with can be found [here](https://pharmcat.org/examples/pharmcat.example.vcf). 


### Input Arguments
- `File vcf_file`: Path to a VCF file or a directory containing VCF files.
- `String sample_ids` (default: `""`): A comma-separated list of sample IDs. Only applicable if you have multiple samples and only want to work on specific ones.
- `File? sample_file` (default: `null`): A file containing a list of samples, one sample per line. Only applicable if you have multiple samples and only want to work on specific ones.

### Preprocessor Arguments
- `Boolean missing_to_ref` (default: `false`): Assume genotypes at absent or unspecified PGx sites are "0/0".  DANGEROUS!
   Equivalent to using both `absent_to_ref` and `unspecified_to_ref`
- `Boolean absent_to_ref` (default: `false`): Assume genotypes at absent PGx sites are "0/0".  DANGEROUS!
- `Boolean unspecified_to_ref` (default: `false`): Assume unspecified genotypes ("./.") are "0/0" when every sample is "./.". DANGEROUS!
- `Boolean no_gvcf_check` (default: `false`): Bypass the gVCF check for the input VCF.

### Named Allele Matcher Arguments
- `Boolean run_matcher` (default: `false`): Run named allele matcher independently.
- `Boolean matcher_all_results` (default: `false`): Return all possible diplotypes, not just top hits.
- `Boolean matcher_save_html` (default: `false`): Save named allele matcher results as HTML.
- `String research_mode` (default: `""`): Comma-separated list of research features to enable: [cyp2d6, combinations].

### Phenotyper Arguments
- `Boolean run_phenotyper` (default: `false`): Run phenotyper independently.

### Reporter Arguments
- `Boolean run_reporter` (default: `false`): Run reporter independently.
- `String reporter_sources` (default: `""`): Comma-separated list of sources to limit recommendations to: [CPIC, DPWG, FDA].
- `Boolean reporter_extended` (default: `false`): Write an extended report (includes all possible genes and drugs, even if no data is available)
- `Boolean reporter_save_json` (default: `false`): Save reporter results as JSON.

### Output Arguments
- `String base_filename` (default: `""`): Prefix for output files. Defaults to the same base name as the input.
- `Boolean delete_intermediate_files` (default: `false`): Delete intermediate PharmCAT files. Defaults to saving all files.

### Concurrency/Memory Arguments
- `Int max_concurrent_processes` (default: `1`): The maximum number of processes to use when concurrent mode is enabled.
- `String max_memory` (default: `"4G"`): The maximum memory PharmCAT should use (e.g. "64G"). This is passed to Java using the -Xmx flag.


## Outputs
- `Array[File] results_all`: The results of the PharmCAT pipeline. These files are saved in the execution directory of the job.


## Execution Platforms
This WDL script can be executed using [Cromwell](https://github.com/broadinstitute/cromwell) or on platforms such as
[Terra](https://support.terra.bio/hc/en-us) and [AnVIL](https://anvil.terra.bio/) that can be launched from
[Dockstore](https://dockstore.org).


### Local Execution
To run the WDL with the PharmCAT-Pipeline locally, ensure that Docker and Cromwell are installed in your execution
environment.
Then, execute the following command in your shell:

```sh
$ java -jar {path}/cromwell-{version}.jar run {path}/PharmCAT_Pipeline.wdl -i {path}/inputs.json
```

Here is an example of how to provide the inputs in a JSON file:

```json
{
  "pharmcat_pipeline.vcf_file": "gs://your-bucket/path/to/your.vcf",
  "pharmcat_pipeline.sample_ids": "",
  "pharmcat_pipeline.sample_file": null,
  "pharmcat_pipeline.missing_to_ref": false,
  "pharmcat_pipeline.absent_to_ref": false,
  "pharmcat_pipeline.unspecified_to_ref": false,
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
  "pharmcat_pipeline.reporter_save_html": true,
  "pharmcat_pipeline.reporter_save_json": false,
  "pharmcat_pipeline.reporter_save_calls_only_tsv": false,
  "pharmcat_pipeline.base_filename": "",
  "pharmcat_pipeline.delete_intermediate_files": false,
  "pharmcat_pipeline.max_concurrent_processes": 1,
  "pharmcat_pipeline.max_memory": "4G"
}
```


## Important Notes
- Runtime parameters are optimized for Google Cloud Platform implementation.
- The provided JSON is a generic ready-to-use example template for the workflow.
  It is the user’s responsibility to correctly set the reference and resource variables for their own particular data
  using the PharmCAT documentation.


## Version Information

### WDL Version
This script is written in WDL (Workflow Description Language) [1.0](https://docs.openwdl.org/en/1.0.0/) to ensure
compatibility with Dockstore.
For more information about WDL, visit the [WDL website](https://openwdl.org/).

### Cromwell Version
Successfully tested on v53 and v87.

### PharmCAT Version

PharmCAT v3.0.0.


## Contact
For technical questions or bug reports, [file an issue](https://github.com/PharmGKB/PharmCAT/issues).

For general questions about the PharmCAT project, contact [pharmcat@pharmgkb.org](mailto:pharmcat@pharmgkb.org).


## License
(C) 2024 PharmGKB | MPL-2

This script, and PharmCAT, is released under the [Mozilla Public License 2.0](https://www.mozilla.org/en-US/MPL/2.0/).
Note, however, that the programs it calls may be subject to different licenses.
Users are responsible for checking that they are authorized to run all programs before running this script.
