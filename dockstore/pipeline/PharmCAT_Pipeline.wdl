version 1.0

workflow pharmcat_pipeline {
  meta {
    author: "ClinPGx"
    email: "pharmcat@clinpgx.org"
    description: "This workflow runs a VCF file through the PharmCAT pipeline."
  }

  parameter_meta {
    # description for this is intentionally different from pipeline script because it's hard to
    # support a file of files on cloud services and directories aren't supported
    vcf_file: "A VCF file (can be gzipped or bgzipped)."
    sample_ids: "A comma-separated list of sample IDs.  Only applicable if you have multiple samples and only want to work on specific ones."
    sample_file: "A file containing a list of sample IDs, one sample ID per line.  Only applicable if you have multiple samples and only want to work on specific ones."

    missing_to_ref: 'Assume genotypes at absent or unspecified PGx sites are "0/0".  DANGEROUS!'
    absent_to_ref: 'Assume genotypes at absent PGx sites are "0/0".  DANGEROUS!'
    unspecified_to_ref: 'Assume unspecified genotypes ("./.") are "0/0" when every sample is "./.". DANGEROUS!'
    no_gvcf_check: "Bypass check if VCF file is in gVCF format."
    # not including retain_specific_regions and reference_regions

    run_matcher: "Run named allele matcher independently."
    matcher_all_results: "Return all possible diplotypes, not just top hits."
    matcher_save_html: "Save named allele matcher results as HTML.'"
    research_mode: "Comma-separated list of research features to enable: [cyp2d6, combinations]"

    run_phenotyper: "Run phenotyper independently."

    run_reporter: "Run reporter independently."
    reporter_sources: "Comma-separated list of sources to limit recommendations to: [CPIC, DPWG, FDA]"
    reporter_extended: "Write an extended report (includes all possible genes and drugs, even if no data is available)"
    reporter_save_html: "Save reporter results as HTML (the default if no format is specified)."
    reporter_save_json: "Save reporter results as JSON."
    reporter_save_calls_only_tsv: "Save call results only as TSV."

    base_filename: "Prefix for output files.  Defaults to the same base name as the input."
    delete_intermediate_files: "Delete intermediate PharmCAT files.  Defaults to saving all files."

    max_concurrent_processes: "The maximum number of processes to use when concurrent mode is enabled."
    max_memory: "The maximum memory PharmCAT should use (e.g. '64G')."
  }


  input {
    File vcf_file
    String sample_ids = ""
    File? sample_file
    Boolean missing_to_ref = false
    Boolean absent_to_ref = false
    Boolean unspecified_to_ref = false
    Boolean no_gvcf_check = false
    Boolean run_matcher = false
    Boolean matcher_all_results = false
    Boolean matcher_save_html = false
    String research_mode = ""
    Boolean run_phenotyper = false
    Boolean run_reporter = false
    String reporter_sources = ""
    Boolean reporter_extended = false
    Boolean reporter_save_html = false
    Boolean reporter_save_json = false
    Boolean reporter_save_calls_only_tsv = false
    String base_filename = ""
    Boolean delete_intermediate_files = false
    Int max_concurrent_processes = 1
    String max_memory = "4G"
  }

  call pharmcat_pipeline_task {
    input:
      vcf_file = vcf_file,
      sample_ids = sample_ids,
      sample_file = sample_file,
      missing_to_ref = missing_to_ref,
      absent_to_ref = absent_to_ref,
      unspecified_to_ref = unspecified_to_ref,
      no_gvcf_check = no_gvcf_check,
      run_matcher = run_matcher,
      matcher_all_results = matcher_all_results,
      matcher_save_html = matcher_save_html,
      research_mode = research_mode,
      run_phenotyper = run_phenotyper,
      run_reporter = run_reporter,
      reporter_sources = reporter_sources,
      reporter_extended = reporter_extended,
      reporter_save_html = reporter_save_html,
      reporter_save_json = reporter_save_json,
      reporter_save_calls_only_tsv = reporter_save_calls_only_tsv,
      base_filename = base_filename,
      delete_intermediate_files = delete_intermediate_files,
      max_concurrent_processes = max_concurrent_processes,
      max_memory = max_memory
  }

  output {
    Array[File] results = pharmcat_pipeline_task.results
  }
}


task pharmcat_pipeline_task {
  meta {
    author: "ClinPGx"
    email: "pharmcat@clinpgx.org"
    description: "This task run a VCF file through the PharmCAT pipeline."
  }

  input {
    File vcf_file
    String sample_ids = ""
    File? sample_file
    Boolean missing_to_ref = false
    Boolean absent_to_ref = false
    Boolean unspecified_to_ref = false
    Boolean no_gvcf_check = false
    Boolean run_matcher = false
    Boolean matcher_all_results = false
    Boolean matcher_save_html = false
    String research_mode = ""
    Boolean run_phenotyper = false
    Boolean run_reporter = false
    String reporter_sources = ""
    Boolean reporter_extended = false
    Boolean reporter_save_html = false
    Boolean reporter_save_json = false
    Boolean reporter_save_calls_only_tsv = false
    String base_filename = ""
    Boolean delete_intermediate_files = false
    Int max_concurrent_processes = 1
    String max_memory = "4G"
  }

  command <<<
    set -x -e -o pipefail
    mkdir -p data
    cp ~{vcf_file} data/

    pharmcat_pipeline data/$(basename ~{vcf_file}) \
    ~{if sample_ids != "" then '-s ' + sample_ids else ''} \
    ~{if defined(sample_file) then '-S ' + sample_file else ''} \
    ~{if missing_to_ref then '-0' else ''} \
    ~{if absent_to_ref then '--absent-to-ref' else ''} \
    ~{if unspecified_to_ref then '--unspecified-to-ref' else ''} \
    ~{if no_gvcf_check then '-G' else ''} \
    ~{if run_matcher then '-matcher' else ''} \
    ~{if matcher_all_results then '-ma' else ''} \
    ~{if matcher_save_html then '-matcherHtml' else ''} \
    ~{if research_mode != "" then '-research ' + research_mode else ''} \
    ~{if run_phenotyper then '-phenotyper' else ''} \
    ~{if run_reporter then '-reporter' else ''} \
    ~{if reporter_sources != "" then '-rs ' + reporter_sources else ''} \
    ~{if reporter_extended then '-re' else ''} \
    ~{if reporter_save_html then '-reporterHtml' else ''} \
    ~{if reporter_save_json then '-reporterJson' else ''} \
    ~{if reporter_save_calls_only_tsv then 'reporterCallsOnlyTsv' else ''} \
    ~{if base_filename != "" then '-bf ' + base_filename else ''} \
    ~{if delete_intermediate_files then '-del' else ''} \
    -cp ~{max_concurrent_processes} -cm ~{max_memory}
  >>>

  output {
    Array[File] results = glob("data/*")
  }

  runtime {
    docker: "pgkb/pharmcat:3.1.1"
    memory: max_memory
    cpu: max_concurrent_processes
  }
}

