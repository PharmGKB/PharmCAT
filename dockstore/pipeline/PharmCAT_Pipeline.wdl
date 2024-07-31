version 1.0

# It is a single task that runs the PharmCAT pipeline on a VCF file
# as a single task the workflow is just a wrapper for the task and get the same name

# The output is an array of files that are the results of the pipeline save in the plataform where the workflow is running

task pharmcat_pipeline {
    input {
        File vcf_file
        String sample_ids = ""
        File? sample_file
        Boolean missing_to_ref = false
        Boolean no_gvcf_check = false
        Boolean retain_specific_regions = false
        File? reference_regions
        Boolean run_matcher = false
        Boolean matcher_all_results = false
        Boolean matcher_save_html = false
        String research_mode = ""
        Boolean run_phenotyper = false
        Boolean run_reporter = false
        String reporter_sources = ""
        Boolean reporter_extended = false
        Boolean reporter_save_json = false
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
        ~{if no_gvcf_check then '-G' else ''} \
        ~{if retain_specific_regions then '-R' else ''} \
        ~{if defined(reference_regions) then '-refRegion ' + reference_regions else ''} \
        ~{if run_matcher then '-matcher' else ''} \
        ~{if matcher_all_results then '-ma' else ''} \
        ~{if matcher_save_html then '-matcherHtml' else ''} \
        ~{if research_mode != "" then '-research ' + research_mode else ''} \
        ~{if run_phenotyper then '-phenotyper' else ''} \
        ~{if run_reporter then '-reporter' else ''} \
        ~{if reporter_sources != "" then '-rs ' + reporter_sources else ''} \
        ~{if reporter_extended then '-re' else ''} \
        ~{if reporter_save_json then '-reporterJson' else ''} \
        ~{if base_filename != "" then '-bf ' + base_filename else ''} \
        ~{if delete_intermediate_files then '-del' else ''} \
        -cp ~{max_concurrent_processes} -cm ~{max_memory}
    >>>

    output {
        Array[File] results = glob("data/*")
    }

    runtime {
        docker: "pgkb/pharmcat:2.13.0"
        memory: max_memory
        cpu: max_concurrent_processes
    }

    meta {
        author: "PharmCAT Team"
        email: "pharmcat@pharmgkb.org"
        description: "Workflow to run the PharmCAT pipeline on a VCF file"
    }
}

workflow pharmcat_pipeline {
    input {
        File vcf_file
        String sample_ids = ""
        File? sample_file
        Boolean missing_to_ref = false
        Boolean no_gvcf_check = false
        Boolean retain_specific_regions = false
        File? reference_regions
        Boolean run_matcher = false
        Boolean matcher_all_results = false
        Boolean matcher_save_html = false
        String research_mode = ""
        Boolean run_phenotyper = false
        Boolean run_reporter = false
        String reporter_sources = ""
        Boolean reporter_extended = false
        Boolean reporter_save_json = false
        String base_filename = ""
        Boolean delete_intermediate_files = false
        Int max_concurrent_processes = 1
        String max_memory = "4G"
    }

    call pharmcat_pipeline {
        input:
            vcf_file = vcf_file,
            sample_ids = sample_ids,
            sample_file = sample_file,
            missing_to_ref = missing_to_ref,
            no_gvcf_check = no_gvcf_check,
            retain_specific_regions = retain_specific_regions,
            reference_regions = reference_regions,
            run_matcher = run_matcher,
            matcher_all_results = matcher_all_results,
            matcher_save_html = matcher_save_html,
            research_mode = research_mode,
            run_phenotyper = run_phenotyper,
            run_reporter = run_reporter,
            reporter_sources = reporter_sources,
            reporter_extended = reporter_extended,
            reporter_save_json = reporter_save_json,
            base_filename = base_filename,
            delete_intermediate_files = delete_intermediate_files,
            max_concurrent_processes = max_concurrent_processes,
            max_memory = max_memory
    }

    output {
        Array[File] results_all = pharmcat_pipeline.results
    }
}

