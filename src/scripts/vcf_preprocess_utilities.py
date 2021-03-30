#! /usr/bin/env python
__author__ = 'BinglanLi'

import os
import re
import sys
import subprocess
import shutil
import urllib.parse
import urllib.request
import gzip
import vcf_preprocess_exceptions as Exceptions
import tempfile
from cyvcf2 import VCF, Writer
import allel

def obtain_vcf_file_prefix(path):
    vcf_file_name = os.path.split(path)[1]
    if re.search('[.]vcf[.]gz$', vcf_file_name):
        return re.search('(.+)[.]vcf[.]gz$', vcf_file_name).group(1)
    else:
        raise Exceptions.InappropriateVCFSuffix(path)


def quit_if_exists(path):
    '''report an error if the file exists'''
    if os.path.exists(path):
        print('File already exists. Delete if you want to proceed: %s' % (path))
        sys.exit()


def download_from_url(url, download_to_dir, save_to_file = None):
    '''download from an url'''

    remote_basename = os.path.basename(urllib.parse.urlparse(url).path)
    if remote_basename:
        # Download to a temp file. If a download succeeds, rename the temp file.
        # If a download fails, the function will throw an exception. The temp file will be removed.

        local_path = os.path.join(download_to_dir, remote_basename) if not save_to_file else save_to_file
        quit_if_exists(local_path)
        temp_download_path = os.path.join(download_to_dir, 'temp_' + remote_basename)
        with tempfile.TemporaryDirectory(dir = download_to_dir) as temp_dir:
            with urllib.request.urlopen(url) as response:
                with open(temp_download_path, 'wb') as out_file:
                    print('Downloading from \"%s\"\n\tto \"%s\"' %(url, local_path))
                    shutil.copyfileobj(response, out_file)
            os.rename(temp_download_path, local_path)
        return local_path
    else:
        raise Exceptions.InvalidURL(url)


def decompress_gz_file(path):
    path_to_decompressed = os.path.splitext(path)[0]
    with gzip.open(path, 'rb') as f_in:
        with open(path_to_decompressed, 'wb') as f_out:
            print('Decompressing %s' %(path))
            shutil.copyfileobj(f_in, f_out)
    return path_to_decompressed


def download_grch38_ref_fasta_and_index(download_to_dir, save_to_file=None):
    '''download the human reference genome sequence GRChh38/hg38 from the NIH FTP site'''

    # download ftp web address (dated Feb 2021)
    url_grch38_fasta = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz'
    url_grch38_fasta_index = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai'
    
    # download and prepare files for vcf normalization using the bcftools
    path_to_ref_seq = download_from_url(url_grch38_fasta, download_to_dir)
    path_to_ref_seq = decompress_gz_file(path_to_ref_seq)
    download_from_url(url_grch38_fasta_index, download_to_dir)
    return path_to_ref_seq


def tabix_index_vcf(tabix_executable_path, vcf_path):
    '''
    index the input vcf using tabix, and the output index file will be written to the working directory

    tabix commands are exclusively "tabix -p vcf <input_vcf>", which generates an index file (.tbi) for an input file (<input_file>) whose file type is specified by "-p vcf".
    .tbi will be output to the current working directory by default.
    '''
    
    try:
        subprocess.run([tabix_executable_path, '-p', 'vcf', vcf_path], cwd = os.path.split(vcf_path)[0])
    except:
        import sys
        print('Error: cannot index the file: %s' %vcf_path)
        sys.exit(1)


def running_bcftools(list_bcftools_command, show_msg = None):
    '''
    run the bcftools following the commands stored in the list_bcftools_command

    "bcftools <common_options> <input_vcf>".
    "-Oz" (capitalized letter O) specifies the output type as compressed VCF (z). "-o" writes to a file rather than to default standard output.
    "--no-version" will cease appending version and command line information to the output VCF header.
    "-s sample_ID(s)" comma-separated list of samples to include or exclude if prefixed with "^".
    '''

    #print("%s [ %s ]" % (show_msg, ' '.join(list_bcftools_command))) if show_msg else print("Running [ %s ]" % (' '.join(list_bcftools_command)))
    print("%s" % (show_msg)) if show_msg else print("Running [ %s ]" % (' '.join(list_bcftools_command)))
    subprocess.run(list_bcftools_command)

def remove_vcf_and_index(path_to_vcf):
    '''remove the compressed vcf as well as the index file'''

    try:
        os.remove(path_to_vcf)
        os.remove(path_to_vcf + '.tbi')
        print("Removed intermediate files:\n\t%s\n\t%s" %(path_to_vcf, path_to_vcf + '.tbi'))
    except OSError as error_remove_tmp:
        print("Error: %s : %s" % (path_to_vcf, error_remove_tmp.strerror))

def get_vcf_pos_min_max(positions, flanking_bp = 100):
    ''' given input positions, return "<min_pos>-<max_pos>"  '''
    return '-'.join([str(min(positions)-flanking_bp), str(max(positions)+flanking_bp)])

def extract_pharmcat_pgx_regions(tabix_executable_path, input_vcf, output_dir, input_ref_pgx_vcf):
    '''
    extract pgx regions in input_ref_pgx_vcf from input_vcf and save variants to path_output
    '''

    print('Modify chromosome names.\nExtract PGx regions based on the input reference PGx position file.')
    path_output = os.path.join(output_dir, obtain_vcf_file_prefix(input_vcf) + '.pgx_regions.vcf.gz')

    input_vcf_cyvcf2 = VCF(input_vcf)
    input_ref_pgx_pos_cyvcf2 = VCF(input_ref_pgx_vcf)

    # get pgx regions in each chromosome
    input_ref_pgx_pos_pandas = allel.vcf_to_dataframe(input_ref_pgx_vcf)
    input_ref_pgx_pos_pandas['CHROM'] = input_ref_pgx_pos_pandas['CHROM'].replace({'chr':''}, regex=True).astype(str).astype(int)
    ref_pgx_regions = input_ref_pgx_pos_pandas.groupby(['CHROM'])['POS'].agg(get_vcf_pos_min_max).reset_index()
    # fix chr names
    chr_name_match = re.compile("^chr")
    if any(chr_name_match.match(line) for line in input_vcf_cyvcf2.seqnames):
        # chromosomes have leading 'chr' characters in the original VCF
        # pgx regions to be extracted
        ref_pgx_regions = ref_pgx_regions.apply(lambda row: ':'.join(row.values.astype(str)), axis=1).replace({'^':'chr'}, regex=True)
    else:
        # chromosomes do not have leading 'chr' characters in the original VCF
        # add chromosome name with leading 'chr' to the VCF header
        for single_chr in input_vcf_cyvcf2.seqnames:
            input_vcf_cyvcf2.add_to_header('##contig=<ID=chr' + single_chr + '>')
        # pgx regions to be extracted
        ref_pgx_regions = ref_pgx_regions.apply(lambda row: ':'.join(row.values.astype(str)), axis=1)

    # write to a VCF output file
    # header
    output_vcf_cyvcf2 = Writer(path_output, input_vcf_cyvcf2, mode="wz")
    # content
    for single_region in ref_pgx_regions:
        for single_variant in input_vcf_cyvcf2(single_region):
            single_variant.CHROM = re.sub(r'^([0-9]+)', r'chr\1', single_variant.CHROM)
            output_vcf_cyvcf2.write_record(single_variant)

    # close pipe
    input_vcf_cyvcf2.close()
    input_ref_pgx_pos_cyvcf2.close()
    output_vcf_cyvcf2.close()

    tabix_index_vcf(tabix_executable_path, path_output)

    return path_output


def merge_vcfs(bcftools_executable_path, tabix_executable_path, input_vcf, input_ref_pgx_vcf):
    '''
    merge the input VCF with a reference VCF of PGx core allele defining positions to enforce the same variant representation format across files

    bcftools merge <options> <input_vcf>". For bcftools common options, see running_bcftools().
    "-m both" allows both SNP and indel records can be multiallelic. But SNP and indel will not be merged as one multiallelic record.
    '''

    path_output = os.path.splitext(os.path.splitext(input_vcf)[0])[0] + '.pgx_merged.vcf.gz'

    bcftools_command_to_merge = [bcftools_executable_path, 'merge', '--no-version', '-m', 'both', '-Oz', '-o', path_output, input_vcf, input_ref_pgx_vcf]
    running_bcftools(bcftools_command_to_merge, show_msg = 'Enforcing the same variant representation as that in the reference PGx variant file')

    tabix_index_vcf(tabix_executable_path, path_output)

    return path_output

def normalize_vcf(bcftools_executable_path, tabix_executable_path, input_vcf, path_to_ref_seq, input_ref_pgx_vcf):
    '''
    normalize the input VCF against the human reference genome sequence GRCh38/hg38

    "bcftools norm <options> <input_vcf>". For bcftools common options, see running_bcftools().
    "-m+" joins biallelic sites into multiallelic records (+).
    "-f <ref_seq_fasta>" reference sequence. Supplying this option turns on left-alignment and normalization.
    "-c ws" when incorrect or missing REF allele is encountered, warn (w) and set/fix(s) bad sites. 's' will swap alleles and update GT and AC acounts. Importantly, s will NOT fix strand issues in a VCF.
    '''

    temp_normalized_vcf = os.path.splitext(os.path.splitext(input_vcf)[0])[0] + '.tmp_normed.vcf.gz'
    # normalize
    bcftools_command_to_normalize_vcf = [bcftools_executable_path, 'norm', '--no-version', '-m+', '-c', 'ws',  '-Oz', '-o', temp_normalized_vcf, '-f', path_to_ref_seq, input_vcf]
    running_bcftools(bcftools_command_to_normalize_vcf, show_msg = 'Normalize VCF')
    tabix_index_vcf(tabix_executable_path, temp_normalized_vcf)

    # extract only PGx positions
    path_output = os.path.splitext(os.path.splitext(input_vcf)[0])[0] + '.normalized.vcf.gz'
    bcftools_command_to_extract_only_pgx = [bcftools_executable_path, 'view', '--no-version', '-U', '-Oz', '-o', path_output, '-R', input_ref_pgx_vcf, temp_normalized_vcf]
    running_bcftools(bcftools_command_to_extract_only_pgx, show_msg = 'Retain only PGx positions in the normalized VCF') # run bcftools to merge VCF files
    tabix_index_vcf(tabix_executable_path, path_output)

    remove_vcf_and_index(temp_normalized_vcf)

    return path_output

def output_pharmcat_ready_vcf(input_vcf, output_dir, output_prefix):
    '''
    iteratively write to a PharmCAT-ready VCF for each sample

    "bcftools view <options> <input_vcf>". For bcftools common options, see running_bcftools().
    "-U" exclude sites without a called genotype, i.e., GT = './.'
    '''

    input_vcf_cyvcf2 = VCF(input_vcf)
    input_vcf_sample_list = input_vcf_cyvcf2.samples
    input_vcf_sample_list.remove('PharmCAT')

    # output each single sample to a separete VCF
    for single_sample in input_vcf_sample_list:
        print('Generating a PharmCAT-ready VCF for ' + single_sample)
        input_vcf_cyvcf2.set_samples(single_sample)

        # write to a VCF output file
        output_file_name = os.path.join(output_dir, output_prefix + '.' + single_sample + '.vcf')
        # header
        output_vcf_cyvcf2 = Writer(output_file_name, input_vcf_cyvcf2, mode = 'w')
        # content
        for single_var in input_vcf_cyvcf2:
            output_vcf_cyvcf2.write_record(single_var)
        output_vcf_cyvcf2.close()

    input_vcf_cyvcf2.close()


def output_missing_pgx_positions(bcftools_executable_path, input_vcf, input_ref_pgx_vcf, output_dir, output_prefix):
    '''
    generate a report VCF of missing PGx positions from the input VCF against a reference PGx VCF

    "bcftools isec <options> <input_vcf>". For bcftools common options, see running_bcftools().
    "-c" controls how to treat records with duplicate positions. "-c indels" means that all indel records are compatible, regardless of whether the REF and ALT alleles match or not. For duplicate positions, only the first indel record will be considered and appear on output.
    "-C" outputs positions present only in the first file but missing in the others.
    "-w" lists input files to output given as 1-based indices. "-w1" extracts and writes records only present in the first file (the reference PGx positions).
    '''

    output_file_name = os.path.join(output_dir, output_prefix + '.missing_pgx_var.vcf.gz')
    bcftools_command_to_report_missing_pgx = [bcftools_executable_path, 'isec', '--no-version', '-c', 'indels', '-w1', '-Oz', '-o', output_file_name, '-C', input_ref_pgx_vcf, input_vcf]
    running_bcftools(bcftools_command_to_report_missing_pgx, show_msg = 'Generating a report of missing PGx core allele defining positions')