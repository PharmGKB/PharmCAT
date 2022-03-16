#! /usr/bin/env python
__author__ = 'BinglanLi'

import allel
import gzip
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.parse
import urllib.request
import copy

import vcf_preprocess_exceptions as Exceptions


# chromosome names
_chr_invalid = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
                "18", "19", "20", "21", "22", "X", "Y", "M", "MT", "chrMT"]
_chr_valid = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
              "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
              "chr22", "chrX", "chrY", "chrM", "chrM", "chrM"]
_chr_valid_sorter = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
                     "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
                     "chr22", "chrX", "chrY", "chrM"]
# check if two chr arrays are of the same length
if len(_chr_invalid) != len(_chr_valid):
    print("Error in internal chromosome mapping arrays")
    sys.exit(1)


def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


def bgzip_file(bgzip_path, vcf_path):
    """
    bgzip the file
    """

    try:
        print("Bgzipping", vcf_path)
        subprocess.run([bgzip_path, '-f', vcf_path], check=True)
    except Exception as e:
        print('Failed to bgzip %s' % vcf_path)
        # comment out this traceback function as subprocess(check = true) should report detailed errors
        # traceback.print_exception(type(e), e, e.__traceback__)
        sys.exit(1)


def bgzipped_vcf(bgzip_path, file):
    """
    make sure file is bgzipped
    """
    if not is_gz_file(file):
        print("Bgzipping VCF!")
        bgzip_file(bgzip_path, file)
        file = file + '.gz'
        if os.path.exists(file + '.tbi'):
            print("Removing pre-existing .tbi")
            os.remove(file + '.tbi')
    return file


def byte_decoder(a):
    """ Decode byte data into utf-8 """
    return a.decode("utf-8")


def download_from_url(url, download_to_dir, save_to_file=None, force_update=False):
    """download from an url"""

    remote_basename = os.path.basename(urllib.parse.urlparse(url).path)
    if remote_basename:
        local_path = os.path.join(download_to_dir, remote_basename) if not save_to_file else save_to_file
        if os.path.exists(local_path):
            if force_update:
                os.remove(local_path)
            else:
                return local_path

        # Download to a temp file. If a download succeeds, rename the temp file.
        # If a download fails, the function will throw an exception. The temp file will be removed.
        with tempfile.TemporaryDirectory(dir=download_to_dir) as temp_dir:
            temp_download_path = os.path.join(temp_dir, 'temp_' + remote_basename)
            with urllib.request.urlopen(url) as response:
                with open(temp_download_path, 'wb') as out_file:
                    print('Downloading from \"%s\"\n\tto \"%s\"' % (url, local_path))
                    shutil.copyfileobj(response, out_file)
            os.rename(temp_download_path, local_path)
        return local_path
    else:
        raise Exceptions.InvalidURL(url)


def get_default_grch38_ref_fasta_and_index(download_to_dir, force_update=False):
    """download the human reference genome sequence GRCh38/hg38"""

    ref_file = os.path.join(download_to_dir, 'reference.fasta.bgz')
    if os.path.exists(ref_file) and not force_update:
        return ref_file

    tar_file = download_from_url('https://zenodo.org/record/5572839/files/GRCh38_reference_fasta.tar?download=1',
                                 download_to_dir, None, force_update)
    with tarfile.open(tar_file, 'r') as tar:
        tar.extractall(path=download_to_dir)
    os.remove(tar_file)

    return ref_file


def tabix_index_vcf(tabix_path, vcf_path):
    """
    index the input vcf using tabix, and the output index file will be written to the working directory

    tabix commands are exclusively "tabix -p vcf <input_vcf>", which generates an index file (.tbi)
        for an input file (<input_file>) whose file type is specified by "-p vcf".
    .tbi will be output to the current working directory by default.
    """

    try:
        print('Generating index (' + vcf_path + '.tbi)')
        subprocess.run([tabix_path, '-p', 'vcf', vcf_path], check=True)
    except Exception as e:
        print('Failed to index %s' % vcf_path)
        # comment out this traceback function as subprocess(check = true) should report detailed errors
        # traceback.print_exception(type(e), e, e.__traceback__)
        sys.exit(1)


def run_bcftools(list_bcftools_command, show_msg=None):
    """
    run the bcftools following the commands stored in the list_bcftools_command

    "bcftools <common_options> <input_vcf>".
    "-Oz" (capitalized letter O) specifies the output type as compressed VCF (z). "-o" writes to a file rather than to
    default standard output.
    "--no-version" will cease appending version and command line information to the output VCF header.
    "-s sample_ID(s)" comma-separated list of samples to include or exclude if prefixed with "^".
    "-r chr|chr:pos|chr:beg-end|chr:beg-[,…]" extracts comma-separated list of regions
    """

    print("%s" % show_msg) if show_msg else print("Running [ %s ]" % (' '.join(list_bcftools_command)))

    p = subprocess.run(list_bcftools_command, stderr=subprocess.PIPE)
    if p.returncode != 0:
        print(p.stderr.decode("utf-8"))
        sys.exit(1)


def obtain_vcf_sample_list(bcftools_path, path_to_vcf):
    """
    obtain a list of samples from the input VCF

    "bcftools query <options> <input_vcf>". For bcftools common options, see running_bcftools().
    "-l" list sample names and exit.
    Samples are delimited by '\\\\n' and the last line ends as 'last_sample_ID\\\\n\\\\n'.
    """

    output = subprocess.check_output([bcftools_path, 'query', '-l', path_to_vcf], universal_newlines=True)
    vcf_sample_list = output.split('\n')[:-1]  # remove the black line at the end
    return vcf_sample_list


def remove_vcf_and_index(path_to_vcf):
    """remove the compressed vcf as well as the index file"""

    try:
        os.remove(path_to_vcf)
        os.remove(path_to_vcf + '.tbi')
        print("Removed intermediate files:\n\t%s\n\t%s" % (path_to_vcf, path_to_vcf + '.tbi'))
    except OSError as error_remove_tmp:
        print("Error: %s : %s" % (path_to_vcf, error_remove_tmp.strerror))


def _get_vcf_pos_min_max(positions, flanking_bp=100):
    """ given input positions, return "<min_pos>-<max_pos>"  """
    return '-'.join([str(min(positions) - flanking_bp), str(max(positions) + flanking_bp)])


def _is_valid_chr(input_vcf):
    chr_status = 'none'
    with gzip.open(input_vcf) as f:
        for line in f:
            try:
                line = byte_decoder(line)
            except:
                line = line
            if line[0] != '#':
                fields = line.rstrip().split()
                if fields[0] in _chr_valid:
                    chr_status = 'true'
                    break
                elif fields[0] in _chr_invalid:
                    chr_status = 'false'
                    break
                else:
                    break
    return chr_status


def extract_regions_from_single_file(bcftools_path, tabix_path, input_vcf, pgx_vcf, output_dir, output_prefix,
                                     sample_list):
    """
    Rename chromosomes in input vcf according to a chr-renaming mapping file.
    Extract pgx regions from input_vcf based on the ref_pgx.

    "bcftools annotate <options> <input_vcf>". For bcftools common options, see running_bcftools().
    "--rename-chrs" renames chromosomes according to the map in file_rename_chrs.
    """

    # output path
    path_output = os.path.join(output_dir, 'PharmCAT_preprocess_' + output_prefix + '.pgx_regions.vcf.gz')

    print("")
    print("Processing", input_vcf)
    # create index if not already existed
    if not os.path.exists(input_vcf + '.tbi'):
        tabix_index_vcf(tabix_path, input_vcf)

    # check whether input if a gVCF
    with gzip.open(input_vcf, 'r') as in_f:
        for line in in_f:
            try:
                line = byte_decoder(line)
            except:
                line = line
            if line[0:2] == '##':
                if ('ALT' in line) and ('ID=NON_REF' in line):
                    print('=============================================================\n'
                          'Proprocessor currently does not support gVCF.\n'
                          '=============================================================\n')
                    sys.exit(1)
            else:
                break

    # obtain PGx regions to be extracted
    df_ref_pgx_pos = allel.vcf_to_dataframe(pgx_vcf)
    df_ref_pgx_pos['CHROM'] = df_ref_pgx_pos['CHROM'].astype("category")
    df_ref_pgx_pos['CHROM'].cat.set_categories(_chr_valid_sorter, inplace=True)
    ref_pgx_regions = df_ref_pgx_pos.groupby(['CHROM'])['POS'].agg(_get_vcf_pos_min_max).reset_index()
    ref_pgx_regions.dropna(axis=0, subset=['POS'], how='any', inplace=True)
    # add a special case for 'chrMT'
    idx_chrM = ref_pgx_regions.index[ref_pgx_regions['CHROM'] == 'chrM']
    ref_pgx_regions = ref_pgx_regions.append(
        ref_pgx_regions.loc[idx_chrM].assign(**{'CHROM': 'chrMT'}), ignore_index=True)

    # generate a temp dir to extract pgx regions and, if necessary, rename chr
    with tempfile.TemporaryDirectory(suffix='extract_pgx_regions', dir=output_dir) as temp_dir:
        # generate temp file of sample list
        file_sample_list = os.path.join(temp_dir, 'sample_list.txt')
        with open(file_sample_list, 'w+') as f:
            for single_sample in sample_list:
                f.write("%s\n" % single_sample)

        # create a temporary chromosome mapping file
        file_chr_rename = os.path.join(temp_dir, 'rename_chr.txt')
        with open(file_chr_rename, 'w+') as f:
            for i in range(len(_chr_invalid)):
                f.write(_chr_invalid[i] + "\t" + _chr_valid[i] + "\n")

        # validate chromosome formats
        if _is_valid_chr(input_vcf) == 'true':
            # format the pgx regions to be extracted
            ref_pgx_regions = ",".join(ref_pgx_regions.apply(lambda row: ':'.join(row.values.astype(str)), axis=1))
        elif _is_valid_chr(input_vcf) == 'false':
            # format pgx regions
            ref_pgx_regions = ",".join(
                ref_pgx_regions.apply(lambda row: ':'.join(row.values.astype(str)), axis=1).replace({'chr': ''},
                                                                                                    regex=True))
        else:
            print("The CHROM column does not conform with either 'chr##' or '##' format.")
            sys.exit(1)

        # extract pgx regions and modify chromosome names if necessary
        bcftools_command = [bcftools_path, 'annotate', '--no-version', '-S', file_sample_list,
                            '--rename-chrs', file_chr_rename, '-r', ref_pgx_regions, '-i', 'ALT="."', '-k',
                            '-Oz', '-o', path_output, input_vcf]
        run_bcftools(bcftools_command,
                     show_msg='Extracting PGx regions and modifying chromosome names for %s.' % input_vcf)

    # index the output PGx file
    tabix_index_vcf(tabix_path, path_output)

    return path_output


def extract_regions_from_multiple_files(bcftools_path, tabix_path, bgzip_path, input_list, ref_pgx,
        output_dir, output_prefix, sample_list):
    """
    iterate through the list of input files
    """

    path_output = os.path.join(output_dir, 'PharmCAT_preprocess_' + output_prefix + '.pgx_regions.vcf.gz')

    with tempfile.TemporaryDirectory(suffix='concat_input_vcfs', dir=output_dir) as temp_dir:
        # process the vcfs in the input list one by one
        preprocessed_file_list = []
        i = 1
        with open(input_list, 'r') as file:
            for line in file:
                line = line.strip()
                if os.path.isfile(line):
                    print("")
                    print("Processing", line)
                    if not os.path.exists(line):
                        print("Cannot find", line)
                        continue
                    line = bgzipped_vcf(bgzip_path, line)

                    temp_output_prefix = output_prefix + '_' + str(i)
                    single_file = extract_regions_from_single_file(bcftools_path, tabix_path, line, ref_pgx,
                                                                   temp_dir, temp_output_prefix, sample_list)
                    preprocessed_file_list.append(single_file)
                    i += 1
                else:
                    print("Warning: Skip %s because the file does not exist." % line)

        # generate a temporary list of files to be concatenated
        temp_file_list = os.path.join(temp_dir, "temp_file_list.txt")
        with open(temp_file_list, 'w+') as f:
            for j in range(len(preprocessed_file_list)):
                f.write(preprocessed_file_list[j] + "\n")

        # concatenate vcfs
        bcftools_command = [bcftools_path, 'concat', '--no-version', '-a', '-f', temp_file_list, '-Oz', '-o',
                            path_output]
        run_bcftools(bcftools_command, show_msg='Concatenating chromosome VCFs.')

    # index the concatenated VCF
    tabix_index_vcf(tabix_path, path_output)
    return path_output


def normalize_vcf(bcftools_path, tabix_path, input_vcf, ref_seq):
    """
    Normalize the input VCF against the human reference genome sequence GRCh38/hg38

    "bcftools norm <options> <input_vcf>". For bcftools common options, see running_bcftools().
    "-m +|-" joins biallelic sites into multiallelic records (+)
        and convert multiallelic records into uniallelic format (-).
    "-f <ref_seq_fasta>" reference sequence. Supplying this option turns on left-alignment and normalization.
    "-c ws" when incorrect or missing REF allele is encountered, warn (w) and set/fix(s) bad sites.  's' will swap
    alleles and update GT and AC acounts. Importantly, s will NOT fix strand issues in a VCF.
    """

    path_output = os.path.splitext(os.path.splitext(input_vcf)[0])[0] + '.normalized.vcf.gz'

    bcftools_command = [bcftools_path, 'norm', '--no-version', '-m-', '-c', 'ws', '-Oz', '-o',
                        path_output, '-f', ref_seq, input_vcf]
    run_bcftools(bcftools_command, show_msg='Normalizing VCF')
    tabix_index_vcf(tabix_path, path_output)

    return path_output


def filter_pgx_variants(bcftools_path, tabix_path, bgzip_path, input_vcf, ref_seq, ref_pgx,
                        missing_to_ref, output_dir, output_prefix):
    """
    Extract specific pgx positions that are present in the reference PGx VCF
    Generate a report of PGx positions that are missing in the input VCF

    "bcftools isec <options> <input_vcf>". For bcftools common options, see running_bcftools().
    "-c none" only records with the same CHR, POS, REF and ALT are considered identical
    "-C" outputs positions present only in the first file but missing in the others.
    "-n=2" positions shared by both inputs
    "-i 'TYPE=snp' " to include only SNPs
    "-w" lists input files to output given as 1-based indices.
        "-w1" extracts and writes records only present in the first file (the reference PGx positions).
    """

    path_output = os.path.splitext(os.path.splitext(input_vcf)[0])[0] + '.multiallelic.vcf.gz'

    with tempfile.TemporaryDirectory(suffix='extract_pgx_variants', dir=output_dir) as temp_dir:
        # convert reference PGx variants to the uniallelic format
        # needed for extracting exact PGx positions and generating an accurate missing report
        ref_pgx_uniallelic = os.path.join(temp_dir, 'temp_ref_pgx_uniallelic.vcf.gz')
        bcftools_command = [bcftools_path, 'norm', '--no-version', '-m-', '-c', 'ws', '-f', ref_seq,
                            '-Oz', '-o', ref_pgx_uniallelic, ref_pgx]
        run_bcftools(bcftools_command, show_msg='Preparing the reference PGx VCF')
        tabix_index_vcf(tabix_path, ref_pgx_uniallelic)

        '''
        extracg PGx positions from input using "bcftools view"
        '''
        # extract any variants with matching positions
        input_pgx_only = os.path.join(temp_dir, 'temp_input_pgx_variants_only.vcf.gz')
        bcftools_command = [bcftools_path, 'view', '--no-version', '-T', ref_pgx_uniallelic,
                            '-Oz', '-o', input_pgx_only, input_vcf]
        run_bcftools(bcftools_command, show_msg='Retaining PGx positions, regardless of alleles')
        tabix_index_vcf(tabix_path, input_pgx_only)

        '''
        extract PGx positions from input VCF
        add in missing multiallelic variants as '0|0'
        bcftools will take care of the rest (sorting, normalization, etc.)
        '''
        # create a dictionary of PharmCAT reference PGx positions
        # the ref_pos_static (dict) will be used to static dictionary of reference PGx positions
        # the ref_pos_dynamic (dict) will be used to remain only PGx pos in the input
        ref_pos_dynamic = {}
        with gzip.open(ref_pgx_uniallelic, 'r') as in_f:
            for line in in_f:
                try:
                    line = byte_decoder(line)
                except:
                    line = line
                # skip headers
                if line[0] == '#':
                    continue
                # read file
                line = line.rstrip('\n')
                fields = line.split('\t')
                # ref_pos_dynamic: a nested dictionary
                # ref_pos_dynamic[(chr,pos)][(ref, alt)] = all fields except GT
                # initiate the dictionary if the first key pair doesn't exist yet
                if (fields[0], fields[1]) not in ref_pos_dynamic.keys():
                    ref_pos_dynamic[(fields[0], fields[1])] = {}
                ref_pos_dynamic[(fields[0], fields[1])][(fields[3], fields[4])] = fields[0:9]
        ref_pos_static = copy.deepcopy(ref_pos_dynamic)

        # use input_pos to record PGx positions that are present in the input VCF
        input_pos = []
        # this list saves genetic variants concurrent at PGx positions
        non_pgx_records = []
        vcf_pgx_only = os.path.join(temp_dir, 'temp_update_pgx_variants_annotations.vcf')
        with open(vcf_pgx_only, 'w') as out_f:
            # get header of samples from merged vcf, add in new contig info
            print('Updating VCF header and PGx position annotations')
            with gzip.open(input_pgx_only) as in_f:
                for line in in_f:
                    try:
                        line = byte_decoder(line)
                    except:
                        line = line
                    # process header lines, skip contig
                    if line[0:8] == '##contig':
                        continue  # skip contig info in the original vcf
                    # process header lines, except lines about contig
                    elif line[0:2] == '##':
                        out_f.write(line)
                    # process the sample line, add PGx-related lines
                    elif line[0:6] == '#CHROM':
                        for single_chr in _chr_valid_sorter:
                            out_f.write('##contig=<ID=' + single_chr + ',assembly=GRCh38.p13,species="Homo '
                                                                       'sapiens">\n')
                        out_f.write('##INFO=<ID=PX,Number=.,Type=String,Description="Gene">\n')
                        out_f.write('##INFO=<ID=POI,Number=0,Type=Flag,Description="Position of Interest but not'
                                    ' part of an allele definition">\n')
                        out_f.write('##FILTER=<ID=PCATxREF,Description="Reference allele does not match PharmCAT '
                                    'reference alleles">\n')
                        out_f.write('##FILTER=<ID=PCATxALT,Description="Alternate alleles do not match PharmCAT '
                                    'alternate alleles">\n')
                        out_f.write(line)
                        # get the number of samples
                        line = line.rstrip('\n')
                        fields = line.split('\t')
                        n_sample = len(fields) - 9
                    # scan the genotype data
                    else:
                        line = line.rstrip('\n')
                        fields = line.split('\t')
                        input_chr_pos = (fields[0], fields[1])

                        # match chromosome positions
                        if input_chr_pos in ref_pos_static:
                            '''
                            match REF and ALT alleles
                            
                            1. positions with matching REF and ALT
                            1.1. update rsID and gene name in FORMAT; bcftools will normalize
                            2. homozygous reference SNPs with matching REF and ALT='.'
                            2.1. PGx SNP: update rsID and gene name in FORMAT; bcftools will normalize
                            '''
                            # list out REF alleles at a position
                            ref_alleles = [x[0] for x in ref_pos_static[input_chr_pos].keys()]
                            alt_alleles = [x[1] for x in ref_pos_static[input_chr_pos].keys()]
                            # check if the reference PGx variant is a SNP
                            len_ref = sum(map(len, ref_alleles))
                            len_alt = sum(map(len, alt_alleles))
                            is_snp = ((len_ref - len_alt) == 0)

                            input_ref_alt = (fields[3], fields[4])
                            # positions with matching REF and ALT
                            if input_ref_alt in ref_pos_static[input_chr_pos]:
                                # record input PGx positions in a list
                                # use this to fill up multiallelic ALT later
                                if input_chr_pos not in input_pos:
                                    input_pos.append(input_chr_pos)

                                # update ID col
                                if fields[2] == '.':
                                    fields[2] = ref_pos_static[input_chr_pos][input_ref_alt][2]
                                elif ref_pos_static[input_chr_pos][input_ref_alt][2] not in fields[2]:
                                    fields[2] = fields[2] + ';' + ref_pos_static[input_chr_pos][input_ref_alt][2]
                                # update INFO field
                                if fields[7] == '.':
                                    fields[7] = ref_pos_static[input_chr_pos][input_ref_alt][7]
                                else:
                                    fields[7] = ref_pos_static[input_chr_pos][input_ref_alt][7] + ';' + fields[7]
                                # concat and write to output file
                                line = '\t'.join(fields)
                                out_f.write(line + '\n')

                                # eliminattion: remove the dictionary item so that it won't be used again
                                if input_chr_pos in ref_pos_dynamic:
                                    ref_pos_dynamic[input_chr_pos].pop(input_ref_alt)
                                    # remove a position if all of its alts are present in the input
                                    if ref_pos_dynamic[input_chr_pos] == {}:
                                        del ref_pos_dynamic[input_chr_pos]
                            # homozygous reference SNPs with matching REF and ALT='.'
                            elif is_snp and (fields[3] in ref_alleles) and (fields[4] == '.'):
                                # record input PGx positions in a list
                                # use this to fill up multiallelic ALT later
                                if input_chr_pos not in input_pos:
                                    input_pos.append(input_chr_pos)

                                # update id when the id is not in the input
                                ref_id = list(set([x[2] for x in ref_pos_static[input_chr_pos].values()]))
                                if fields[2] == '.':
                                    fields[2] = ';'.join(ref_id)
                                else:
                                    input_id = fields[2].split(';')
                                    fields[2] = ';'.join(set(ref_id + input_id))

                                # update info
                                info_col = ';'.join(set([x[7] for x in ref_pos_static[input_chr_pos].values()]))
                                fields[7] = ';'.join([fields[7], info_col]) if fields[7] != '.' else info_col

                                alt_alleles = [x[1] for x in ref_pos_static[input_chr_pos].keys()]
                                for i in range(len(alt_alleles)):
                                    # update contents, concatenate cols from input and from the ref PGx VCF
                                    fields[3] = ref_alleles[i]
                                    fields[4] = alt_alleles[i]
                                    # concat and write to output file
                                    line = '\t'.join(fields)
                                    out_f.write(line + '\n')
                                    # eliminattion: remove the dictionary item so that it won't be used again
                                    if input_chr_pos in ref_pos_dynamic:
                                        ref_pos_dynamic[input_chr_pos].pop((ref_alleles[i], alt_alleles[i]))
                                        # remove a position if all of its alts are present in the input
                                        if ref_pos_dynamic[input_chr_pos] == {}:
                                            del ref_pos_dynamic[input_chr_pos]
                            # flag if the variant doesn't match PharmCAT ALT
                            elif fields[3] in ref_alleles:
                                for i in range(len(ref_alleles)):
                                    print('=============================================================\n\n'
                                          'Warning: \"%s:%s REF=%s ALT=%s\" does not match PharmCAT expectation '
                                          'of ALT at "%s:%s REF=%s ALT=%s"\n\n'
                                          '=============================================================\n'
                                          % (fields[0], fields[1], fields[3], fields[4],
                                             fields[0], fields[1], ref_alleles[i], alt_alleles[i]))
                                # update filter
                                fields[6] = 'PCATxALT'
                                # output the  line
                                line = '\t'.join(fields)
                                non_pgx_records.append(line)
                                # out_f.write(line + '\n')
                            # flag if the variant doesn't match PharmCAT REF nor ALT
                            else:
                                for i in range(len(ref_alleles)):
                                    print('=============================================================\n\n'
                                          'Warning: \"%s:%s REF=%s ALT=%s\" does not match PharmCAT expectation '
                                          'of REF at "%s:%s REF=%s ALT=%s"\n\n'
                                          '=============================================================\n'
                                          % (fields[0], fields[1], fields[3], fields[4],
                                             fields[0], fields[1], ref_alleles[i], alt_alleles[i]))
                                # update filter
                                fields[6] = 'PCATxREF'
                                # output the  line
                                line = '\t'.join(fields)
                                non_pgx_records.append(line)
                                # out_f.write(line + '\n')
            # output genetic variants that concurred at PGx positions
            for line in non_pgx_records:
                out_f.write(line + '\n')
            # if missing_to_ref is true, output all missing positions
            if missing_to_ref:
                for input_ref_alt in ref_pos_dynamic.values():
                    for val in input_ref_alt.values():
                        line = '\t'.join(val + ['0|0'] * n_sample)
                        out_f.write(line + '\n')
            # otherwise, only complete lines for multiallelic loci
            else:
                for single_pos in input_pos:
                    if single_pos in ref_pos_dynamic:
                        for val in ref_pos_dynamic[single_pos].values():
                            line = '\t'.join(val + ['0|0'] * n_sample)
                            out_f.write(line + '\n')

        # sort vcf
        sorted_vcf = os.path.join(temp_dir, 'temp_sorted.vcf.gz')
        bcftools_command = [bcftools_path, 'sort', '-Oz', '-o', sorted_vcf, vcf_pgx_only]
        run_bcftools(bcftools_command, show_msg='Sorting file')
        tabix_index_vcf(tabix_path, sorted_vcf)

        # This normalization enforces the output to comply with PharmCAT format
        bcftools_command = [bcftools_path, 'norm', '--no-version', '-m+', '-c', 'ws', '-f', ref_seq,
                            '-Oz', '-o', path_output, sorted_vcf]
        run_bcftools(bcftools_command,
                     show_msg='Enforcing the variant representation per PharmCAT')
        tabix_index_vcf(tabix_path, path_output)

        # report missing positions in the input VCF
        missing_report = os.path.join(output_dir, output_prefix + '.missing_pgx_var.vcf')
        with open(missing_report, 'w') as out_f:
            print('Generating a report of missing PGx allele defining positions')
            # get VCF header from the reference PGx VCF
            with gzip.open(ref_pgx, 'r') as in_f:
                for line in in_f:
                    try:
                        line = byte_decoder(line)
                    except:
                        line = line
                    # process header lines
                    if line[0:2] == '##':
                        out_f.write(line)
                    # process the sample line, add PGx-related lines
                    elif line[0:6] == '#CHROM':
                        line = line.rstrip('\n')
                        fields = line.split('\t')
                        fields[-1] = 'Missing'
                        out_f.write('\t'.join(fields) + '\n')
                    else:
                        break
            # output positions that were not detected in the input VCF
            for input_ref_alt in ref_pos_dynamic.values():
                for val in input_ref_alt.values():
                    line = '\t'.join(val + ['0|0'])
                    out_f.write(line + '\n')
        # bgzip the missing report
        bgzipped_vcf(bgzip_path, missing_report)
    return path_output


def output_pharmcat_ready_vcf(bcftools_path, input_vcf, output_dir, output_prefix, sample_list):
    """
    iteratively write to a PharmCAT-ready VCF for each sample

    "bcftools view <options> <input_vcf>". For bcftools common options, see running_bcftools().
    "-U" exclude sites without any called genotype, i.e., all GT = './.'
    "--force-samples" only warn about unknown subset samples
    """

    for single_sample in sample_list:
        output_file_name = os.path.join(output_dir, output_prefix + '.' + single_sample + '.vcf')
        bcftools_command = [bcftools_path, 'view', '--no-version', '--force-samples', '-U', '-s', single_sample,
                            '-Ov', '-o', output_file_name, input_vcf]
        run_bcftools(bcftools_command,
                     show_msg='Generating a PharmCAT-ready VCF for ' + single_sample)
