#!/usr/bin/env python3

__author__ = 'BinglanLi'

import sys
from glob import glob
from pathlib import Path
from timeit import default_timer as timer

import numpy as np
import pandas as pd

import utilities as util


if __name__ == "__main__":
    import argparse

    # describe the tool
    parser = argparse.ArgumentParser(description='Identify overlapping diplotypes in unphased genetic data.')

    # input arguments
    parser.add_argument("-i", "--input",
                        type=str,
                        required=False,
                        help="File name or the directory to the PharmCAT JSONs of allele definitions.")
    # input arguments
    parser.add_argument("-m", "--missing",
                        type=bool,
                        required=False,
                        default=False,
                        help="To evaluate diplotypes with missing positions, True or False.")

    # output args
    parser.add_argument("-o", "--output-dir",
                        type=str,
                        required=False,
                        metavar='<dir>',
                        help="(Optional) directory for outputs.  Defaults to the current working directory.")
    parser.add_argument("-bf", "--base-filename",
                        type=str,
                        required=False,
                        default='predicted_pharmcat_calls',
                        help="(Optional) base filename of the output file.")

    parser.add_argument('-c', '--clinical-calls',
                        action='store_true',
                        required=False,
                        help='Generate data for clinical calls')

    # parse arguments
    args = parser.parse_args()

    try:
        # define output file name
        m_output_basename: str = args.base_filename
        # define output file path
        m_output_dir: Path = Path(args.output_dir) if args.output_dir else Path.cwd()
        if not args.output_dir:
            print(f'Output directory is not specified.\n'
                  f'\tUse the current working directory: {m_output_dir}')

        # define input
        m_input: Path
        if args.input:
            m_input = Path(args.input).absolute()
        else:
            script_dir: Path = Path(globals().get("__file__", "./_")).absolute().parent
            is_repo: bool = script_dir.absolute().as_posix().endswith("scripts/diplotype_comparison")
            if is_repo:
                m_input = script_dir.parent.parent / 'main/resources/org/pharmgkb/pharmcat/definition/alleles'
                print(f'Running in repo.')
            else:
                m_input = Path.cwd()
                print(f'No input provided. '
                      f'Looking for the allele definition JSON files under the current working directory: {m_input}')
        allele_definition_jsons: list = []
        if m_input.is_file():
            allele_definition_jsons = [m_input]
        elif m_input.is_dir():
            print(f'Looking for the allele definition JSON files under {m_input}')
            # get the list allele definition files
            allele_definition_jsons = glob(str(m_input.joinpath('*_translation.json')))
        # check whether there are any json files
        if len(allele_definition_jsons) == 0:
            print(f'No allele definition JSON files are found under {m_input}')
            sys.exit(1)

        if args.missing:
            m_output_file: Path = m_output_dir / (m_output_basename + '_missing.tsv')
        else:
            m_output_file: Path = m_output_dir / (m_output_basename + '.tsv')
        # delete file
        m_output_file.unlink(missing_ok=True)
        print()
        if args.clinical_calls:
            print(f'Generating clinical call data')
            if args.base_filename == 'predicted_pharmcat_calls':
                m_output_file: Path = m_output_dir / 'clinical_calls_diff.tsv'
        print(f'Saving to {m_output_file}\n')

        # process each gene
        for json_file in allele_definition_jsons:
            # read the json file
            start = timer()
            filename: str = json_file.split('/')[-1]
            # get the gene name
            gene: str = filename.split('_')[0]

            if args.clinical_calls and (gene == 'CYP2D6' or gene == 'DPYD'):
                continue

            # read files
            print(f'Processing {filename}')
            json_data: dict = util.read_json(json_file)

            # get necessary information of the allele-defining positions, including hgvs names, positions, references
            allele_defining_variants: dict = util.get_allele_defining_variants(json_data)

            # read in the list of alleles, allele-defining positions, and defining genotypes
            allele_definitions: dict = util.get_allele_definitions(json_data, allele_defining_variants)

            # get an allele-position array to identify alleles that do not share definitions with others
            allele_arrays: np.ndarray = util.convert_dictionary_to_numpy_array(allele_definitions)
            allele_position_array: np.ndarray = np.sum(allele_arrays, axis=1).astype(bool).astype('int8')

            # identify reference alleles
            allele_names: np.ndarray = np.array([*allele_definitions])
            bool_reference_alleles: np.ndarray = np.all(allele_position_array, axis=1)
            reference_alleles: list[str] = list(allele_names[bool_reference_alleles])

            # identify the 'busy' positions that define only one allele besides the reference allele
            idx_busy_positions: np.ndarray = np.where(np.sum(allele_position_array, axis=0) > 2)[0]

            # identify alleles that share definitions with others besides the reference allele
            idx_gregarious_alleles: np.ndarray = np.unique(np.where(allele_position_array[:, idx_busy_positions])[0])
            gregarious_alleles: list[str] = list(allele_names[idx_gregarious_alleles])

            # skip the gene if none of its alleles shares any allele-defining positions with others
            if len(gregarious_alleles) == 0:
                print(f'\tSkipping - no alleles share definitions with others')
                continue

            # find all possible outcomes of alternative calls for each diplotype
            dict_predicted_calls = dict()
            if args.missing:
                if gene == 'CYP2D6':
                    print(f'\tSkipping - too complex')
                    continue

                # identify the combinations of missing positions to evaluate
                hgvs_names = [*allele_defining_variants]
                missing_combinations = util.find_missingness_combinations(allele_position_array, list(allele_names),
                                                                          hgvs_names, reference_alleles,
                                                                          gregarious_alleles)

                # if the gene does not have positions whose absence will cause ambiguous pharmcat calls, then skip
                if not missing_combinations:
                    print(f'\tSkipping - no ambiguous calls')
                    continue
                else:
                    print(f'\tNumber of combinations = {len(missing_combinations)}')
                    for m in missing_combinations:
                        # find possible calls for each missing position
                        dict_predicted_calls_m = util.predict_pharmcat_calls(allele_defining_variants,
                                                                             allele_definitions, m)

                        # append to dict_predicted_calls
                        for key in dict_predicted_calls_m:
                            dict_predicted_calls[key] = dict_predicted_calls.get(key, []) + dict_predicted_calls_m[key]
            else:
                dict_predicted_calls = util.predict_pharmcat_calls(allele_defining_variants,
                                                                   allele_definitions,
                                                                   gregarious_alleles)

            # convert the python dictionary to a pandas data frame for output
            rez = pd.DataFrame.from_dict(dict_predicted_calls)
            # add gene name to the data frame for readability in the output
            rez['gene'] = gene
            # write to an output
            if len(rez):
                if args.clinical_calls:
                    # keep only diplotypes that have multiple top-scored calls
                    mask = rez[rez['actual'].str.contains(';')]
                    # remove duplicates
                    rows = mask.drop_duplicates(subset=['gene', 'actual'])
                    cols = ['gene', 'actual']
                else:
                    rows = rez
                    cols = rows.columns
                if m_output_file.is_file():
                    mode = 'a'
                    header = False
                else:
                    mode = 'w'
                    header = True
                rows.to_csv(m_output_file.absolute(), mode=mode, sep="\t", columns=cols, header=header, index=False)
            else:
                print(f'\tNo alternative calls')

            # show processing time
            end = timer()
            print(f'\tFinished in %.2f seconds' % (end - start))
    except Exception as e:
        print(e)
        sys.exit(1)
