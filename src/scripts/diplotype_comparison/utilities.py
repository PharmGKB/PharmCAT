#! /usr/bin/env python
__author__ = 'BinglanLi'

import json
import sys
from pathlib import Path
from typing import List, Set, Optional

import numpy as np


# todo: convert sys.exit(1) to reportableExceptions

# define the list of wobble genotypes
_wobble_genotype_list: list = ['S', 'Y', 'M', 'K', 'R', 'W', 'V', 'H', 'D', 'B', 'N']
# build a dictionary of wobble genotypes and their corresponding base pairs
_wobble_match_table: dict[str, Set[str]] = {
    'M': {'A', 'C'},
    'R': {'A', 'G'},
    'W': {'A', 'T'},
    'S': {'C', 'G'},
    'Y': {'C', 'T'},
    'K': {'G', 'T'},
    'V': {'A', 'C', 'G'},
    'H': {'A', 'C', 'T'},
    'D': {'A', 'G', 'T'},
    'B': {'C', 'G', 'T'},
    'N': {'A', 'C', 'G', 'T'}
}


def read_json(file: Path) -> dict[str, str]:
    """
    Read json files
    :param file: path to the input json file
    :return: json_data: a dictionary of json data
    """
    try:
        # open and read json files
        with open(file, 'r') as json_file:
            json_data: dict = json.load(json_file)
        return json_data

    except FileNotFoundError:
        print(f"File '{file}' not found.")
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")


def get_allele_defining_variants(json_data: dict) -> dict[str, dict[str, str]]:
    """
    Extract information of allele-defining positions
    :param json_data: a dictionary of allele-defining variants and named alleles from PharmCAT JSON
    :return: a dictionary of {hgvs_name: {
        position: genomic_coordinate,
        ref: reference_genotype} }
    """
    # extract desirable information of allele-defining positions from the json data
    try:
        hgvs_names: List[str] = [entry['chromosomeHgvsName'] for entry in json_data['variants']]
        positions: List[str] = [entry['position'] for entry in json_data['variants']]
        reference_genotypes: List[str] = [entry['ref'] for entry in json_data['variants']]

        # check whether there are any empty values
        if None in hgvs_names:
            print('One of the HGVS names is empty.')
        if None in positions:
            print('Genomic coordinate is empty for one of the allele definition positions.')
        if None in reference_genotypes:
            print('Reference genotype is empty for one of the allele definition positions.')

        # generate a dictionary of all allele-defining positions and its reference genotypes
        allele_defining_variants: dict[str, dict[str, str]] = {
            h: {
                'position': p,
                'ref': r
            } for h, p, r in zip(hgvs_names, positions, reference_genotypes)}

        return allele_defining_variants

    except Exception as e:
        # print warnings if any of the attributes does not exist in the json file
        print(f'An error occurred in one of the {e}.'
              f'Check whether the attributes exist for all allele-defining positions.')


def get_allele_definitions(json_data: dict,
                           variants: dict[str, dict[str, str]]) -> dict[str, dict[str, str]]:
    """
    Extract allele definitions
    :param json_data: a dictionary of allele-defining variants and named alleles from PharmCAT JSON
    :param variants: a dictionary of all allele-defining positions
    :return: a dictionary of {allele_name: {
            allele-defining position (by hgvs names): {allele-defining genotype} }
        }
=    """
    # extract named allele definitions
    try:
        alleles: list[str] = [entry['name'] for entry in json_data['namedAlleles']]
        genotypes: list[list[str]] = [entry['alleles'] for entry in json_data['namedAlleles']]
        hgvs_names: list[str] = [*variants]
        n_variants = len(hgvs_names)
        n_alleles = len(alleles)

        # check whether there are any missing values
        if None in alleles:
            print('One of the allele names is missing')
        # check whether any allele definition is empty
        if any(x.count(None) == len(x) for x in genotypes):
            print('One of the alleles is not defined by any genotypes at any positions.')
        # check whether each allele is defined by the same number of positions in the json 'variants' attribute
        if any(len(g) != n_variants for g in genotypes):
            print('One of the alleles has a different number of allele-defining positions.')

        # generate a dictionary of allele definitions with the allele name and allele-defining genotypes
        allele_definitions: dict[str, dict[str, str]] = {}
        for i in range(n_alleles):
            allele_definitions[alleles[i]] = {}
            for j in range(n_variants):
                allele_definitions[alleles[i]][hgvs_names[j]] = genotypes[i][j]

        return allele_definitions

    except Exception as e:
        # print warnings if any of the attributes does not exist in the json file
        print(f'An error occurred in one of the {e}.'
              f'Check whether the attributes exist for all allele-defining positions.')


def count_allele_scores(dict_allele_definitions: dict[str, dict[str, str]]) -> dict[str, int]:
    """
    calculate the haplotype scores by counting the number of non-empty allele-defining positions
    :param dict_allele_definitions: a dictionary of allele definitions
        {allele_name: {hgvs_name: genotype}}
    :return: hap_score: dict[str, int] = { allele_name: allele_score }
    """
    # initialize the return variable
    allele_scores: dict[str, int] = dict()

    # iterate over alleles
    for allele, dict_defining_genotypes in dict_allele_definitions.items():
        # count the number of allele-defining genotypes
        genotypes = list(dict_defining_genotypes.values())
        n_defining_genotypes: int = len(genotypes) - genotypes.count(None)
        # add the allele score to allele_scores
        allele_scores[allele] = n_defining_genotypes

    return allele_scores


def count_diplotype_scores(allele_scores: dict[str, int]) -> dict[str, int]:
    """
    Based on a dictionary of allele/haplotype scores, get all the diplotype combinations and their scores
    :param allele_scores: a dictionary {allele_name: allele_score}
    :return: dip_score: a dictionary {diplotype_name: diplotype_score}
    """
    # calculate the diplotype scores
    diplotype_scores: dict[str, int] = {k1 + '/' + k2: v1 + v2
                                        for k1, v1 in allele_scores.items()
                                        for k2, v2 in allele_scores.items()}

    return diplotype_scores


def fill_definitions_with_references(dict_allele_definitions: dict[str, dict[str, str]],
                                     dict_allele_defining_variants: dict[str, dict[str, str]]) \
        -> dict[str, dict[str, str]]:
    """
    fill up empty allele-defining positions with the reference genotype
    :return: a list of genotypes where empty/None genotypes have been filled with reference genotypes at the position
    """
    # initialize the return variable
    dict_allele_definitions_filled: dict[str, dict[str, str]] = dict()

    # iterate over each allele's definition
    for allele, dict_defining_genotypes in dict_allele_definitions.items():
        # initialize an empty dictionary for an allele
        dict_allele_definitions_filled[allele] = {}
        # find the allele-defining genotype at each position
        for position, genotype in dict_defining_genotypes.items():
            # if the genotype at this position is None, use the reference genotype
            updated_genotype: str = dict_allele_defining_variants[position]['ref'] if genotype is None else genotype
            # add the update genotype to dict_allele_definitions_filled
            dict_allele_definitions_filled[allele][position] = updated_genotype

    return dict_allele_definitions_filled


def replace_wobble(wobl_genotype: str) -> Set[str]:
    """
    replace wobble genotypes with basic base pairs A/T/C/G
    :param wobl_genotype: an allele-defining genotype
    :return: a list of genotypes with only A, T, C, G, or indels
    """
    g_flat = _wobble_match_table[wobl_genotype]
    return g_flat


def get_unique_combinations(g1: str, g2: str) -> set[str]:
    """
    get all possible combinations of genotypes at an allele-defining positions
    :param g1: an allele-defining genotypes from allele 1
    :param g2: an allele-defining genotypes from allele 2
    :return: a set of unique combinations of genotypes
    """
    # replace wobbles
    g1_: list[str] = replace_wobble(g1) if g1 in _wobble_genotype_list else [g1]
    g2_: list[str] = replace_wobble(g2) if g2 in _wobble_genotype_list else [g2]

    # get all possible genotype combinations at this positions
    if None in g1_ and None in g2_:
        uniq_comb: Set[str] = set()
    elif None in g1_:
        uniq_comb: Set[str] = {x for x in g2_}
    elif None in g2_:
        uniq_comb: Set[str] = {x for x in g1_}
    else:
        uniq_comb: Set[str] = {'/'.join(sorted([x, y])) for x in g1_ for y in g2_}

    return uniq_comb


def get_diplotype_definition_dictionary(dict_allele_definitions: dict[str, dict[str, str]],
                                        alleles_to_test: list[str]) -> dict[str, dict[str, Set[str]]]:
    """
    obtain a dictionary of diplotype definitions
    :param dict_allele_definitions: a dictionary of allele definitions
    :param alleles_to_test: list of alleles that share allele-defining positions with other alleles
    :return: a dictionary {diplotype_name: {position: {unique combinations of genotypes at a position} } }
    """
    # initialize variables
    # set up an empty dictionary that stores each diplotype's definition
    diplotype_definitions: dict[str, dict[str, set[str]]] = {}

    # get the list of alleles
    allele_names: list[str] = [*dict_allele_definitions]
    n_alleles: int = len(allele_names)

    # iterate over the allele list
    for i in range(n_alleles):
        # get the allele name
        a1: str = allele_names[i]

        for j in range(i, n_alleles):
            # get the allele names
            a2: str = allele_names[j]

            # skip if both alleles are defined by exclusive positions
            if (a1 not in alleles_to_test) and (a2 not in alleles_to_test):
                continue

            # specify the diplotype's name, which will be used as the key in dic_mat
            diplotype_name = a1 + '/' + a2

            # third, create a dictionary variable to store the unique genotype combinations at each position
            single_diplotype_definition: dict[str, set[str]] = dict()
            for hgvs_name in dict_allele_definitions[a1]:
                g1 = dict_allele_definitions[a1][hgvs_name]
                g2 = dict_allele_definitions[a2][hgvs_name]
                single_diplotype_definition[hgvs_name] = get_unique_combinations(g1, g2)

            # finally, append the diplotype and its genotypes to the dictionary
            diplotype_definitions[diplotype_name] = single_diplotype_definition

    return diplotype_definitions


def find_powerset(s):
    """
    return the power set (all combination subsets) of input s
    :param s: a list of elements
    :return: power set of s
    """
    x = len(s)
    masks = [1 << i for i in range(x)]
    # use range(1, 1 << x) to exclude the empty set
    # otherwise, range(1 << x) to include the empty set
    for i in range(1, 1 << x):
        yield {ss for mask, ss in zip(masks, s) if i & mask}


def find_missingness_combinations(allele_array: np.ndarray,
                                  allele_names: list[str],
                                  hgvs_names: list[str],
                                  reference_alleles: list[str],
                                  alleles_to_test: list[str]) -> list[list[str]]:
    """
    Find the combinations of positions to be set to missing for autogenerated tests
    :param allele_array: an allele-position array (A, P) where each cell denotes whether
        the allele is defined by a specific genotype at a position
        A = number of alleles in a gene
        P = allele-defining positions
    :param allele_names: list of allele names of length A
    :param hgvs_names: list of hgvs names of length P
    :param reference_alleles: list of reference alleles
    :param alleles_to_test: list of alleles that share allele-defining positions with others
    :return: power set (all combinations) of missing positions
        if missing_position_combinations is empty, then there is no autogenerate test with missing positions
        if missing_position_combinations is not empty, it lists all missing combinations except the empty set
    """
    p_separator: str = ';;'

    # identify missing position combinations for each diplotype
    missing_position_combinations: list[list[str]] = []
    for a1 in allele_names:
        for a2 in alleles_to_test:
            # skip reference/reference
            if a1 in reference_alleles and a2 in reference_alleles:
                continue
            # otherwise, get the index of allele-defining positions
            # this can include positions exclusive to either alleles
            elif a1 in reference_alleles:
                i: int = allele_names.index(a2)
                p_idx: np.ndarray = np.where(allele_array[i, :])[0]
            elif a2 in reference_alleles:
                i: int = allele_names.index(a1)
                p_idx: np.ndarray = np.unique(np.where(allele_array[i, :])[0])
            else:
                i: int = allele_names.index(a1)
                j: int = allele_names.index(a2)
                p_idx: np.ndarray = np.unique(np.where(allele_array[[i, j], :])[1])

            # get the positions' hgvs names
            p_names: np.ndarray = np.array(hgvs_names)[p_idx]

            # find indices of positions that need to be missing together for an allele to be mistaken as another
            u, indices = np.unique(allele_array[:, p_idx], axis=1, return_inverse=True)
            # get the list of missing positions based on indices
            p_to_permute: set[str] = set()
            for k in range(indices.max() + 1):
                p_to_permute.add(p_separator.join(p_names[indices == k]))

            # find the powerset of missing positions
            combs: list[set[str]] = list(find_powerset(p_to_permute))

            # add combinations to missing_position_combinations
            for comb in combs:
                comb_split: list[str] = [s for x in comb for s in x.split(sep=p_separator)]
                if comb_split not in missing_position_combinations:
                    missing_position_combinations.append(comb_split)

    return missing_position_combinations


def convert_dictionary_to_numpy_array(dict_definitions: dict[str, dict]) -> np.ndarray:
    """
    convert the diplotype definition dictionary to a numpy array
    :param dict_definitions:
        {allele_name: {position_hgvs_name: nucleotide} }
        or
        {diplotype_name: {position_hgvs_name: nucleotide_combinations} }
    :return: a 3D numpy array [alleles/diplotypes, genotypes, positions]
    """
    # initialize variables
    # set up an empty set to record all the possible genotypes across positions in a gene
    unique_genotypes: set[str] = set()
    for dict_defining_genotypes in dict_definitions.values():
        # get the specific genotype(s) at a certain position
        for genotypes in dict_defining_genotypes.values():
            # for diplotypes, genotypes are a set of nucleotide combinations, like 'A/A' or 'delT/T'
            if isinstance(genotypes, set):
                # add each genotype to the list and skip None
                for x in genotypes:
                    if x is not None:
                        unique_genotypes.add(x)
            # for alleles, genotype(s) are a single nucleotide that defines an allele at a position
            elif isinstance(genotypes, str):
                unique_genotypes.add(genotypes)
            # for alleles, besides the specific allele-defining positions, others are reference and left as None
            elif genotypes is None:
                continue
            else:
                print(f'\tUnexpected data type was found in convert_dictionary_to_numpy_array')

    # convert the genotype set to a list which has an order
    genotype_list: list[str] = list(unique_genotypes)

    # get the list of allele or diplotype names
    definition_names: list[str] = [*dict_definitions]

    # get the list of hgvs names for all positions
    n_definitions: int = len(definition_names)
    hgvs_names: list[str] = [*dict_definitions[definition_names[0]]] if n_definitions > 0 else []

    # initialize an empty numpy array with rows and columns
    n_positions: int = len(hgvs_names)
    n_genotypes: int = len(genotype_list)
    definition_arrays: np.ndarray = np.zeros((n_definitions, n_genotypes, n_positions), dtype='int8')

    # fill the numpy arrays with 1 where position(key)-genotype(value) combinations exist
    for definition, dict_one_definition in dict_definitions.items():
        # get the allele or diplotype index for the numpy array
        idx_definition = definition_names.index(definition)

        for position, genotypes in dict_one_definition.items():
            # get the position index for the numpy array
            position_idx = hgvs_names.index(position)

            # skip None
            if genotypes is None:
                continue
            # otherwise, denote the corresponding genotype-position as '1' in the matrix
            else:
                for g in genotypes:
                    # get the genotype index for the numpy array
                    genotype_idx = genotype_list.index(g)

                    # fill the corresponding cell in the numpy array with 1 to denote the definition
                    definition_arrays[idx_definition, genotype_idx, position_idx] = 1

    return definition_arrays


def is_sharing_definitions(g1: np.ndarray, g2: np.ndarray, axis: Optional[int] = 1) -> np.ndarray[bool]:
    """
    to identify whether g1 (a diplotype's definition) shares definitions with another or any other diplotypes
    :param g1: a numpy array of a diplotype's definition (M, P).
                M = possible genotypes at a position for a diplotype
                P = allele-defining position for a gene
    :param g2: a numpy array for one diplotype (m,p) or multiple diplotypes (D, M, P)
                D = the total number of diplotypes for a gene
    :param axis: the axis at which g1 and g2 will be compared
            if g2 is a single diplotype, use axis = 0
            if g2 is a d x m x p array of multiple diplotypes, use axis = 1
    :return: a numpy array with boolean elements of length (1, ) or (D, )
    """

    # calculate product of g1 and g2
    g_prod = g1 * g2

    # check whether g1 shares any definitions with another diplotype by identifying non-empty columns
    # a non-empty column of g_prod means that g1 and this other diplotype shared a definition at a position
    non_empty_cols = np.any(g_prod != 0, axis=axis)

    # find diplotypes that share definitions with g1 at every allele-defining positions
    status = np.all(non_empty_cols, axis=axis)

    return status


def is_discrepant(g1: np.ndarray, g2: np.ndarray, axis: Optional[int] = 1) -> np.ndarray[bool]:
    """
    to identify whether g1 (a diplotype's definition) is different from another or any other diplotypes
    :param g1: a numpy array of a diplotype's definition (M, P).
                M = possible genotypes at a position for a diplotype
                P = allele-defining position for a gene
    :param g2: a numpy array for one diplotype (m,p) or multiple diplotypes (D, M, P)
                D = the total number of diplotypes for a gene
    :param axis: the axis at which g1 and g2 will be compared
            if g2 is a single diplotype, use axis = 0
            if g2 is a d x m x p array of multiple diplotypes, use axis = 1
    :return: a numpy array with boolean elements of length (1, ) or (D, )
\
    """
    # calculate product of g1 and g2
    g_prod = g1 * g2

    # check whether there is any zero column
    # an empty column means that g1 does not share defining genotypes with a certain diplotype
    # np.all(axis=axis) checks whether a position column has only empty entries
    empty_cols = np.all(g_prod == 0, axis=axis)

    # identify diplotypes that differ from g1 by identifying empty position columns
    status = np.any(empty_cols, axis=axis)

    return status


def is_equivalent(g1: np.ndarray, g2: np.ndarray) -> bool:
    """
    to identify whether g1 and g2 are equivalent
    :param g1: a numpy array of a diplotype's definition (M, P)
    :param g2: a numpy array of a diplotype's definition (M, P)
    :return: True or False
    """
    # compare d1 with d2 to find equivalent diplotypes
    # 'np.all(g1 == g2, axis=(1, 2))' is faster than '[numpy.array_equal(g1, g) for g in g2]'
    status = np.all(g1 == g2)

    return status


def is_included(g1: np.ndarray, g2: np.ndarray) -> bool:
    """
    to identify whether g1 is included in g2
    in other words, whether g1 can be considered as a sub-allele of g2

    :param g1: a numpy array of a diplotype's definition (M, P)
                M = possible genotypes at a position for a diplotype
                P = allele-defining position for a gene
    :param g2: a numpy array of a diplotype's definition (M, P)
    :return: True or False
    """
    # use numpy vectorized manipulation to check whether g1 is included as a sub-allele by any other diplotypes
    g_subs = g1 - g2

    # find diplotypes that includes g1
    status = set(np.unique(g_subs)) == {0, -1}

    return status


def is_inclusive(g1: np.ndarray, g2: np.ndarray) -> bool:
    """
    to identify whether g1 includes g2
    in other words, whether g1 can be considered as a super-allele of g2

    :param g1: a numpy array of a diplotype's definition (M, P)
                M = possible genotypes at a position for a diplotype
                P = allele-defining position for a gene
    :param g2: a numpy array of a diplotype's definition (M, P)
    :return: True or False
    """
    # use numpy vectorized manipulation to check whether g1 is included as a sub-allele by any other diplotypes
    g_subs = g1 - g2

    # find diplotypes that includes g1
    status = set(g_subs.flatten()) == {0, 1}

    return status


def is_overlapping(g1: np.ndarray, g2: np.ndarray) -> bool:
    """
    to identify whether g1 overlaps with g2
    g1 and g2 must have intersecting definitions, but g1 and g2 cannot be equivalent or a sub-allele of one the other

    :param g1: a numpy array of a diplotype's definition (M, P)
                M = possible genotypes at a position for a diplotype
                P = allele-defining position for a gene
    :param g2: a numpy array of a diplotype's definition (M, P)
    :return: True or False
    """
    # calculate g1*g2 and g1-g2
    g_prod = g1 * g2
    g_subs = g1 - g2

    # make sure g_prod does not have any empty columns
    # an empty column means that g1 and g2 do not share defining genotypes at a position
    non_empty_cols = np.any(g_prod != 0, axis=0)
    # if overlapping, two diplotypes must share defining genotypes at every position and there cannot be an empty column
    status_sharing_definitions = np.all(non_empty_cols)

    # find diplotypes that is not equivalent to, inclusive of, or included in g1
    status_overlapping = set(g_subs.flatten()) == {1, 0, -1}

    # check whether each diplotype have exclusive defining genotypes
    status = status_sharing_definitions * status_overlapping

    return status


def find_pairwise_relationship(d1: str, d2: str, g1: np.ndarray, g2: np.ndarray) -> str:
    """
    a wrapper function to determine the relationship between two diplotypes
    the input d1 and d2 must share definitions

    :param d1: the name of diplotype 1
    :param d2: the name of diplotype 2
    :param g1: the numpy array of diplotype 1 (M, P)
                M = possible genotypes at a position for a diplotype
                P = allele-defining position for a gene
    :param g2: the numpy array of diplotype 1 (M, P)
    :return: a string specifying the relationship between d1 and d2,
            including "equivalent", "overlapping", "included", and "inclusive"
    """
    # compare g1 with g2
    status_equivalent = is_equivalent(g1, g2)
    status_included = is_included(g1, g2)
    status_inclusive = is_inclusive(g1, g2)
    status_overlapping = is_overlapping(g1, g2)

    # g1 and g2 must have only one True status
    status_sum = sum([status_equivalent, status_included, status_inclusive, status_overlapping])
    relationship: str = ''
    if status_sum != 1:
        print(f'Something is wrong between {d1} and {d2}.')
        sys.exit(1)
    elif status_equivalent:
        relationship = 'equivalent'
    elif status_included:
        relationship = d1 + ' is included in ' + d2
    elif status_inclusive:
        relationship = d1 + ' is inclusive of ' + d2
    else:
        relationship = 'overlapping'

    return relationship


def find_possible_calls(g1: np.ndarray,
                        definition_arrays: np.ndarray,
                        diplotype_names: list[str]) -> list[str]:
    """
    find all alternative calls for one diplotype
    :param g1: the numpy array of diplotype 1 (M, P)
                M = possible genotypes at a position for a diplotype
                P = allele-defining position for a gene
    :param definition_arrays: a numpy array of all diplotype definitions (D, M, P)
                D = the total number of diplotypes for a gene
    :param diplotype_names: a list of diplotype names
    :return: a list of diplotypes that share definitions with g1 and can show up as alternative calls for g1
    """
    # check whether g1 shares definitions with any other diplotype definitions
    status_sharing_definitions: np.ndarray = is_sharing_definitions(g1, definition_arrays)
    # # check whether each of the diplotypes differs from g1
    # status_discrepant: np.ndarray = is_discrepant(g1, definition_arrays)
    # # sanity check to make sure no diplotypes can be discrepant but share definitions at the same time
    # if (status_sharing_definitions == status_discrepant).any():
    #     print(f'Something went wrong. Some diplotypes are different but also the same at the same time.')
    #     sys.exit(1)
    # # continue to the next diplotype if the current one is discrepant from all other diplotypes
    # if all(status_discrepant):
    #     print(f'Something went wrong. A diplotype should at least share the allele definition with itself.')
    #     sys.exit(1)

    # get the names of the diplotypes that share definitions with g1
    diplotypes_sharing_definitions: list[str] = [d for d, s in zip(diplotype_names, status_sharing_definitions) if s]

    return diplotypes_sharing_definitions


def find_wobble_subsets(genotypes: np.ndarray, diplotypes: np.ndarray, wobble_subsets: set[frozenset[str]]) \
        -> set[frozenset[str]]:
    """
    find subsets of possible calls for provided diplotypes
    :param genotypes: a numpy array of multiple diplotypes (D, M, P)
            D = the total number of diplotypes for a gene
            M = all possible genotypes at a position for a diplotype
            P = allele-defining position for a gene
    :param diplotypes: an 1D array of diplotypes (D, )
    :param wobble_subsets: a set that saves all possible calls based on the input genotype matrix
    :return: subsets of diplotypes that share unique genotype definitions
    """
    # establish the terminating condition
    # if there are only two or less diplotypes, no need to find subsets
    if len(diplotypes) < 3:
        return wobble_subsets
    elif genotypes.shape[0] != len(diplotypes):
        print(f'genotype matrix does not match the number of diplotypes.')
    else:
        # find the places in genotype_matrix where only a subset of diplotypes share a certain definition
        genotype_idx, position_idx = np.where(np.any(genotypes, axis=0) ^ np.all(genotypes, axis=0))

        # identify unique subsets
        # (D, S) array where S = # of unique subsets
        # specify the data type of subsets_idx as boolean to properly extract the corresponding columns later
        subsets_idx = np.unique(genotypes[:, genotype_idx, position_idx], axis=1).astype(bool)
        n_subsets = subsets_idx.shape[1]

        # process each subset of diplotype
        for i in range(n_subsets):
            diplotype_idx: np.ndarray[bool] = subsets_idx[:, i]
            diplotype_subset: np.ndarray[int] = diplotypes[diplotype_idx]

            # continue to the next if the current subset is already recorded in the wobble subsets
            if frozenset(diplotype_subset) in wobble_subsets:
                continue
            else:
                # add the current subset to the summary wobble_subsets
                wobble_subsets.add(frozenset(diplotype_subset))

                # continue to find granular subsets
                genotypes_subset: np.ndarray[int] = genotypes[diplotype_idx, :, :]
                find_wobble_subsets(genotypes_subset, diplotype_subset, wobble_subsets)

        return wobble_subsets


def find_all_possible_calls(definition_arrays: np.ndarray, diplotype_names: list[str],
                            possible_call_sets: list[set[str]]) -> list[set[str]]:
    # initialize variables
    n_positions = definition_arrays.shape[2]
    # initialize an empty list to store the sets of overlapping diplotypes based on the current position
    current_possible_call_sets: list[set[str]] = []

    for idx_position in range(n_positions):
        # get the total number of possible_call_sets
        n_combinations = len(possible_call_sets)

        # update the overlapping diplotype calls based on the current position
        for i in range(n_combinations):
            # get the possible calls
            possible_calls: np.ndarray = np.array(list(possible_call_sets[i]))
            idx_possible_calls: list[int] = [diplotype_names.index(x) for x in possible_calls]
            genotypes_of_possible_calls = definition_arrays[idx_possible_calls, :, idx_position]

            # find genotypes that define more than one diplotype
            busy_genotypes = np.where(np.sum(genotypes_of_possible_calls, axis=0) > 1)[0]
            if len(busy_genotypes):
                for idx_genotype in busy_genotypes:
                    idx_diplotypes = np.unique(np.where(genotypes_of_possible_calls[:, idx_genotype])[0])
                    current_possible_calls = set(possible_calls[idx_diplotypes])
                    if current_possible_calls not in current_possible_call_sets:
                        current_possible_call_sets.append(current_possible_calls)

        # update the possible_call_sets
        possible_call_sets = current_possible_call_sets
        # reset current_possible_call_sets
        current_possible_call_sets = []

    # return the summary dictionary of alternative calls for all diplotypes
    return possible_call_sets


def find_predicted_calls(possible_call_sets: list[set[str]],
                         diplotype_scores: dict[str, int],
                         missing_positions: list[str] = None) -> dict[list[str], list[str], list[str]]:
    """
    for a diplotype that have alternative calls, find the call with the highest score that PharmCAT will return
    :param possible_call_sets: (list) sets of possible calls.
        Each set comprises diplotypes with overlapping definitions.
    :param diplotype_scores: dictionary of diplotype scores based on the number of allele-defining positions
    :param missing_positions: list of positions that are presumed missing and not considered in this round of comparison
    :return: dictionary of three lists of the same length
        expected calls: purported genotypes for a diplotype in a VCF
        actual calls: diplotypes with the highest score based on the number of allele-defining positions
        alternative calls: other diplotypes that share definitions with the actual calls but have lower scores
    """
    # initialize the return dictionary
    predict_calls_cols = ['expected', 'actual', 'alternative', 'missing_positions']
    predicted_calls: dict = {key: [] for key in predict_calls_cols}

    # iterate over possible_call_set
    for possible_call_set in possible_call_sets:
        if len(possible_call_set) == 0:
            continue
        # initialize lists to concatenate
        actual_calls: list[str] = []
        alternative_calls: list[str] = []

        # get diplotype scores for the possible call
        scores = np.array([diplotype_scores[d] for d in possible_call_set])
        # get the diplotypes with the max scores
        idx_max_score = np.where(scores == np.max(scores))[0]

        # separate the actual and alternative calls based on diplotype scores
        for i, d in enumerate(possible_call_set):
            (alternative_calls, actual_calls)[i in idx_max_score].append(d)

        # add results to the predicted_calls
        for d in possible_call_set:
            # add to the summary dictionary
            predicted_calls['expected'].append(d)
            predicted_calls['actual'].append(';'.join(actual_calls))
            predicted_calls['alternative'].append(';'.join(alternative_calls))
            predicted_calls['missing_positions'].append(';'.join(missing_positions) if missing_positions else None)

    # return the expected, actual, and alternative calls
    return predicted_calls


def predict_pharmcat_calls(dict_allele_defining_variants: dict[str, dict[str, str]],
                           dict_allele_definitions: dict[str, dict[str, set[str]]],
                           alleles_to_test: list[str],
                           missing_positions: Optional[set[str]] = None) -> dict[list[str], list[str], list[str]]:
    # initialize values
    definitions = dict()
    defining_variants = dict()

    # remove missing positions from the list of allele defining positions
    for position in dict_allele_defining_variants:
        # skip if position is one of the missing positions
        if missing_positions:
            if position in missing_positions:
                continue
        # add the position to the new dictionary
        defining_variants[position] = {
            'position': dict_allele_defining_variants[position]['position'],
            'ref': dict_allele_defining_variants[position]['ref']
        }

    # remove missing positions from allele definitions
    for allele, dict_defining_genotypes in dict_allele_definitions.items():
        if missing_positions:
            # get the hgvs names of allele-defining positions, excluding missing positions
            p: list[str] = [k for k in dict_defining_genotypes if k not in missing_positions]
            # get genotypes at each allele-defining position, excluding missing positions
            g: list[str] = [dict_defining_genotypes[hgvs_name] for hgvs_name in p]
            # if definitions are not empty for this allele, add to the definitions dictionary
            if any(g):
                definitions[allele] = {k: v for k, v in zip(p, g)}
        else:
            definitions[allele] = dict_defining_genotypes

    # update the allele_definitions and fill empty cells with reference genotypes
    allele_definitions_filled: dict[str, dict[str, str]] = fill_definitions_with_references(definitions,
                                                                                            defining_variants)

    # get a dictionary of diplotype definitions for comparison
    # {diplotype_name: {position: {unique combinations of genotypes at a position} } }
    diplotype_definitions = get_diplotype_definition_dictionary(allele_definitions_filled, alleles_to_test)

    # get diplotype definition arrays
    # vectorized computation by numpy arrays speeds up diplotype comparison
    diplotype_definition_arrays: np.ndarray = convert_dictionary_to_numpy_array(diplotype_definitions)
    # get diplotype names
    diplotype_names: list[str] = [*diplotype_definitions]
    # initialize a possible call sets
    possible_call_sets: list[set[str]] = [set(diplotype_names)]
    # find possible calls for each diplotype
    possible_call_sets = find_all_possible_calls(diplotype_definition_arrays, diplotype_names, possible_call_sets)

    # get a dictionary of calculated diplotype scores based on the number of core allele-defining positions
    # first, calculate the haplotype scores
    allele_scores = count_allele_scores(definitions)
    # then, calculate the diplotype scores by adding up the allele scores
    diplotype_scores = count_diplotype_scores(allele_scores)
    # summarize expected vs actual calls
    dict_predicted_calls = find_predicted_calls(possible_call_sets, diplotype_scores, missing_positions)

    # compare diplotype pairs and identify pairwise relationships
    # pairwise_comparisons = compare_diplotype_pairs(diplotype_definitions, allele_defining_variants)

    return dict_predicted_calls


def compare_diplotype_pairs(dict_diplotype_definitions: dict[str, dict[str, set[str]]]) -> dict:
    """
    compare unphased diplotype definitions
    G = a MxP matrix where rows are all possible nucleotide pairs in a gene (M) and columns are positions (P)

    Possible outcomes of diplotype comparison
    1. Sharing definitions
        - between any diplotypes (G1, G2)
        - criteria:
            1) G1 * G2 = G3. G3 cannot have empty columns.
                empty column = G1 and G2 do not share any definitions at a position

        Four possible relationships between diplotypes that share definitions
        1. equivalent
            - between any diplotypes (G1, G2)
            - criteria:
                1) G1 shares a definition with G2
                2) G1 - G2 = G3 = {0}
            - alternative criteria:
                1) G1 * G2 = G3. G3 cannot have empty columns.
                2) G1 - G2 = G4 = {0}

        2. included
            - G1 (a wobble or a non-wobble diplotype) is included in G2 (another wobble diplotype) as if a sub-allele
            - criteria:
                1) G1 shares a definition with G2
                2) G1 != G2 (not equivalent)
                3) G1 <= G2
            - alternative criteria
                1) G1 * G2 = G3. G3 cannot have empty columns.
                2) G1 - G2 = G4 where G4 must be made up by {0, -1}

            criteria:
            1) (integer operation) G1 != G2 and G1 <= G2
            1) (alternative) G1 - G2 = G3 with only {0,-1}. G3 cannot have empty columns.

        3. inclusive
            - G1 (a wobble diplotype) includes G2 (other wobble and non-wobble diplotypes)
            - criteria:
                1) G1 shares a definition with G2
                2) G1 != G2 (not equivalent)
                3) G1 >= G2
            - alternative criteria
                1) G1 * G2 = G3. G3 cannot have empty columns.
                2) G1 - G2 = G4 where G4 must be made up by {0, -1}

        4. overlapping
            - G1 (a wobble diplotype) overlaps with G2 (another wobble diplotype)
            - criteria:
                1) G1 shares a definition with G2
                2) G1 != G2 (not equivalent)
                3) G1 is not greater than, nor less than, G2
            - alternative criteria
                1) G1 * G2 = G3. G3 cannot have empty columns.
                2) G1 - G2 = G4 where G4 must be made up by {1, 0, -1}

    2. discrepant
        - between any diplotypes (G1, G2)
        - criteria:
            1) G1 - G2 = G3 where G3 must be made up by {1, 0, -1}
            2) G1 * G2 = G4. G4 has empty columns
                empty column = G1 and G2 do not share any defining genotypes at a position


    :param dict_diplotype_definitions:
    :return: a dictionary {
        'diplotype_1': {
            score: 123,
            diplotype_name: relation
        } }
    """
    # set up an empty dictionary to store diplotype comparison results
    pairwise_comparison = dict()

    # convert the diplotype definition dictionaries to numpy arrays of genotypes for each diplotype
    # vectorized computation by numpy arrays speeds up diplotype comparison
    definition_arrays: np.ndarray = convert_dictionary_to_numpy_array(dict_diplotype_definitions)

    # get the list of diplotype names
    diplotype_names: list[str] = [*dict_diplotype_definitions]
    # get the number of diplotypes
    n_diplotypes = len(dict_diplotype_definitions)

    # iterate over diplotypes
    for i in range(n_diplotypes):

        # get the genotype array of the target diplotype
        d1 = diplotype_names[i]
        g1 = definition_arrays[i, :, :]

        # check whether g1 shares definitions with any other diplotype definitions
        status_sharing_definitions: np.ndarray = is_sharing_definitions(g1, definition_arrays)
        # check whether each of the diplotypes differs from g1
        status_discrepant = is_discrepant(g1, definition_arrays)
        # sanity check to make sure no diplotypes can be discrepant but share definitions at the same time
        if (status_sharing_definitions == status_discrepant).any():
            print(f'Something went wrong. Some diplotypes are different but also the same at the same time.')
            sys.exit(1)
        # continue to the next diplotype if the current one is discrepant from all other diplotypes
        if status_discrepant.all():
            print(f'Something went wrong. A diplotype should at least share the allele definition with itself.')
            sys.exit(1)

        # specify the relationship between diplotypes that share definitions
        for j, s in enumerate(status_sharing_definitions):
            if i == j:  # skip self
                continue
            elif s:  # s == True meaning two diplotypes share definitions
                # get the definition array of another diplotype
                d2 = diplotype_names[j]
                g2 = definition_arrays[j, :, :]

                # find the relationship between d1 and d2
                relationship = find_pairwise_relationship(d1, d2, g1, g2)

                # add diplotype and its score to the output dictionary
                if d1 not in pairwise_comparison:
                    pairwise_comparison[d1] = {}

                # add the relationship between g1 and g2 to a dictionary
                pairwise_comparison[d1][d2] = relationship

            else:
                continue

    return pairwise_comparison


def order_list(list_to_be_ordered: list[str], reference_list: list[str]) -> list[str]:
    """
    order a list based on another list
    :param list_to_be_ordered:
    :param reference_list:
    :return:
    """
    return [x for x in reference_list if x in list_to_be_ordered]
