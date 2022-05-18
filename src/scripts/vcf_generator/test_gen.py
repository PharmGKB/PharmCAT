#!/usr/bin/env python3

"""
Usage: test_gen.py <definition.json> <vcfOutputDir>
Authors Alex Frase, Scott Dudek

This tool generates test cases for PharmCAT by reading a single gene definition JSON file
and the PharmCAT positions file. It then writes a set of synthetic VCF files with their
expected haplotype calls for a variety of scenarios.

Tests include:
* Homozygous (two copies of the same allele) complete (all basepair positions included and
  non-missing) for each allele; *1/*1,*2/*2, ...
* Heterozygous complete cases for each pair of alleles; *1/*2, *1/*3, *2/*3, ...
* Alleles with wobble (ambiguous nucleic codes) are expanded into separate tests for each
  of the cases above.

Output files are VCFv4.3 and have a single sample each. The call for the sample is
reflected in the filename and the "##PharmCATnamedAlleles" meta information header). If
an allele contains ambiguous codes it, multiple files will be generated for that call with
a different variation of the wobble in each. Such files will have "wobble" appended to the
name.

Example filenames:
* no wobble for CYP3A5 *1/*7
CYP3A5_s1_s7_t1.vcf

* CYP3A5 *3 has multiple ambiguous codes resulting in 16 files when paired with *7
CYP3A5_s3_s7_wobble1.vcf, CYP3A5_s3_s7_wobble2.vcf, ... CYP3A5_s3_s7_wobble16.vcf

"""

# ask: definition specifies genomeBuild "b38", vcf seems to expect "hg38"
# ask: does the VCF need 'source' meta info? old test files have it

# ASK: VCFv4.3 ALT is defined as
#          "base Strings made up of the bases A,C,G,T,N (case insensitive)
#          or the '*'  symbol (allele missing due to overlapping deletion)
#          or a MISSING value '.' (no variant)"
#      some gene definitions have deletions with a length that would run
#      over subsequent variants in the definition; how is this handled?
#      will PharmCAT parse the '*' symbol correctly in that context?

import collections
import datetime
import itertools
import json
import os
import sys

import test_gen_utilities as util


# TODO: a proper argparse handler!
if len(sys.argv) != 3:
    sys.exit("ERROR: Expecting 2 arguments\n")

if not os.path.isfile(sys.argv[1]):
    sys.exit(f"ERROR: Cannot find JSON file '{sys.argv[1]}'\n")

# load the named allele definitions
defPath = sys.argv[1]

quietMode = os.environ.get('PHARMCAT_TEST_QUIET') == "true"

if not quietMode:
    print(f"Loading '{defPath}'...")
with open(defPath, 'r') as defFile:
    definition = json.load(defFile)
# with defFile

if not quietMode:
    print(f"done: {len(definition['variants'])} variants, {len(definition['namedAlleles'])} named alleles\n")

# scan named alleles and identify the referent named allele
if not quietMode:
    print("Scanning named alleles...")
namedalleles = list()
refNamedallele = None
for namedallele in definition['namedAlleles']:
    if len(namedallele['alleles']) != len(definition['variants']):
        sys.exit("ERROR: named allele '%s' has %d variants, expected %d" % (
            namedallele['name'], len(namedallele['alleles']), len(definition['variants'])
        ))
    namedallele['_num'] = len(namedalleles)
    namedallele['_mindef'] = {i: a for i, a in enumerate(namedallele['alleles']) if a is not None}
    namedallele['_indecies'] = set(namedallele['_mindef'].keys())
    if not namedallele['_mindef']:
        print(f"  WARNING: null definition for named allele '{namedallele['name']}'")
        continue
    namedalleles.append(namedallele)
    # check if it's the referent named allele
    if namedallele['reference']:
        refNamedallele = namedallele
    # if ref?
# for namedallele
if not quietMode:
    print("done\n")


# validate referent allele
if not quietMode:
    print("Checking nucleic code notations...")
for i, a in refNamedallele['_mindef'].items():
    if a.startswith('del'):
        if a != 'del':
            sys.exit("ERROR: named allele '%s' invalid deletion #%d '%s'" % (refNamedallele['name'], i + 1, a))
        refNamedallele['_mindef'][i] = '.'
    elif a.startswith('ins'):
        sys.exit("ERROR: named allele '%s' invalid insertion #%d '%s'" % (refNamedallele['name'], i + 1, a))
    else:
        if any((s not in util.symbolBases) for s in a):
            sys.exit("ERROR: named allele '%s' invalid allele #%d '%s'" % (refNamedallele['name'], i + 1, a))
    # if del/ins/swap
# for i,a in _mindef

# validate other alleles
for namedallele in namedalleles:
    if namedallele is refNamedallele:
        continue
    for i, a in namedallele['_mindef'].items():
        if a.startswith('del'):
            if (a != 'delGene') and (a[3:] != refNamedallele['_mindef'][i]):
                sys.exit("ERROR: named allele '%s' invalid deletion #%d '%s'" % (refNamedallele['name'], i + 1, a))
            namedallele['_mindef'][i] = '.'
        else:
            if a.startswith('ins'):
                if refNamedallele['_mindef'][i] != '.':
                    sys.exit("ERROR: named allele '%s' invalid insertion #%d '%s'" % (refNamedallele['name'], i + 1, a))
                a = a[3:]
            # if ins
            if any((s not in util.symbolBases) for s in a):
                sys.exit("ERROR: named allele '%s' invalid allele #%d '%s'" % (refNamedallele['name'], i + 1, a))
            namedallele['_mindef'][i] = a
        # if del/ins/swap
    # for i,a in _mindef
# for namedallele

# check referent allele against defined reference in json and correct where needed for ambiguous
util.checkRefNamedAllele(definition, refNamedallele)

# generate tests
if not quietMode:
    print("Generating test cases...")
for namedallele in namedalleles:
    namedallele['_tests'] = list()
    # generate tests for each namedallele with every position defined and referent used
    # for any position that is NULL in the defintion. Wobbles are expanded so that
    # a namedallele may have multiple tests
    namedallele['_maxdef'] = {
        i: namedallele['_mindef'].get(i, refNamedallele['_mindef'][i]) for i in range(len(definition['variants']))
    }
    for alleles in util.expandAlleles(namedallele['_maxdef'].items(), util.optNumAlleleExpansions):
        namedallele['_tests'].append(frozenset(alleles))
    # for alleles
# for namedallele

basepath = sys.argv[2]
if not os.path.exists(basepath):
    os.makedirs(basepath)
if not quietMode:
    print("Writing files...")
numFiles=0
for i,na1 in enumerate(namedalleles):
    test1alleles = list(dict(test) for test in na1["_tests"])
    matches1 = (i,)

    for j,na2 in itertools.islice(enumerate(namedalleles),i, None):
        test2alleles = list(dict(test) for test in na2["_tests"])
        matches2 = (j,)

        combos = collections.defaultdict(list)
        testtype="t"
        if len(test1alleles) > 1 or len(test2alleles) > 1:
            testtype="wobble"
        numTest=1
        for alleles1, alleles2 in itertools.product(test1alleles, test2alleles):
            tall = tuple(sorted(alleles2.items()))
            # skip if combination already done
            if any(x == alleles1 for x in combos[tuple(sorted(alleles2.items()))]):
                continue
            combos[tuple(sorted(alleles1.items()))].append(alleles2)

            outPath = os.path.join(basepath, util.percentEncode("%s_%s_%s_%s%s.vcf" % (
                definition['gene'],
                util.joinMatches(matches1, namedalleles),
                util.joinMatches(matches2, namedalleles),
                testtype,
                numTest
            )))
            numTest+=1

            with open(outPath, 'w') as outFile:
                numFiles += 1
                # write meta-info header
                outFile.write("##fileformat=VCFv4.3\n")
                outFile.write("##fileDate=%s\n" % (datetime.date.today().strftime("%Y%m%d"),))
                # TODO: proper hg# notation?
                outFile.write("##reference=%s\n" % (definition['genomeBuild'].replace("b", "hg"),))
                outFile.write("##PharmCATnamedAlleles=%s/%s\n" % (
                    ('_'.join(namedalleles[na]['name'] for na in matches1) or "?"),
                    ('_'.join(namedalleles[na]['name'] for na in matches2) or "?"),
                ))
                outFile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
                # based on the spec there should be no need to define this explicitly
                # outFile.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")

                row = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT","TEST" + str(numTest)]
                outFile.write("\t".join(row) + "\n")
                for v, variant in enumerate(definition['variants']):
                    refAllele = refNamedallele['_maxdef'].get(v)
                    if not refAllele or refAllele == ".":
                        continue
                    varalleles = [refNamedallele['_maxdef'].get(v)]
                    for alleles in [alleles1, alleles2]:
                        allele = alleles.get(v,None)
                        alt = None if (allele is None) or allele == "." else allele
                        if alt and (alt not in varalleles):
                            varalleles.append(alt)
                    # for alleles
                    if(len(varalleles)==1 and definition['variants'][v]['alts']):
                        varalleles.extend(definition['variants'][v]['alts'])
                    row = [
                        str(variant.get('chromosome') or "."),
                        str(variant.get('position') or "."),
                        str(variant.get('rsid') or "."),
                        refAllele,
                        (",".join(varalleles[1:]) if (len(varalleles) > 1) else "."),
                        ".",  # QUAL
                        "PASS",  # FILTER
                        ".",  # INFO
                        "GT",  # FORMAT
                    ]
                    a1 = alleles1.get(v)
                    a2 = alleles2.get(v)
                    row.append("%s/%s" % (
                        (varalleles.index(a1) if (a1 and a1 != ".") else "."),
                        (varalleles.index(a2) if (a2 and a2 != ".") else "."),
                    ))
                    outFile.write("\t".join(row) + "\n")
                # for v,variant
            # with outFile
        # for alleles1, alleles2
    # for j,na2
# for i,na1

if not quietMode:
    print(f"Done: {numFiles} files")
