#!/usr/bin/env python3

"""
Usage: test_gen_missing.py <definition.json> <positions.vcf> <vcfOutputDir>
Author Scott Dudek

This tool generates test cases for PharmCAT by reading a single gene definition JSON file
and the PharmCAT positions file. It then writes a set of synthetic VCF files with their
expected haplotype calls for a variety of scenarios.

Tests include:
* Homozygous (two copies of the same allele) for all alleles which are defined by more
  than one variant
* Heterozygous with all other alleles for each allele defined by more than one
  position
* In all cases above, each defining variant set as missing in a series of test files. For
alleles with more than two defining position, create a test for each combination
where all members of the combination are missing.
* Alleles with ambiguous codes will have additional tests containing the wobble

Note:
This script checks the number of files that will be produced and restricts the missing
position combination size to keep the number of files <= 25k or uses only single
positions if no combination size meets that criterion. This restriction currently only
affects CYP2D6 which would produce 313,478,400 files if unrestricted. Therefore, it only
produces VCF files for CYP2D6 with a single position missing and still produces 295,220
files.

Output:

Files have the alleles in the VCF file and the term wobble if at least one of the alleles
has a wobble. The missing position(s) are identified by the indexes in the VCF.

Sample output:
NUDT15 test file homozygous for *2 with one position defining *2 set to missing (index 4)
NUDT15_s2_s2_t1_m4.vcf

CYP3A5 test file heterozygous for *3/*5 with a wobble and 2 positions identifying
*3 set to missing ( indexes 0 and 5 )
CYP3A5_s3_s5_wobble8_m0.5.vcf

"""

import datetime
import itertools
import json
import operator as op
import os
import sys
from functools import reduce

import test_gen_utilities as util


# TODO: a proper argparse handler!
if len(sys.argv) != 4:
    sys.exit("ERROR: Expecting 3 arguments\n")

if not os.path.isfile(sys.argv[1]):
    sys.exit(f"ERROR: Cannot find JSON file '{sys.argv[1]}'\n")

if not os.path.isfile(sys.argv[2]):
    sys.exit(f"ERROR: Cannot find VCF file '{sys.argv[2]}'\n")

# maximum missing combination size used in producing test files
comboLimit = 16
# maximum file number
fileLimit = 250000

# functions for calculating file total
def ncr(n, r):
    """
    Takes the number of objects and the sample size
    Returns number of combinations
    """
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer // denom
# ncr()

def calcFileTotal(nas, ref, totalTests, maxcombo):
    """
    Takes the number of objects and the sample size.
    Total files = sum (number of tests for each allele > 1 position *
    combinations of missing) * total tests for all alleles
    Returns number of files that will be written
    """
    totalFiles = 0
    for na in nas:
        npos = len(na['_mindef'])
        if npos < 2 or na == ref:
            continue
        combosum = 0
        climit = min(npos-1, maxcombo)
        ntests = 0
        for sz in range(1,climit+1):
            combosum += ncr(npos,sz)
        ntests += len(na['_tests']) * combosum
        totalFiles += ntests * totalTests
    return totalFiles
# calcFileTotal()

# load the named allele definitions
defPath = sys.argv[1]

print(f"Loading '{defPath}'...")
with open(defPath, 'r') as defFile:
    definition = json.load(defFile)
# with defFile

# update indel variants with information from positions VCF
vcfPath = sys.argv[2]
print(f"Loading '{vcfPath}'...")
vcfreference = util.parseVCF(vcfPath)
util.updateVarAlleles(definition,vcfreference)


print(f"done: {len(definition['variants'])} variants, {len(definition['namedAlleles'])} named alleles\n")

# scan named alleles and identify the referent named allele
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
print("done\n")


# validate referent allele
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

# generate tests
totalTests = 0
print("Generating test cases...")
for namedallele in namedalleles:
    namedallele['_tests'] = list()
    # generate tests for each namedallele with every position defined and referent used
    # for any position that is NULL in the defintion. Wobbles are expanded so that
    # a namedallele may have multiple tests
    namedallele['_maxdef'] = {
        i: namedallele['_mindef'].get(i, refNamedallele['_mindef'][i]) for i in range(len(definition['variants']))
    }
    for alleles in util.expandAlleles(namedallele['_maxdef'].items(), 16 ):
        namedallele['_tests'].append(frozenset(alleles))
    totalTests += len(namedallele['_tests'])
    # for alleles
# for namedallele

# check referent allele against VCF positions file and correct where needed for ambiguous
util.checkRefNamedAllele(definition, vcfreference, refNamedallele)

# select maximum combination size within limits
# will always use at least single missing position for tests (as for CYP2D6)
maxsize = comboLimit
for sz in reversed(range(1,comboLimit+1)):
    sumtests = calcFileTotal(namedalleles, refNamedallele, totalTests, sz)
    maxsize = sz
    if sumtests < fileLimit:
        break

numFiles = 0
numMultiPos = 0
testtype = "t"

basepath = sys.argv[3]
if not os.path.exists(basepath):
    os.makedirs(basepath)
print("Writing files...")

for i,na1 in enumerate(namedalleles):
    # if a namedallele has defining positions > 1, then it will have output files
    # skip the ref alllele
    if len(na1['_mindef']) < 2 or na1 == refNamedallele:
        continue
    numMultiPos += 1
    sizes = [i for i in range(1, min(maxsize, len(na1["_mindef"])))] if maxsize > 1 else [1]

    test1alleles = list(dict(test) for test in na1["_tests"])
    matches1 = {i,}

    # create VCF files by pairing with all alleles
    for j,na2 in enumerate(namedalleles):
        test2alleles = list(dict(test) for test in na2["_tests"])
        matches2 = (j,)
        if len(test1alleles) > 1 or len(test2alleles) > 1:
            testtype="wobble"
        else:
            testtype="t"
        numTest = 1

        for alleles1, alleles2 in itertools.product(test1alleles, test2alleles):
            numCombos = 0
            for size in sizes:
                # missing combinations are in reference to the first named allele
                for combo in itertools.combinations(na1['_mindef'].keys(),size):
                    numCombos += 1

                    outPath = os.path.join(basepath, util.percentEncode("%s_%s_%s_%s%s_m%s.vcf" % (
                        definition['gene'],
                        util.joinMatches(matches1, namedalleles),
                        util.joinMatches(matches2, namedalleles),
                        testtype,
                        numTest,
                        ".".join([str(definition['variants'][i]['position']) for i in sorted(combo)])
                    )))
                    # skip test cases if all variants will be missing from second allele
                    if set(na2['_mindef'].keys()).issubset(combo):
                      continue

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
                            if(len(varalleles)==1 and "ALT" in vcfreference[definition["gene"]][variant['_vcfidx']]):
                                varalleles.append(vcfreference[definition["gene"]][variant['_vcfidx']]["ALT"])
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
                                (varalleles.index(a1) if not (any(v == idx for idx in combo)) else "."),
                                (varalleles.index(a2) if not (any(v == idx for idx in combo)) else ".")
                            ))
                            outFile.write("\t".join(row) + "\n")
                        # for v,variant
                    # with outFile
                # for combo
            # for size
            numTest += 1
        #for alleles1, alleles2
    #for j,na2
#for i,na1


print(f"Alleles defined by multiple positions: {numMultiPos}")
print(f"Done: {numFiles} files")

