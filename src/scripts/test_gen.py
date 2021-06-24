#!/usr/bin/env python3

"""
Usage: test_gen.py <definition.json> <positions.vcf> <vcfOutputDir>
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

# TODO: once available, use definition's specified referent named allele rather than just the first fully-defined named allele

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
import csv
import datetime
import itertools
import json
import os
import re
import sys


# TODO: a proper argparse handler!
if len(sys.argv) != 4:
    sys.exit("ERROR: Expecting 3 arguments\n")

if not os.path.isfile(sys.argv[1]):
    sys.exit(f"ERROR: Cannot find JSON file '{sys.argv[1]}'\n")

if not os.path.isfile(sys.argv[2]):
    sys.exit(f"ERROR: Cannot find VCF file '{sys.argv[2]}'\n")

# set options
optNumAlleleExpansions = 128

# build basepair symbol lookup tables
symbolBases = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
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
    'N': {'A', 'C', 'G', 'T'},
    '.': {'.'},
}
symbolBaselist = {s: tuple(sorted(symbolBases[s])) for s in symbolBases.keys()}
basesSymbol = {frozenset(bases): symbol for symbol, bases in symbolBases.items()}


# define VCF parsing and JSON update functions

def parseVCF(vcf):
    """
    Takes the path to the pharmCAT VCF file containing all positions
    Returns a
     with gene names (from 'INFO' column) as key
    and a list of dicts (from rows) as values
    """
    start=0
    with open(vcf) as vcffile:
        for index,line in enumerate(vcffile):
            if line.startswith('#CHROM'):
                start=index
                break
        # for index,line
    # open vcf

    vcfRef=collections.defaultdict(list)

    pxpattern=re.compile('PX=(\w+)')
    with open(vcf) as vcffile:
        reader = csv.DictReader(itertools.islice(vcffile,start,None), delimiter='\t')
        for row in reader:
            m=pxpattern.match(row['INFO'])
            if m is not None:
                vcfRef[m.group(1)].append(row)
        # for row
    # open vcf
    return vcfRef
# parseVCF()


def updateVarAlleles(jsondef, vcfref):
    """
    Takes JSON dict and dict from parsed VCF
    Updates JSON dict with corrected indel information obtained from the VCF dict
    """
    gene=jsondef['gene']
    if gene not in vcfref:
        sys.exit("ERROR: %s is not present in VCF file" % gene)
    if len(vcfref[gene]) != len(jsondef["variantAlleles"]):
        sys.exit("ERROR: number of rows ({0}) for {1} in VCF file doesn't match number of variantAlleles ({2}) in JSON file".format(
        len(vcfref[gene]), gene, len(jsondef["variantAlleles"])))

    allpattern = re.compile('ins|del')
    delpattern = re.compile('del')
    inspattern = re.compile('ins')
    for i,variant in enumerate(jsondef["variants"]):
        variant['_jsonidx'] = i

    for i,variant in enumerate(sorted(jsondef["variants"], key= lambda i: i['position'])):
       variant['_vcfidx'] = i
       if variant['rsid'] is not None and variant['rsid'] != vcfref[gene][i]["ID"]:
           sys.exit("ERROR: '%s' doesn't match any ID in VCF" % variant['rsid'])

       # look for ins/del and update from VCF reference
       if any(allpattern.match(allele) for allele in jsondef["variantAlleles"][variant['_jsonidx']]):
           variant["position"] = vcfref[gene][i]["POS"]
           vmap={}
           prebase=vcfref[gene][i]["REF"][0]
           # update indels with base before event (from VCF reference)
           for j,allele in enumerate(jsondef["variantAlleles"][variant['_jsonidx']]):
               if delpattern.match(allele):
                   vmap[allele] = prebase
               elif inspattern.match(allele):
                   vmap[allele] = prebase + allele[3:]
               else:
                   vmap[allele] = prebase+allele
               jsondef['variantAlleles'][variant['_jsonidx']][j] = vmap[allele]
           # for j, allele
           for na in jsondef['namedAlleles']:
               if na["alleles"][variant['_jsonidx']]:
                       na["alleles"][variant['_jsonidx']] = vmap[na["alleles"][variant['_jsonidx']]]
           # for na
       # isindel?
    # for i, variant
# updateVarAlleles()

def checkRefNamedAllele(jsondef, vcfref, refna):
    """
    Takes JSON dict, dict from parsed VCF and reference named allele
    Updates alleles in reference to match VCF where different
    """
    gene=jsondef['gene']
    for variant in jsondef["variants"]:
        if refna['alleles'][variant['_jsonidx']] != vcfref[gene][variant['_vcfidx']]['REF']:
            jsondef['variantAlleles'][variant['_jsonidx']] = [vcfref[gene][variant['_vcfidx']]["REF"] if v == refna['alleles'][variant['_jsonidx']] else v for v in jsondef['variantAlleles'][variant['_jsonidx']]]
            refna['alleles'][variant['_jsonidx']] = vcfref[gene][variant['_vcfidx']]['REF']
            refna['_maxdef'][variant['_jsonidx']] = vcfref[gene][variant['_vcfidx']]['REF']
    # for variant
# checkRefNamedAllele()


# define helper functions

symbolOffset = {}
def expandSymbol(symbol):
    """
    Takes a single basepair symbol or IUPAC code and yields each possible actual basepair.
    "A" -> "A" ; "M" -> "A","C" ; ...
    Also maintains a global state of the last yielded basepair for each code and begins
    each iteration from there, i.e. "M" -> "A","C" for the first caller, but "C","A" for
    the second. This avoids codes only ever expanding to the same (first) possibility
    if (for example) each caller consumes only the first returned value.
    """
    # resume from the previous iteration point so that combinatorial expansion of several multi-symbols doesn't start with A,A,A every time
    number = len(symbolBaselist[symbol])
    for n in range(number):
        offset = symbolOffset.get(symbol, 0)
        symbolOffset[symbol] = (offset + 1) % number
        yield symbolBaselist[symbol][offset % number]
# expandSymbol()

def expandString(string, limit=None):
    """
    Takes a string of one or more basepairs or IUPAC codes and returns an iterator over
    each possible actual basepair string after combinatorically expanding each symbol.
    "A" -> "A" ; "AC" -> "AC" ; "AM" -> "AA","AC" ; "MR" -> "AA","AG","CA","CG" ; ...
    """
    return itertools.islice(
        (''.join(expansion) for expansion in itertools.product(*(expandSymbol(symbol) for symbol in string))),
        0, limit, 1
    )
# expandString()

def expandAlleles(alleles, limit=None):
    """
    Takes an iterable of (index,string) pairs, for example defining a named gene allele,
    and returns an iterator over each possible actual genotype after combinatorically
    expanding each string (basepair(s) or IUPAC code(s)) at each index.
    ((1,"A"),(2,"C")) -> ((1,"A"),(2,"C")) ; ((1,"A"),(2,"M")) -> ((1,"A"),(2,"A")),((1,"A"),(2,"C")) ; ...
    """
    return itertools.islice(
        itertools.product(*(tuple((i, expansion) for expansion in expandString(string)) for i, string in alleles)),
        0, limit, 1
    )
# expandAlleles()

def symbolIntersect(string1, string2):
    """
    Takes two strings of basepairs or IUPAC codes and returns a string representing
    their pairwise intersection, or None if they don't intersect at one or more positions.
    ("AM","AR") -> "AA" ; ("AM","AK") -> None ; ...
    """
    try:
        intersection = ''.join(
            basesSymbol[frozenset(symbolBases[symbol1] & symbolBases[symbol2])]
            for symbol1, symbol2 in itertools.zip_longest(string1, string2, fillvalue='.')
        )
    except KeyError:
        return None
    return intersection
# symbolIntersect()

def percentEncode(string):
    """
    Takes a string and replaces all characters except numbers, letters, underscores, dashes and
    periods with "%" followed by their hexadecimal character code.
    "A B" -> "A%20B" ; ...
    """
    return re.sub(r"[^0-9A-Za-z_.-]", (lambda match: '%' + match[0].encode('utf-8').hex().upper()), string)
# percentEncode()

def joinMatches(matches):
    if len(matches) != 1:
        return "noCall%s" % (".".join(str(m) for m in sorted(matches)))
    else:
        return namedalleles[next(iter(matches))]['name'].replace('*', 's').replace('>', '-')
# joinMatches()


# load the named allele definitions
defPath = sys.argv[1]

print(f"Loading '{defPath}'...")
with open(defPath, 'r') as defFile:
    definition = json.load(defFile)
# with defFile

# update indel variants with information from positions VCF
vcfPath = sys.argv[2]
print(f"Loading '{vcfPath}'...")
vcfreference = parseVCF(vcfPath)
updateVarAlleles(definition,vcfreference)


if len(definition['variants']) != len(definition['variantAlleles']):
    sys.exit("ERROR: variants / variantAlleles length mismatch")
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
    # (eventually this should be explicit; for now we assume only the referent allele is defined for all positions,
    # except for total-gene deletion alleles which will be defined as 'delGene' for all positions)
    if not (None in namedallele['alleles']) and any((not a.startswith('del')) for a in namedallele['alleles']):
        if not refNamedallele:
            refNamedallele = namedallele
        elif len(definition['variants']) > 1:
            # if there's only one variant position, all named alleles must be defined for that one position
            # so we shouldn't warn about it and just have to trust that the correct "referent " allele is first
            print("  WARNING: two possible referent named alleles, '%s' and '%s'" %
                  (refNamedallele['name'], namedallele['name']))
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
        if any((s not in symbolBases) for s in a):
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
            if any((s not in symbolBases) for s in a):
                sys.exit("ERROR: named allele '%s' invalid allele #%d '%s'" % (refNamedallele['name'], i + 1, a))
            namedallele['_mindef'][i] = a
        # if del/ins/swap
    # for i,a in _mindef
# for namedallele

# generate tests
print("Generating test cases...")
# tests = collections.defaultdict(set)
for namedallele in namedalleles:
    namedallele['_tests'] = set()
    # generate tests for each namedallele with every position defined and referent used
    # for any position that is NULL in the defintion. Wobbles are expanded so that
    # a namedallele may have multiple tests
    namedallele['_maxdef'] = {
        i: namedallele['_mindef'].get(i, refNamedallele['_mindef'][i]) for i in range(len(definition['variants']))
    }
    for alleles in expandAlleles(namedallele['_maxdef'].items(), optNumAlleleExpansions):
        namedallele['_tests'].add(frozenset(alleles))
    # for alleles
# for namedallele

# check referent allele against VCF positions file and correct where needed for ambiguous
checkRefNamedAllele(definition, vcfreference, refNamedallele)

basepath = sys.argv[3]
if not os.path.exists(basepath):
    os.makedirs(basepath)
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

            outPath = os.path.join(basepath, percentEncode("%s_%s_%s_%s%s.vcf" % (
                definition['gene'],
                joinMatches(matches1),
                joinMatches(matches2),
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
                        (varalleles.index(a1) if (a1 and a1 != ".") else "."),
                        (varalleles.index(a2) if (a2 and a2 != ".") else "."),
                    ))
                    outFile.write("\t".join(row) + "\n")
                # for v,variant
            # with outFile
        # for alleles1, alleles2
    # for j,na2
# for i,na1

print(f"Done: {numFiles} files")

