#!/usr/bin/env python3

"""
Usage: test_gen.py <definition.json> <positions.vcf> <vcfOutputDir>
Author: Alex Frase

This tool generates test cases for PharmCAT by reading a single gene definition JSON file
and the PharmCAT positions file. It then writes a set of synthetic VCF files with their
expected haplotype calls for a variety of scenarios.

Tests include:
* Homozygous (two copies of the same allele) complete (all basepair positions included and
  non-missing) and incomplete (only positions relevant to the allele included, all others
  missing) cases for each allele; *1/*1, *2/*2, ...
* Heterozygous complete and incomplete cases for each pair of alleles; *1/*2, *1/*3, *2/*3, ...
* A case with one referent allele and one uncallable allele due to a non-referent base
  that is not part of any minor allele definition
* Cases with one referent allele and one uncallable allele due to missing data that
  causes ambiguity between two or more possible minor alleles; several tests are generated
  for each cluster of minor alleles with overlapping definitions: some which are ambiguous
  only between two minor alleles in the cluster (one test for each such pair), and some
  which are ambiguous between all but one minor allele in the cluster (one test which
  rules out each allele)

Tests are grouped by expected haplotype call, such that each VCF file covers a single
call (reflected in the filename and the "##PharmCATnamedAlleles" meta information header)
but may contain multiple "samples", one for each generated test case that yielded the same
final call.
"""

# TODO: once available, use definition's specified referent named allele rather than just the first fully-defined named allele
# TODO: add maxdef test cases for cluster ambiguity

# ask: definition specifies genomeBuild "b38", vcf seems to expect "hg38"
# ask: does the VCF need 'source' meta info? old test files have it

# ASK: VCFv4.3 ALT is defined as
#          "base Strings made up of the bases A,C,G,T,N (case insensitive)
#          or the '*'  symbol (allele missing due to overlapping deletion)
#          or a MISSING value ‘.’ (no variant)"
#      some gene definitions have deletions with a length that would run
#      over subsequent variants in the definition; how is this handled?
#      will PharmCAT parse the '*' symbol correctly in that context?

import collections
import csv
import datetime
import itertools
import json
import math
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
optNumAlleleExpansions = 2

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
    Returns a collection with gene names (from 'INFO' column) as key
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
        yield symbolBaselist[symbol][(n + offset) % number]
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
unalleles = set()
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
        # initialize a pool of "unknown alleles" starting with all possible non-referent alleles
        if len(a) == 1:
            unalleles.update((i, b) for b in ('A', 'C', 'G', 'T') if b not in symbolBases[a])
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
    # update the pool of unknown alleles, removing any that are part of this definition
    unalleles.difference_update(*(expandAlleles(namedallele['_mindef'].items())))
# for namedallele
unalleles = tuple(sorted(unalleles))
unalleleOffset = 0
print(f"done: {len(unalleles)} possible unknown alleles\n")

# generate basic tests
print("Generating test cases...")
tests = collections.defaultdict(set)
first = True
for namedallele in namedalleles:
    # generate tests for both the minimum and maximum definition of the named allele.
    # the minimum definition is the way it's notated in the original JSON, with only
    # the non-referent positions defined for the named allele and all others left null,
    # implying that the named allele could be identified by variants at those positions
    # alone; the maximum definition fills in all those missing positions from the referent
    # named allele, leaving no positions missing.
    namedallele['_maxdef'] = {
        i: namedallele['_mindef'].get(i, refNamedallele['_mindef'][i]) for i in range(len(definition['variants']))
    }
    if namedallele['_mindef'] != namedallele['_maxdef']:
        for alleles in expandAlleles(namedallele['_mindef'].items(), optNumAlleleExpansions):
            test = frozenset(alleles)
            tests[test].add(namedallele['name'])
            if first:
                # for the first non-referent allele, also generate a test with an unexpected allele
                # (pull one from the pool which is at a position not already part of this allele's mindef)
                first = False
                indecies = set(i for i, a in test)
                o = unalleleOffset
                while unalleles[o][0] in indecies:
                    o = (o + 1) % len(unalleles)
                # no-op just to populate the defaultdict key
                tests[frozenset(test.union(unalleles[o:(o + 1)]))].discard(None)
                unalleleOffset = (o + 1) % len(unalleles)
            # if first
        # for alleles
    # if min!=max
    for alleles in expandAlleles(namedallele['_maxdef'].items(), optNumAlleleExpansions):
        tests[frozenset(alleles)].add(namedallele['name'])
# for namedallele

# identify clusters of overlapping alleles
clusters = list()
naCluster = dict()
n = m = 0
for na1, namedallele1 in enumerate(namedalleles):
    # for each named allele that's not yet part of a cluster, we start with the variants
    # in that named allele's mindef and then iteratively expand the cluster to include
    # all other named alleles that share mindef variants with any that are in the cluster
    # (here, "variants" means both the position and the allele
    if (namedallele1 is refNamedallele) or (na1 in naCluster):
        continue
    cluster = {na1}
    naCluster[na1] = len(clusters)
    clusters.append(cluster)
    union = set()
    for alleles in expandAlleles(namedallele1['_mindef'].items()):
        union.update(alleles)
    expand = True
    while expand:
        expand = False
        for na2, namedallele2 in enumerate(namedalleles):
            if (namedallele2 is refNamedallele) or (na2 in naCluster):
                continue
            if any(any(ia in union for ia in alleles) for alleles in expandAlleles(namedallele2['_mindef'].items())):
                expand = True
                cluster.add(na2)
                naCluster[na2] = naCluster[na1]
                for alleles in expandAlleles(namedallele2['_mindef'].items()):
                    union.update(alleles)
            # if overlap
        # for na2,namedallele2
    # while expand
    # just a little fun, tallying up the number of possible subsets of this cluster in
    # terms of named alleles and of specific variants. I added this when debugging an
    # early draft that tried to generate tests for *all* possible subsets of each cluster
    # and was running out of memory; turns out there are a few BIG clusters (see below)
    for i in range(len(cluster) + 1):
        n += (math.factorial(len(cluster)) / (math.factorial(i) * math.factorial(len(cluster) - i)))
    for i in range(len(union) + 1):
        m += (math.factorial(len(union)) / (math.factorial(i) * math.factorial(len(union) - i)))
# for na1,namedallele1
# sys.exit("%d / %d combinations\n" % (n,m))

# generate ambiguity tests for each cluster
for cluster in clusters:
    # for some subsets of each cluster, test for the alleles they all have in common
    # plus the referent alleles that distinguish them from the rest of the cluster.
    # (limit to subsets of 1,2,n-1,n; CYP2D6 has a cluster of 91 named alleles with
    # 2,475,880,078,570,759,450,286,620,672 or ~2.5 octillion possible subsets)
    indecies = set.union(*(namedalleles[na]['_indecies'] for na in cluster))
    sizes = {1, min(2, len(cluster)), max(1, len(cluster) - 1), len(cluster)}
    for size in sizes:
        for subcluster in itertools.combinations(cluster, size):
            # build a _mindef as the intersection of all subcluster members
            subdef = dict(namedalleles[subcluster[0]]['_mindef'])
            for na in itertools.islice(subcluster, 1, None, 1):
                for i in tuple(subdef.keys()):  # copy so we can delete
                    subdef[i] = symbolIntersect(subdef[i], namedalleles[na]['_mindef'].get(i, '.'))
                    if not subdef[i]:
                        del subdef[i]
                # for i
            # for na
            if subdef:
                # add referent alleles for variants present in the cluster but not the subcluster
                referent = dict(next(expandAlleles(refNamedallele['_maxdef'].items(), 1)))
                subdef.update((i, referent[i])
                              for i in indecies.difference(*(namedalleles[na]['_indecies'] for na in subcluster)))
                # generate test case(s) for this subcluster mindef
                for alleles in expandAlleles(subdef.items(), optNumAlleleExpansions):
                    tests[frozenset(alleles)].update(namedalleles[na]['name'] for na in subcluster)
                # add referent alleles for all remaining positions and generate test case(s) for the maxdef
                subdef.update((i, referent[i]) for i in range(len(definition['variants'])) if (i not in indecies))
                for alleles in expandAlleles(subdef.items(), optNumAlleleExpansions):
                    tests[frozenset(alleles)].update(namedalleles[na]['name'] for na in subcluster)
            # if intersection
        # for subcluster
    # for size
# for cluster
print(f"done: {len(tests)} tests\n")

# analyze test cases
for test in tests.keys():
    alleles = dict(test)

    # identify candidate named alleles
    candidates = collections.OrderedDict()
    for na, namedallele in enumerate(namedalleles):
        match = collections.Counter(
            bool(symbolIntersect(alleles[i], namedallele['_maxdef'][i])) for i in alleles.keys())
        if match[True] and not match[False]:
            candidates[na] = match[True]
    # for na,namedallele

    # drop superseded candidates
    if candidates:
        best = max(candidates.values())
        for na in list(candidates):  # copy keys to a list so we can delete from the original as we go
            if candidates[na] < best:
                del candidates[na]
        # for candidate
    elif tests[test]:
        print("WARNING: test '%s' expected %s, matched none" % (alleles, tests[test]))
    # if candidates

    # store result
    tests[test] = (tuple(candidates.keys()), -len(test))
# for test

# collate tests by expected calls
alltests = list()
matchesTests = collections.defaultdict(list)
for test, testmatches in tests.items():
    alltests.append(test)
    matchesTests[testmatches[0]].append(test)
# for test,matches
alltests.sort(key=tests.get)
for matches in matchesTests.keys():
    matchesTests[matches].sort(key=tests.get)
# for matches

# output results
if False:
    # old CSV format, for initial manual inspection
    outFile = csv.writer(sys.stdout)
    outFile.writerow(variant['position'] for variant in definition['variants'])
    outFile.writerow(variant['rsid'] for variant in definition['variants'])
    outFile.writerow(variant['chromosomeHgvsName'] for variant in definition['variants'])
    outFile.writerow(variant['geneHgvsName'] for variant in definition['variants'])
    for test in alltests:
        alleles = dict(test)
        row = list(alleles.get(i, '') for i in range(len(definition['variants'])))
        row.append('; '.join(namedalleles[na]['name'] for na in tests[test][0]))
        outFile.writerow(row)
    # for test
    sys.exit(0)
# if CSV

basepath = sys.argv[3]
if not os.path.exists(basepath):
    os.makedirs(basepath)
print("Writing files...")
numFiles = 0
for matches1, matches2 in itertools.combinations_with_replacement(sorted(matchesTests.keys()), 2):
    # always put uncallable haplotypes second
    if len(matches1) != 1:
        matches1, matches2 = matches2, matches1
    # only generate uncallables paired with *1/...
    if (len(matches1) != 1) or (len(matches2) != 1 and matches1 != (refNamedallele['_num'],)):
        continue
    outPath = os.path.join(basepath, percentEncode("%s_%s_%s.vcf" % (
        definition['gene'],
        joinMatches(matches1),
        joinMatches(matches2),
    )))
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

        # write variants header
        match1tests = matchesTests[matches1]
        match2tests = matchesTests[matches2]
        test1alleles = list(dict(test) for test in match1tests)
        test2alleles = list(dict(test) for test in match2tests)
        row = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        row.extend(("TEST%d" % (t + 1,)) for t in range(len(match1tests) * len(match2tests)))
        outFile.write("\t".join(row) + "\n")

        # write variants data
        for v, variant in enumerate(definition['variants']):
            refAllele = refNamedallele['_maxdef'].get(v)
            if not refAllele or refAllele == ".":
                continue
            varalleles = [refNamedallele['_maxdef'].get(v)]
            # collate alt alleles used by these tests
            for t, alleles in itertools.chain(enumerate(test1alleles), enumerate(test2alleles)):
                alt = (None if ((v not in alleles) or alleles[v] == ".") else alleles[v])
                if alt and (alt not in varalleles):
                    varalleles.append(alt)
            # for t, alleles
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
            for alleles1, alleles2 in itertools.product(test1alleles, test2alleles):
                a1 = alleles1.get(v)
                a2 = alleles2.get(v)
                row.append("%s/%s" % (
                    (varalleles.index(a1) if (a1 and a1 != ".") else "."),
                    (varalleles.index(a2) if (a2 and a2 != ".") else "."),
                ))
            # for t,alleles
            outFile.write("\t".join(row) + "\n")
        # for variant
    # with outFile
# for matches1,matches2
print(f"Done: {numFiles} files")
