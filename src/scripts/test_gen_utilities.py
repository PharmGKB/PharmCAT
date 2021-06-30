#! /usr/bin/env python
__author__ = 'ScottDudek'

import collections
import re
import csv
import itertools

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

def joinMatches(matches, na):
    if len(matches) != 1:
        return "noCall%s" % (".".join(str(m) for m in sorted(matches)))
    else:
        return na[next(iter(matches))]['name'].replace('*', 's').replace('>', '-')
# joinMatches()
