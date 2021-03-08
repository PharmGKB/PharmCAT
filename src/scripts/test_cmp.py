#!/usr/bin/env python3

"""
Usage: test_cmp.py <olddir> <newdir>
Author: Alex Frase

This tool compares old and new test cases for PharmCAT by scanning all samples of all
test VCF files to determine how many of the old samples (tests) are represented in the
new tests.

The comparison is inexact; because the named allele definitions have changed, old and new
tests for the same gene may have a few different variant positions. Only the intersection
of positions present in both files are compared; this isn't ideal since every included
variant is necessary to distinguish at least one named allele, but suffices to get an
initial sense of coverage.

A futue option would be to add a variant position mapping table for each gene, so that
variants which changed positions can still be matched.
"""

import collections
import glob
import os
import re
import sys


if len(sys.argv) < 3:
	exit("error: two paths required")
oldpath = sys.argv[1]
newpath = sys.argv[2]

def getVCFSampleAlleles(vcfPath):
	"""
	Takes a path to a VCF file and returns an ordered dict of dicts:
	for each sample, for each variant (chm,pos) tuple, a sorted tuple
	of alleles as basepair strings.
	"""
	sampleAlleles = collections.OrderedDict()
	with open(vcfPath,'rU') as vcfFile:
		header = None
		for line in vcfFile:
			if line.startswith("#CHROM"):
				if header:
					exit("repeat header in %s" % (vcfName,))
				header = line.rstrip("\r\n").split() # these SHOULD be tab-delimited just like the data rows, but aren't always :/
				if len(header) < 9:
					exit("invalid header in %s: %s" % (vcfPath,line))
				# sampleAlleles is an OrderedDict, so just by initializing it in
				# column order now, we can iterate over it in that same order later
				for sample in header[9:]:
					if sample in sampleAlleles:
						exit("duplicate sample in %s: %s" % (vcfPath,sample))
					sampleAlleles[sample] = dict()
				#for sample
			elif header:
				row = line.rstrip("\r\n").split("\t")
				if len(row) != len(header):
					exit("invalid row in %s: %s" % (vcfPath,line))
				chmpos = (row[0], int(row[1]))
				# combine the ref and alt allele(s) into one list to look up alleles by index
				ref = row[3]
				alts = (list() if (row[4] == ".") else row[4].split(','))
				alleles = [ref] + alts
				for s,sample in enumerate(sampleAlleles.keys()):
					# store the alleles in sorted order even for phased data so that comparisons work
					# even if the ref/alt are swapped between files
					sampleAlleles[sample][chmpos] = tuple(sorted(("." if (a == ".") else alleles[int(a)]) for a in re.split(r'/|\|', row[9+s])))
				#for s,sample
			#if header
		#for line
	#with vcfFile
	return sampleAlleles
#getVCFSampleAlleles()

# load all the old tests into memory
vcfSampleAlleles = dict()
nf0,nt0 = 0,0
for vcfPath0 in glob.iglob(os.path.join(oldpath, "*.vcf")):
	vcfSampleAlleles[vcfPath0] = getVCFSampleAlleles(vcfPath0)
	nf0 += 1
	nt0 += len(vcfSampleAlleles[vcfPath0])
#for vcfPath
sys.stderr.write("%s contains %d tests in %d files\n" % (oldpath,nt0,nf0))

# scan all the new tests and check for coverage of old tests
nf1,nt1,nm = 0,0,0
for vcfPath1 in glob.iglob(os.path.join(newpath, "*.vcf")):
	nf1 += 1
	for sample1,alleles1 in getVCFSampleAlleles(vcfPath1).items():
		nt1 += 1
		# loop over all old files and samples looking for matches;
		# this makes it an n^2 operation but that's fast enough, and would be tricky to
		# improve because of the inexact matching between covered variants of each gene
		for vcfPath0,sampleAlleles0 in vcfSampleAlleles.items():
			matches = list()
			for sample0,alleles0 in sampleAlleles0.items():
				if alleles0:
					# identify the common variant positions (if any) and match only on those
					variants = alleles0.keys() & alleles1.keys()
					if variants and all(alleles0[v] == alleles1[v] for v in variants):
						matches.append(sample0)
				#if not yet matched
			#for sample0,alleles0
			nm += len(matches)
			for match in matches:
				vcfSampleAlleles[vcfPath0][match] = None
		#for vcfPath0,sampleAlleles0
	#for sample1,alleles1
#for vcfPath1
sys.stderr.write("%s contains %d tests in %d files\n" % (newpath,nt1,nf1))

# report coverage
sys.stderr.write("%d (%d%%) old tests matched, %d (%d%%) unmatched\n" % (nm,(100*nm/nt0),(nt0-nm),(100*(nt0-nm)/nt0)))

# TODO: report specific old tests not yet covered, to investigate what they test for and how to generate matching tests
