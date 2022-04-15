#!/usr/bin/env python3
import sys
import os

input_path = sys.argv[1] # summary.txt file from StellarPGx
output_dir = sys.argv[2] # output directory

with open(input_path, "r") as infile:
	entries = infile.readlines() # read all the StellarPGx file lines into a list

for entry in entries: # for each sample in the StellarPGx file
	split = entry.split("\t") # split the lines into columns
	sample = split[0] # extract the sample name
	call = split[1].strip("[] \n") # extract the sample call and strip excess spaces and brackets

	if not os.path.exists(output_dir): # create the output directory if it doesn't exist yet
		os.makedir(output_dir)

	with open("%s/%s_cyp2d6.txt" % (output_dir,sample), "w") as outfile:
		# write the sample to a text file in the specified output directory
		outfile.write("CYP2D6\t%s" % call)