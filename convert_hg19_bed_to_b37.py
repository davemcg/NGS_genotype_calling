#!/usr/local/Anaconda/envs/py3.4.3/bin/python


"""
Converts hg19 bed files to b37

e.g. 
chr1 100 110 will be converted to 1 100 110 
chrM 100 110 to MT 100 110
"""

import sys

counter = 0
for line in open(sys.argv[1]):
	line = line.split()
	if line[0][0:3] != "chr":
		counter += 1
	elif line[0] == "chrM":
		line[0] = "MT"
	elif counter > 1:
		print("ERROR, multiple non chr lines")
		break
	else:
		line[0] = line[0][3:]
	print('\t'.join(line))
