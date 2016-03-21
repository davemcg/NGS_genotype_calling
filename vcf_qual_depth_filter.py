#!/usr/local/Anaconda/envs/py3.4.3/bin/python

"""
Filters vcf filter for qual and depth, tagging the filter field
"""

import sys

for line in sys.stdin:
	if line[0]=='#':
		print(line[:-1])
	elif float(line.split()[5]) > 20 and int(line.split()[9].split(':')[1]) > 4:
		line = line.split()
		line[6] = "PASS"
		print('\t'.join(line))
	else:
		line = line.split()
		line[5] = "FAIL"
	#	print('\t'.join(line))
		
