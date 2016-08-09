#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import argparse

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description= \
	"""
	Uses tabix to grab overlapping fields (bed format, user defined) against a user-given \
	vcf file. The output is then parsed to give: \
	1. The number of variants overlapping the window \
	2. The number of variants with a recorded MAF \
	3. The average MAF of the above variants (highest MAF used)
	""", formatter_class=RawTextHelpFormatter)
parser.add_argument('-v','--vcf', required=True, \
	help = \
	'Give vcf file from ensembl (http://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/, \
	\"Homo_sapiens_incl_consequences.vcf.gz\") to search for variation across positions')
parser.add_argument('-b','--bed', required=True, \
	help =\
	'Coordinates to search for variation in bed format. Must be formatted the same as vcf')
parser.add_argument('-w','--window', default=50, type=int, \
	help = \
	'Window size to search up and down from position. Default is 50 (so 100bp covered total)')

def call_tabix(vcf, bed, window):
	# example tabix
	# tabix -p vcf the_file.vcf.gz 1:10000-11000
	bed_data = open(bed)
	for line in bed:
		sline = line.split()
		region = sline[0] + ':' + sline[1] + '-' + sline[2]
		tabix_query = 'tabix -p vcf ' + vcf + ' ' + region
		tabix = subprocess.check_output(tabix_query, shell=True).decode('utf-8')
		tabix = tabix.split('\n')
	return(tabix)

def parse_tabix(tabix_data):
	num_of_variant = len(tabix_data)
	regex=re.compile(r'.*MAF=0\.\d+', re.I)
	for line in tabix_data:
		sline = line.split('\t') 
		x=[m.group(0) for l in tabix[10].split(';') for m in [regex.search(l)] if m]
		print(x)
def main():
	args = parser.parse_args()
	vcf = args.vcf
	bed = args.bed
	window = args.window
