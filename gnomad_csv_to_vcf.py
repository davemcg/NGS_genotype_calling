#!/usr/bin/env python3

import argparse
import subprocess

parser = argparse.ArgumentParser(description= \
	"Takes unmodified ensembl human variation VCF, region to extract, and gnomad CSV output \
	to create a new VCF merging ensembl and gnomad data for a region (gene) \
	\
	Example: gnomad_csv_to_vcf.py Homo_sapiens.vcf.gz chr1:215,796,236-216,596,738 gnomad.csv | bgzip > Homo_sapiens.ensembl.gnomad.USH2A.vcf.gz")

parser.add_argument('ensembl_vcf', help = 'Unmodified ensembl vcf file from http://www.ensembl.org/info/data/ftp/index.html')
parser.add_argument('region', help = 'chr:start-stop to extract. chr and commas will be parsed out')
parser.add_argument('gnomad_csv', help = 'gnomad csv exported from gnomad.broadinstitute.org for the gene/region')

args = parser.parse_args()
ensembl = args.ensembl_vcf
region = args.region
gnomad = args.gnomad_csv

# write temp region file for tabix
chr = region.split('chr')[1].split(':')[0]
start = region.split(':')[1].split('-')[0].replace(',','')
stop = region.split(':')[1].split('-')[1].replace(',','')
output = chr + '\t' + start + '\t' + stop
with open('tabix_region.temp', 'w') as out:
	out.write(output)

# use tabix to grab region of interest from Ensembl VCF
tabix_call = 'tabix -hR tabix_region.temp ' + ensembl
region_info = subprocess.check_output(tabix_call, shell = True)

# pull header
header = [x.split('\t') for x in region_info.decode('utf-8').split('\n')[:-1] if x.split('\t')[0][0:2]=='##']
header.append(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'HET'])
# pull dbSNP ids
ids=[x.split('\t')[2] for x in region_info.decode('utf-8').split('\n')[:-1] if x.split('\t')[0][0]!='#']
# pull variants
variation = [x.split('\t')[0:4] + ['100'] + x.split('\t')[5:] + ['GT','0/1'] \
			for x in region_info.decode('utf-8').split('\n')[:-1] \
			if x.split('\t')[0][0]!='#']

# append gnomad variation to ensembl, skipping ones with matching rs IDs
for line in open(gnomad):
	line = line.replace('\"','').split(',')
	if line[2] not in ids and line[0] != 'Chrom':
		gnomad_line = [\
			line[0], \
			line[1], \
			line[2], \
			line[3], \
			line[4], \
			"100", "PASS", ".", "GT", "0/1"]
		variation.append(gnomad_line)


# sort by start position
variation.sort(key=lambda x: int(x[1]))
# print out new file to stdout
for l in header:
	print('\t'.join(l))
for l in variation:
	print('\t'.join(l))
