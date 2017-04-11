#!/bin/env python2
from __future__ import print_function

# Multi-step script that:
# Takes HGVS in GeneName:c. notation (e.g. CHD7:c.3404C>A)
# Converts gene name to RefSeq NM in two steps:
# Grabs the highest appris rating from the GencodeGenes GRCh37 Basic GTF
#  ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/GRCh37_mapping/gencode.v26lift37.basic.annotation.gtf.gz
# Uses the Ensembl transcript name (e.g. ENST00000423902.6) 
# Then searches for the accompanying RefSeq name in the GencodeGenes RefSeq metadata file
#  ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/GRCh37_mapping/gencode.v26lift37.metadata.RefSeq.gz
# The validate the new HGVS NM.c notation with biocommons (invitae) HGVS
#  https://github.com/biocommons/hgvs
# and convert to GRCh37 and 38 g. and output the name
# also indicate failure at any stage (variant will then require hand conversion)

import sys
import re

hgvs = open(sys.argv[1], 'r')
gtf = open(sys.argv[2], 'r')
refseq = open(sys.argv[3], 'r')

hgvs_gene_gtf = {}

# Collapse tx info to the gene name
for line in gtf:
	if line[0] == '#':
		continue
	line = line[:-1]
	if 'transcript' in line.split()[2] and 'gene_name' in line:
		gene_index = re.split('[;\s+]',line).index('gene_name')
		gene = re.split('[;\s+]',line)[gene_index+1].replace('\"','')
		#tx = line.split()[11].split(';')[0][1:-1]
		if gene not in hgvs_gene_gtf:
			hgvs_gene_gtf[gene] = line
		else:
			line_old = hgvs_gene_gtf[gene]
			line_new = line_old + '___' + line
			hgvs_gene_gtf[gene] = line_new

appris_ranking = \
	['appris_principal_1',
	'appris_principal_2',
	'appris_principal_3',
	'appris_principal_4',
	'appris_principal_5',
	'appris_alternative_1',
	'appris_alternative_2',
	'appris_principal',
	'appris_candidate',
	'appris_candidate_ccds',
	'appris_candidate_highest_score',
	'appris_candidate_longest',
	'appris_candidate_longest_ccds',
	'appris_candidate_longest_seq']

gene_best_tx = {}
for gene in hgvs:
	gene = gene[:-1]
	if gene in hgvs_gene_gtf.keys():
		for rank in appris_ranking:
			if rank in hgvs_gene_gtf[gene]:
				tx_index = re.split('[;\s+]',hgvs_gene_gtf[gene]).index('transcript_id')
				gene_best_tx[gene] = re.split('[;\s+]',hgvs_gene_gtf[gene])[tx_index+1].replace('\"','').split('.')[0]
				break
			else:
				gene_best_tx[gene] = 'No appris'
	else:
		gene_best_tx[gene] = 'Not in Gencode GTF'

tx_refseq = {}
for line in refseq:
	tx_refseq[line.split('.')[0]] = line.split()[1]

for k,v in gene_best_tx.items():
	
