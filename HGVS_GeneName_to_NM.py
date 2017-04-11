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


# ensembl tx as key, refseq as value
tx_refseq = {}
for line in refseq:
	tx_refseq[line.split('.')[0]] = line.split()[1]

hgvs_gene_gtf = {}
# Collapse tx info to the gene name, skipping tx who aren't a key in tx_refseq
for line in gtf:
	if line[0] == '#':
		continue
	line = line[:-1]
	if 'transcript' in line.split()[2] and 'gene_name' in line:
		gene_index = re.split('[;\s+]',line).index('gene_name')
		gene = re.split('[;\s+]',line)[gene_index+1].replace('\"','')
		tx_index = re.split('[;\s+]',line).index('transcript_id')
		tx = re.split('[;\s+]',line)[tx_index+1].replace('\"','').split('.')[0]
		if tx not in tx_refseq.keys():
			continue
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

gene_best_tx = [] 
for line in hgvs:
	line = line[:-1]
	gene = line.split()[0]
	if gene in hgvs_gene_gtf.keys():
		for rank in appris_ranking:
			if rank in hgvs_gene_gtf[gene]:
				info = hgvs_gene_gtf[gene].split('___')
				rank_truth = [rank in s for s in info]
				rank_index = [i for i, x in enumerate(rank_truth) if x][0]
				tx_index = re.split('[;\s+]',info[rank_index]).index('transcript_id')
				out = [line, gene, re.split('[;\s+]',info[rank_index])[tx_index+1].replace('\"','').split('.')[0]]
				gene_best_tx.append(out)
				break
			# case where no appris info, but there's still a refseq tx, just take the first one
			else: 
				info = hgvs_gene_gtf[line].split('___')[0]
				tx_index = re.split('[;\s+]',info).index('transcript_id')
				out = [line, gene, re.split('[;\s+]',info)[tx_index+1].replace('\"','').split('.')[0]]
				gene_best_tx.append(out)
	else:
		out = [line, gene, 'Not in Gencode GTF/RefSeq']
		gene_best_tx.append(out)


for k,v in gene_best_tx.iteritems():
	if v in tx_refseq.keys():
		print(k,v,tx_refseq[v])
	else:
		print(k,v,'Not in refseq?')
