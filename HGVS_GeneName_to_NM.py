#!/bin/env python2
from __future__ import print_function

# Multi-step script that:
# Takes HGVS in GeneName:c. notation (e.g. CHD7:c.3404C>A)
# Converts gene name to RefSeq NM in two steps:
# First creates a dict with Ensembl TX - NCBI TX | key - value
#  ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/GRCh37_mapping/gencode.v26lift37.metadata.RefSeq.gz
# Grabs the highest appris rating from the GencodeGenes GRCh37 Basic GTF, skipping those which don't match to NCBI
#  ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/GRCh37_mapping/gencode.v26lift37.basic.annotation.gtf.gz
# Uses the Ensembl transcript name (e.g. ENST00000423902.6) then map over
# If the two step approach above fails to find a gene, then try the full RefSeq tx file (from UCSC) and grab the longest tx
#  http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqCurated.txt.gz
# The validate the new HGVS NM.c notation with biocommons (invitae) HGVS
#  https://github.com/biocommons/hgvs
# and convert to GRCh37 and 38 g. and output the name
# also indicate failure at any stage (variant will then require hand conversion)

import sys
import re
import hgvs.validator
import hgvs.exceptions
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.assemblymapper

hgvs_file= open(sys.argv[1], 'r')
gtf = open(sys.argv[2], 'r')
refseq_gencode = open(sys.argv[3], 'r')
refseq_ucsc = open(sys.argv[3], 'r')

# ensembl tx as key, refseq as value
tx_gencode_refseq = {}
for line in refseq_gencode:
	tx_gencode_refseq[line.split('.')[0]] = line.split()[1]

hgvs_gene_gtf = {}
# Collapse tx info to the gene name, skipping tx who aren't a key in tx_gencode_refseq
for line in gtf:
	if line[0] == '#':
		continue
	line = line[:-1]
	if 'transcript' in line.split()[2] and 'gene_name' in line:
		gene_index = re.split('[;\s+]',line).index('gene_name')
		gene = re.split('[;\s+]',line)[gene_index+1].replace('\"','')
		tx_index = re.split('[;\s+]',line).index('transcript_id')
		tx = re.split('[;\s+]',line)[tx_index+1].replace('\"','').split('.')[0]
		if tx not in tx_gencode_refseq:
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

# Roll through hgvs_filefile line by line
# extract out gene name, match with best
# guess for RefSeq NM transcript
gene_best_tx = [] 
for line in hgvs_file:
	line = line[:-1]
	gene = line.split(':')[0]
	if gene in hgvs_gene_gtf:
		for rank in appris_ranking:
			if rank in hgvs_gene_gtf[gene]:
				info = hgvs_gene_gtf[gene].split('___')
				rank_truth = [rank in s for s in info]
				rank_index = [i for i, x in enumerate(rank_truth) if x][0]
				tx_index = re.split('[;\s+]',info[rank_index]).index('transcript_id')
				ensembl_tx = re.split('[;\s+]',info[rank_index])[tx_index+1].replace('\"','').split('.')[0]
				out = [line, gene, tx_gencode_refseq[ensembl_tx]]
				gene_best_tx.append(out)
				break
			# case where no appris info, but there's still a refseq tx, just take the first one
			# else: 
			#	info = hgvs_gene_gtf[line].split('___')[0]
			#	tx_index = re.split('[;\s+]',info).index('transcript_id')
			#	out = [line, gene, re.split('[;\s+]',info)[tx_index+1].replace('\"','').split('.')[0]]
			#	gene_best_tx.append(out)
	else:
		out = [line, gene, 'Not in Gencode GTF/RefSeq']
		gene_best_tx.append(out)

# Now run through the list again and try to fill in the missing with the longest RefSeq transcript
# But first, build the dict, using the longest tx
tx_refseq_name = {}
for line in refseq_ucsc:
	line = line[:-1]
	refseq_tx = line.split()[1]
	gene = line.split()[12]
	start = int(line.split()[4])
	stop = int(line.split()[5])
	size = abs(stop-start)
	if gene not in tx_refseq_name:
		tx_refseq_name[gene] = [refseq_tx, size]
	else:
		old_size = tx_refseq_name[gene][1]
		if size > old_size:
			tx_refseq_name[gene] = [refseq_tx, size]
# ok, now run through gene_best_tx and try to fill in some more missing
for line in gene_best_tx:
	if line[2]=='Not in Gencode GTF/RefSeq':
		if line[1] in tx_refseq_name:
			line[2] = tx_refseq_name[line[1]][0]

# setup biommons hgvs tooling
hp = hgvs.parser.Parser()
hdp = hgvs.dataproviders.uta.connect()
vr = hgvs.validator.Validator(hdp=hdp)
vm37 = hgvs.assemblymapper.AssemblyMapper(
    hdp, assembly_name='GRCh37', alt_aln_method='splign')
vm38 = hgvs.assemblymapper.AssemblyMapper(
    hdp, assembly_name='GRCh38', alt_aln_method='splign')
# Run through each line, build new HGVS and validate
converted_hgvs = []
for variant in gene_best_tx:
	original_hgvs = variant[0] 
	tx = variant[2]
	hgvs_right = variant[0].split(':')[1]
	new_hgvs = tx + ':' + hgvs_right
	if tx != 'Not in Gencode GTF/RefSeq':
		try:
			variant = hp.parse_hgvs_variant(new_hgvs)
			vr.validate(variant)
			# on successful validation, try to convert
			try:
				p_dot = str(vm37.c_to_p(variant))
				g37_dot = str(vm37.c_to_g(variant))
				g38_dot = str(vm38.c_to_g(variant))
				out = [original_hgvs, new_hgvs, p_dot, g37_dot, g38_dot, 'Success']
			except:
				out = [original_hgvs, new_hgvs, 'null', 'null', 'null', 'Liftover failure']
			converted_hgvs.append(out)
		except hgvs.exceptions.HGVSError as e:
			out = [original_hgvs, new_hgvs, 'null', 'null', 'null', e]
			converted_hgvs.append(out)
	else:
		converted_hgvs.append([original_hgvs, 'null', 'null', 'null', 'null', 'Gene not in Gencode GTF or RefSeq'])






