#!/usr/local/Anaconda/envs_app/hgvs/1.0.0/bin/python 
from __future__ import print_function

import sys
import re
import hgvs.validator
import hgvs.exceptions
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.assemblymapper
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description = \
"""
This script takes HGVS in GeneName:c. notation (e.g. CHD7:c.3404C>A), with
one variant per line, and outputs a tab separated text file with the original
HGVS notation, a proper HGVS with the RefSeq transcript name, then several 
different conversion and liftovers to p. and g. For example if given CHD7:c.3404C>A,
the following is returned:
Original_HGVS   Validated_HGVS_c.       HGVS_p. HGVS_g._hg19    HGVS_g._hg38    Status
CHD7:c.3404C>A  NM_017780.3:c.3404C>A   NP_060250.2:p.(Thr1135Asn)      NC_000008.10:g.61741247C>A      NC_000008.11:g.60828688C>A      Success

The Status field will indicate failure (and the reason).

The four most common errors will relate to:
1. Failure to find a transcript for a gene name. 
2. Intronic variants cannot be validated (fuzzyness on where exactly they lie)
3. HGVS reference does not agree with genome reference (either because the HGVS is wrong or the wrong transcript was picked)
4. Liftover fails (p. cannot always be generated, occasionally cannot match hg19 to hg38)

Run this script with -h to get a listing of the types of input files required (and some optional ones)
""",
formatter_class=RawTextHelpFormatter)

parser.add_argument('hgvs_file', help = 'line separated HGVS file in GeneName:c. format for conversion and liftover\n\n')
parser.add_argument('gtf', help = 'GTF file from ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/GRCh37_mapping/gencode.v26lift37.basic.annotation.gtf.gz\n\n')
parser.add_argument('refseq_gencode', help = 'Ensembl RefSeq transcript match file from ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_26/GRCh37_mapping/gencode.v26lift37.metadata.RefSeq.gz\n\n')
parser.add_argument('refseq_ucsc', help = 'Curated RefSeq transcript info from UCSC http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqCurated.txt.gz\n\n')
parser.add_argument('--manual_conversion', help = 'User can give conversion table in with the GeneName and the accompanying RefSeq transcript. One set (space separated) per line\n\n')
parser.add_argument('--exhaustive', default = 'off', help = 'Default is \'off\'. Setting to \'on\' will query each variant with every available RefSeq transcript\n\n')
parser.add_argument('output', help = 'Output file name. Output will be tsv\n\n')

args = parser.parse_args()
hgvs_file = open(args.hgvs_file, 'r')
gtf = open(args.gtf, 'r')
refseq_gencode = open(args.refseq_gencode, 'r')
refseq_ucsc = open(args.refseq_ucsc, 'r')
if args.manual_conversion is not None:
	manual_conversion = open(args.manual_conversion, 'r')
exhaustive_flag = args.exhaustive
output_file = open(args.output, 'w')

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
	line = line.strip()
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
	else:
		out = [line, gene, 'Not in Gencode GTF/RefSeq']
		gene_best_tx.append(out)

# But first, build the dict, using the longest tx as the value with gene as the key
# also build a second dictionary for the --exhaustive search with all tx as the values
tx_refseq_name = {}
#alltx_refseq_name = {}
for line in refseq_ucsc:
	line = line[:-1]
	refseq_tx = line.split()[1]
	gene = line.split()[12]
	start = int(line.split()[4])
	stop = int(line.split()[5])
	size = abs(stop-start)
	if gene not in tx_refseq_name:
		tx_refseq_name[gene] = [refseq_tx, size]
#		alltx_refseq_name[gene] = [refseq_tx]
	else:
#		orig_tx = alltx_refseq_name[gene]
#		new_tx = orig_tx.append(refseq_tx)
#		alltx_refseq_name[gene] = new_tx

		old_size = tx_refseq_name[gene][1]
		if size > old_size:
			tx_refseq_name[gene] = [refseq_tx, size]
# ok, now run through gene_best_tx and try to fill in some more missing
for line in gene_best_tx:
	if line[2]=='Not in Gencode GTF/RefSeq':
		if line[1] in tx_refseq_name:
			line[2] = tx_refseq_name[line[1]][0]

# now use user given conversion table to override choices made in assignment RefSeq tx to gene name
if args.manual_conversion is not None:
	# build dict
	manual_con = {}
	for line in manual_conversion:
		if line.split('\t') != 2:
			continue
		manual_con[line.split()[0]] = line.split()[1]
	for line in gene_best_tx:
		if line[1] in manual_con:
			line[2] = manual_con[line[1]]

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
	try:
		hgvs_right = variant[0].split(':')[1]
	except:
		converted_hgvs.append([original_hgvs, 'null', 'null', 'null', 'null', 'Missing c dot'])
	try:
		hgvs_right.decode('ascii')
	except:
		converted_hgvs.append([original_hgvs, 'null', 'null', 'null', 'null', 'Improper formatting. Ascii in input'])
		continue
	new_hgvs = tx + ':' + hgvs_right
	if tx != 'Not in Gencode GTF/RefSeq':
		try:
			variant = hp.parse_hgvs_variant(new_hgvs)
			vr.validate(variant)
			# on successful validation, try to convert
			error = 'Success'
			try:
				p_dot = str(vm37.c_to_p(variant))
			except:
				error = "Cannot make p.: " + str(sys.exc_info()[0])
				p_dot = 'null'
			try:
				g37_dot = str(vm37.c_to_g(variant))
			except:
				error = "Cannot make GRCh37 g.: " + str(sys.exc_info()[0])
				g37_dot = 'null'
			try:
				g38_dot = str(vm38.c_to_g(variant))
			except:
				error = "Cannot make GRCh38 g.: " + str(sys.exc_info()[0])
				g38_dot = 'null'
			out = [original_hgvs, new_hgvs, p_dot, g37_dot, g38_dot, error]
			converted_hgvs.append(out)
		except hgvs.exceptions.HGVSError as e:
			out = [original_hgvs, new_hgvs, 'null', 'null', 'null', str(e)]
			converted_hgvs.append(out)
	else:
		converted_hgvs.append([original_hgvs, 'null', 'null', 'null', 'null', 'Gene not in Gencode GTF or RefSeq'])

# output!
output_file.write('Original_HGVS\tHGVS_c.\tHGVS_p.\tHGVS_g._hg19\tHGVS_g._hg38\tStatus\n')
for line in converted_hgvs:
	output_file.write('\t'.join(line))
	output_file.write('\n')

output_file.close()


