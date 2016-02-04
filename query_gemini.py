#!/usr/local/Anaconda/envs/py3.4.3/bin/python


"""
Runs autosomal_recessive, de_novo, mendel_errors, comp_hets, autosomal_dominant, and roh
tests on select families (trios) in a gemini database

Returns results in a xlsx file for the clinicians/genetic counselors
"""


import argparse
import subprocess
import xlsxwriter

#########PARSER##############
parser = argparse.ArgumentParser(description=\
    "Queries gemini (v0.18) database to identify variants matching \
	models of automosomal recessive, de novo, mendelian error, \
	compound hets and autosomal dominant. \
	\
	Returns an xlsx file of the results.\
	\
	Input: \
		database to query \
		family to analyze (use a '-' if skipping, max one family) \
		name for output xlsx file.\
	\
	Examples (no need for sbatch): \
		query_gemini.py CCG0.gemini.db CCG0_800042 CCGO_800042.variants.xlsx \
		query_gemini.py CCGO.gemini.db - all.variants.xlsx")

parser.add_argument('--database')
parser.add_argument('--family')

#########CODE#############
def autosomal_recessive(db, family):
	filter = " --filter \"max_aaf_all < 0.01 AND (is_coding=1 OR is_splicing=1) \
				AND filter IS NULL\" --gt-pl-max 10 -d 5 --min-gq 20 "
	if family=='-':
		ar_query = "gemini autosomal_recessive" + columns + db + " " + filter
	else:
		ar_query = "gemini autosomal_recessive" + columns + \
					"--families " + family + " " + db + " " + filter
	ar = subprocess.check_output(ar_query,shell=True).decode('utf-8')
	ar = ar.split('\n')
	return(ar)

def de_novo(db, family):
	filter = " --filter \"max_aaf_all < 0.005 AND (is_coding=1 OR is_splicing=1) \
				AND filter IS NULL\" --gt-pl-max 10 -d 5 --min-gq 20 "
	if family=="-":
		dn_query = "gemini de_novo" + columns + db + " " + filter
	else:
		dn_query = "gemini de_novo" + columns + \
					"--families " + family + " " + db + " " + filter
	dn = subprocess.check_output(dn_query,shell=True).decode('utf-8')	
	dn = dn.split('\n')
	return(dn)

def mendel_errors(db, family):
	# gemini v0.18 has two bugs with this call:
		# Can't parse by family
		# Can't call exac numbers
		# I can work around the first, but not the second (easily)
	filter = " --filter \"max_aaf_all < 0.005 AND (is_coding=1 OR is_splicing=1) \
				AND filter IS NULL\" --gt-pl-max 1 -d 5 --min-gq 20"
	columns = 	" --columns \"chrom, start, end, codon_change, aa_change, type, impact, \
				impact_severity, gene, vep_hgvsp, aaf_1kg_all, aaf_exac_all, \
				gerp_bp_score, polyphen_score, \
				cadd_raw, sift_pred, sift_score, vep_grantham, (gt_depths).(*) \" "
	
	me_query = "gemini mendel_errors" + columns + db + " " + filter 	
	me = subprocess.check_output(me_query,shell=True).decode('utf-8')
	me = me.split('\n')
	# get the header in
	me_out = []
	me_out.append(me[0]) 
	# filter for only the family we want
	for line in me:
		s_line = line.split('\t')
		if line and s_line[46]==family:
			me_out.append(line)
	return(me_out)
	
def output_to_xlsx(data,sheet_name):
	worksheet = workbook.add_worksheet(sheet_name)
	row = 0
	col = 0
	for line in data:
		line = line.split('\t')
		for unit in line: 
			worksheet.write(row, col, unit)
			col += 1
		col = 0
		row += 1

def main():
	args = parser.parse_args()
	db = args.database
	family = args.family
	
	# output time
	ar = autosomal_recessive(db, family)
	output_to_xlsx(ar, "Autosomal Recessive")	

	dn = de_novo(db, family)
	output_to_xlsx(dn, "De Novo")	
	
	me = mendel_errors(db, family)
	output_to_xlsx(me, "Mendelian Errors")
	workbook.close()
	


# global stuff
workbook = xlsxwriter.Workbook('test.xlsx')
columns = 	" --columns \"chrom, start, end, codon_change, aa_change, type, impact, \
			impact_severity, gene, vep_hgvsp, aaf_1kg_all, aaf_exac_all, \
			exac_num_hom_alt, exac_num_het, gerp_bp_score, polyphen_score, \
			cadd_raw, sift_pred, sift_score, vep_grantham, (gt_depths).(*) \" "

# run it!
main()
