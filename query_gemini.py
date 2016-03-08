#!/usr/local/Anaconda/envs/py3.4.3/bin/python

# Yes, this uses python3. 3.4.3, to be specific. Should work on other
# versions of python3, but I haven't tested 



import argparse
from argparse import RawTextHelpFormatter
import subprocess
import xlsxwriter

#########PARSER##############
parser = argparse.ArgumentParser(description=\
    """
	Queries gemini (v0.18) database to identify variants matching 
	models of automosomal recessive, de novo, mendelian error, 
	compound hets and autosomal dominant. 
	
	Returns an xlsx file of the results.
	
	Input: 
		database to query 
		family to analyze (optional) 
		name for output xlsx file.
	
	Examples (no need for sbatch): 
		query_gemini.py --database CCG0.gemini.db --family CCG0_800042 CCGO_800042.variants.xlsx 
		query_gemini.py --database CCGO.gemini.db --family all.variants.xlsx""", 

	formatter_class=RawTextHelpFormatter)

parser.add_argument('-d','--database',required=True)
parser.add_argument('-f','--family')
parser.add_argument('-o','--output_name', required=True)
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
	return(ar,ar_query)

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
	return(dn, dn_query)

def mendel_errors(db, family):
	# gemini v0.18 has a bug with this call:
		# Can't parse by family
		# Hence my workaround
	filter = " --filter \"max_aaf_all < 0.005 AND (is_coding=1 OR is_splicing=1) \
				AND filter IS NULL\" --gt-pl-max 1 -d 5 --min-gq 20 "
	me_query = "gemini mendel_errors" + columns + db + " " + filter 	
	me = subprocess.check_output(me_query,shell=True).decode('utf-8')
	me = me.split('\n')
	if family == '-':
		me_out = me
	else:
		# Get header in, unformatted (will happen later)
		me_out = []
		me_out.append(me[0])
		# find family_id index 
		header = me[0].split('\t')	
		family_id_index = header.index('family_id') 
		# filter for only the family we want
		for line in me:
			s_line = line.split('\t')
			if line and s_line[family_id_index]==family:
				me_out.append(line)
	return(me_out, me_query)

def comp_hets(db, family):
	# can't call exac numbers in v0.18 (reported bug, fixed in next release)
	columns = 	" --columns \"chrom, start, end, codon_change, aa_change, type, impact, \
				impact_severity, gene, vep_hgvsp, aaf_1kg_all, aaf_exac_all, \
				gerp_bp_score, polyphen_score, \
				cadd_raw, sift_pred, sift_score, vep_grantham, (gt_depths).(*) \" "
	
	filter = " --filter \"max_aaf_all < 0.01 AND (is_coding=1 OR is_splicing=1) \
				AND filter IS NULL\" --gt-pl-max 10 -d 5 --min-gq 20 --max-priority 2 "
	if family == "-":
		ch_query = "gemini comp_hets" + columns + db + " " + filter
	else:
		ch_query = "gemini comp_hets" + columns + \
					"--families " + family + " " + db + " " + filter
	ch = subprocess.check_output(ch_query,shell=True).decode('utf-8')
	ch = ch.split('\n')
	return(ch, ch_query)

def autosomal_dominant(db, family):
	filter = " --filter \"max_aaf_all < 0.0001 AND (is_coding=1 OR is_splicing=1) \
				AND filter IS NULL\" --gt-pl-max 10 -d 5 --min-gq 20 "
	if family == "-":
		ad_query = "gemini autosomal_dominant" + columns + db + " " + filter
	else:
		ad_query = "gemini autosomal_dominant" + columns + \
					"--families " + family + " " + db + " " + filter
	ad = subprocess.check_output(ad_query,shell=True).decode('utf-8')
	ad = ad.split('\n')
	return(ad, ad_query)

def overview(db, queries):
	# summary stats, queries used, ped file
	stats_query = "gemini stats --gts-by-sample " + db
	stats = subprocess.check_output(stats_query,shell=True).decode('utf-8')
	stats = stats.split('\n')
	ped_query = "gemini query --header -q \"SELECT * FROM samples\" " + db
	ped = subprocess.check_output(ped_query,shell=True).decode('utf-8')
	ped = ped.split('\n')
	output = []
	output.append('Overall genotypes by sample')
	output.extend(stats)
	output.append('PED information used for calls')
	output.append('Gender: 1=male, 2=female, 0=unknown')
	output.append('Phenotype: 1=unaffected, 2=affected, 0=unknown')
	output.append("Any column after 'phenotype' is not used in these queries")
	output.extend(ped)
	output.append('Gemini Queries')
	output.extend(queries)
	return(output)

def output_to_xlsx(data,sheet_name):
	worksheet = workbook.add_worksheet(sheet_name)
	row = 0
	col = 0
	# Handling for nothing found. Don't want anyone thinking a messup happened
	# Will print first bit of info if this logic screws up
	if len(data) < 2:
		worksheet.write(0,0, "No variants found")	
		worksheet.write(1,0, data[0])
	else:		
		for line in data:
			line = line.split('\t')
			for unit in line: 
				worksheet.write(row, col, unit)
				col += 1
			col = 0
			row += 1

def main():
	db = args.database
	if args.family:
		family = args.family
	else:
		family = '-'
	# output time
	ar, ar_query = autosomal_recessive(db, family)
	output_to_xlsx(ar, "Autosomal Recessive")	

	dn, dn_query = de_novo(db, family)
	output_to_xlsx(dn, "De Novo")	
	
	me, me_query = mendel_errors(db, family)
	output_to_xlsx(me, "Mendelian Errors")

	ch, ch_query = comp_hets(db, family)
	output_to_xlsx(ch, "Compound Hets")

	ad, ad_query = autosomal_dominant(db, family)
	output_to_xlsx(ad, "Autosomal Dominant")

	# get all queries in one list
	queries = []
	queries.append(ar_query.replace('\t',' ')), queries.append(dn_query.replace('\t',' '))
	queries.append(me_query.replace('\t',' ')), queries.append(ch_query.replace('\t',' '))
	queries.append(ad_query.replace('\t',' '))

	# Create the info worksheet
	overview_info = overview(db, queries)
	output_to_xlsx(overview_info, "Info")
	workbook.close()
	


# global stuff
args = parser.parse_args()
workbook = xlsxwriter.Workbook(args.output_name)
columns = 	" --columns \"chrom, start, end, codon_change, aa_change, type, impact, \
			impact_severity, gene, vep_hgvsp, aaf_1kg_all, aaf_exac_all, \
			exac_num_hom_alt, exac_num_het, gerp_bp_score, polyphen_score, \
			cadd_scaled, sift_pred, sift_score, vep_grantham, (gt_depths).(*) \" "

# run it!
main()
