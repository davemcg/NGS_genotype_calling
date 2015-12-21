

"""
Takes two files, creates a key based on multiple columns (0-based, separated by commas),
builds a dictionary, merges the two files on the key, and prints output
"""

import sys

file1 = open(sys.argv[1])
file2 = open(sys.argv[2])
columns = str(sys.argv[3])

def key_columns(columns_input):
	cols = columns_input.split(',')
	col_list = [int(i) for i in cols]
	return(col_list)

def build_key(col_list, sline):
	key_list = [sline[i] for i in col_list]
	key = '_'.join(key_list)
	return(key)

def build_dict(file, col_list):
	file1_dict = {}
	for line in file.readline():
		sline = line.split('\t')
		key = build_key(col_list, sline)
		if key in file1_dict:
			print("Non-unique key in first file")
			break
		else:
			file1_dict[key] = line
	return(file1_dict)

	
def main():
	col_list = key_columns(columns)
	file1_dict = build_dict(file1, col_list)
	for line in file2.readline():
		sline = line.split('\t')
		key = build_key(col_list, sline)	
		if key in file1_dict:
			print(file1_dict[key], line, sep="\t")
		

main()
