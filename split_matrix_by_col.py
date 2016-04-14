#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import argparse
import os.path
import sys

parser = argparse.ArgumentParser()
parser.add_argument('file', nargs='?',\
	help = 
	"Takes in a file (intended for use for the NISC lane bam \
	exome to GVCF pipeline matrix) and splits it into separate \
	files by the first column. \
	\
	Input File contents: \
	\
	Sample1 131129_YOSHI_C2P4RACXX.3.5489791 \
	Sample1 131129_YOSHI_C2P4RACXX.3.5489802 \
	Sample2 131129_YOSHI_C2P4RACXX.1.548979 \
	\
	This script will output two files: \
	Sample1.laneBam.matrix \
		Sample1 131129_YOSHI_C2P4RACXX.3.5489791 \
		Sample1 131129_YOSHI_C2P4RACXX.3.5489802 \
	Sample2.laneBam.matrix \
		Sample2 131129_YOSHI_C2P4RACXX.1.548979")

args = parser.parse_args()
# if file given
if args.file:
	file = open(args.file)
#if no file given, then use stdin 
else:
	file = []
	file.append(sys.stdin.read())
	
for line in file:
	s_line = line.split()
	sample = s_line[0]
	sample_file_name = sample + ".laneBam.matrix"
	# if file doesn't exist, then create it and write the line
	if not os.path.isfile(sample_file_name):
		sample_file = open(sample_file_name, 'w')
		sample_file.write(line)
	# if file already exists, then just write line
	else:
		sample_file.write(line)
