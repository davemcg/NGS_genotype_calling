#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import argparse
import os.path

parser = argparse.ArgumentParser()
parser.add_argument('file', help = 
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
file = args.file

for line in open(file):
	s_line = line.split()
	sample = s_line[0]
	sample_file_name = sample + ".laneBam.matrix"
	if not os.path.isfile(sample_file_name):
		sample_file = open(sample_file_name, 'w')
		sample_file.write(line)
	else:
		sample_file.write(line)
