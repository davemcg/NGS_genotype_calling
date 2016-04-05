#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import argparse
import subprocess
import collections
from collections import defaultdict
import sys

parser = argparse.ArgumentParser()
parser.add_argument('filelist', help= 
	"Takes in list of NISC lane bams, transfer the bams over \
	from Trek, extracts fastq, realigns fastq with bwa-mem \
	to 1000g phaseII hs37d5 with read group info extracted \
	from the lane bam \
	\
	Invoke script with sbatch --mem=20G --cpus-per-task 10 \
	\
	Example (of the input list): \
		CCGO_800014	150223_OPTIMUS_C6HGHANXX.8.9477663 \
		CCGO_800014	150306_YOSHI_C6HDMANXX.5.9477663 \
		CCGO_800014	150223_OPTIMUS_C67T9ANXX.7.9477663 \
		CCGO_800014	150223_OPTIMUS_C67T9ANXX.8.9477663 \
		CCGO_800014	150223_OPTIMUS_C67T9ANXX.6.9477663 \
		CCGO_800015	150223_OPTIMUS_C67T9ANXX.7.9477645 \
		CCGO_800015	150223_OPTIMUS_C67T9ANXX.8.9477645 \
		CCGO_800015	150223_OPTIMUS_C67T9ANXX.6.9477645 \
		CCGO_800015	150306_YOSHI_C6HDMANXX.5.9477645 \
		CCGO_800015	150223_OPTIMUS_C6HGHANXX.8.9477645 \
	\
	Example (exome): \
   		sbatch --mem=20G --cpus-per-task 10 this_script.py list_of_laneBams.txt")
args = parser.parse_args()
bamlist = args.filelist

# loop through the bamlist, making a dict
sample_laneBam = collections.defaultdict(list)
for line in open(bamlist):
	line = line.split()
	if len(line) != 2:
		print("Input list formatting issue\n\n", line)
		sys.exit(0)
	sample_laneBam[line[0]].append(line[1])

# create subfolders for each sample
# Use scp to move lane Bams from Trek to biowulf2
samples = []
[samples.append(k) for k,v in sample_laneBam.items()]
samples.sort()
for one_sample in samples:
	mkdir_call = "mkdir " + one_sample
	subprocess.check_call(mkdir_call, shell=True)
	# create directory structure
	base_dir = "/cluster/ifs/projects/solexa/reads/"
	# loop through each laneBam and do the scp
	for laneBam in sample_laneBam[one_sample]:
		folder = laneBam.split('.')[0]
		full_dir = base_dir + folder + '/' + laneBam + '*'
		scp_call = 'scp trek.nhgri.nih.gov:' + full_dir + ' ' + one_sample + '/'
		subprocess.check_call(scp_call, shell=True)
		
		# extract read group info from bam
		samtools_header_call = 'samtools view -H ' + laneBam + '.bam ' +  '| grep ^@RG'
		read_group = subprocess.check_call(samtools_header_call, shell=True)
		# example read_group:
		# '@RG\tID:3\tSM:AMS0013\tLB:AMS0013_1\tPU:160129_OPTIMUS_C7NP2ANXX.3.11457313\tCN:NISC\tDT:2016-02-04T20:49:40-05:00\tPL:HiSeq2000\n'
		# replace SM with sample name
		read_group = read_group.split('\t')
		read_group[2] = 'SM:' + one_sample
		# https://github.com/IARCbioinfo/BAM-tricks
		# https://github.com/samtools/samtools/issues/532#issuecomment-205877064
		# now shuffle, extract interleaved fastq, and run bwa-mem
		big_run = 'samtools collate -uOn 128 ' + laneBam + '.bam' + 'TMP' + laneBam + \
				  '| samtools fastq - | bwa mem -M -t 10 -B 4 -O 6 -E 1 -M -R ' + \
				  "\"" + '\t'.join(read_group) + "\"" + \
			 	  ' /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta \
				  - | samtools view -1 - > ' + laneBam + '.bwa-mem.b37.bam' 

# merge laneBams into one bam for GATK workflow
