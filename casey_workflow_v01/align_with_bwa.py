#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('file', help= 
	"Takes in fastq.gz files from Casey panels,
	aligns with bwa-mem to 1000g phaseII hs37d5 with \
	read group info  \
	Requires samtools/1.3 and bwa to be loaded \
	\
	Invoke script with sbatch --mem=30G --cpus-per-task 10 \
	\
	Example (exome): \
   		sbatch --mem=30G --cpus-per-task 10 this_script.py A_BAM_FILE_001.bam")
args = parser.parse_args()
bamfile = args.file
bam_name = bamfile.split('.bam')[0]

# ID is basically the NISC file name for the fastq (minus the .fq at the end)
ID = 'ID:' + fastq_file[0:-3]
# SM is the bam file name or sample name. Needs to the same for each sample!
SM = 'SM:' + bam_name
# LB is the library
LB = 'LB:' + fastq_file.split('.')[2]
PL = 'PL:Illumina\\" \\'

Output = SM + '.bwa-mem.b37.bam'
# Joins all together
RG_core = '\\\\t'.join(['\\"\@RG',ID, SM, LB, PL])


# bwa alignment
print("BWA run beginning")
run_bwa =   ('/home/mcgaugheyd/bin/exome_workflow_v02/run_bwa-mem_hg37d5.sh ' +
            bam_name + '_1.fastq.gz ' + bam_name + '_2.fastq.gz ' +
            '\\@RG\\\\t' + ID + '\\\\t' + SM + '\\\\t' + LB + '\\\\t' + 'PL:Illumina ' +
            bam_name + '.bwa-mem.b37.bam')
print(run_bwa)
subprocess.check_call(run_bwa, shell=True)
print("All done!")
