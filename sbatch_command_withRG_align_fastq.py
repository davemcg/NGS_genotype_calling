#!/usr/local/Anaconda/envs/py3.4.3/bin/python

"""
Creates bwa-mem hg19 command to align fastq files
with proper RG info
"""
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('file', help='Input bam file to generate read group information')
args = parser.parse_args()
bamfile = args.file

# Runs samtools view -h. Needs samtools(1.2) already loaded in bash environment)
samtools_input = 'samtools view -h ' + bamfile + '| head -n 100 | grep ^@RG'
samtools_view = (subprocess.check_output(samtools_input, shell=True)).decode('utf-8')

info = samtools_view.split('\t')

# Builds the new RG from file name and NISC provided info from their bam
ID = 'ID:' + info[4].split(':')[1]
SM = 'SM:' + bamfile.split('.')[0]
LB = 'LB:' + info[4].split(':')[1].split('.')[2]
PL = 'PL:Illumina\\" \\'
Output = SM + '.bwa-mem.hg19.bam'
# Joins all together
RG_core = '\\\\t'.join(['\\"\@RG',ID, SM, LB, PL])

core = bamfile.split('.')[0]

# runs alignment!
run_bwa = 	('sbatch --mem=50G --cpus-per-task=10 run_bwa-mem_hg19.sh ' +
			' ' + core + '_1.fastq' + ' ' +  core + '_2.fastq ' +
			'"@RG\\t' + ID + '\\t' + SM + '\\t' + LB + '\\t' + 'PL:Illumina\" ' +
			core + '.bwa-mem.hg19.bam')
print(run_bwa, sep='')
subprocess.check_output(run_bwa, shell=True).decode('utf-8')

