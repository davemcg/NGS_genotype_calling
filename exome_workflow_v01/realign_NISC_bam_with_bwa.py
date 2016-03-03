#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('file', help= 
	"Takes in NISC bam file, extracts fastsq, \
	realigns fastq with bwa-mem to hg19, with \
	read group info extracted from NISC bam \
	Requires bam2fastq/1.1.0 and samtools/1.2 to be loaded \
	\
	Invoke script with sbatch --mem=50G --cpus-per-task 10 \
	\
	Example: \
   		sbatch --mem=50G --cpus-per-task this_script.py A_BAM_FILE_001.bam")
args = parser.parse_args()
bamfile = args.file
bam_name = bamfile.split('.')[0]

# Extract fastq
print("Beginning fastq extraction")
bam2fastq_call = "bam2fastq -o " + bam_name + "#.fastq " + bamfile
subprocess.check_call(bam2fastq_call, shell=True)
print("Done")

# Runs samtools view -h and captures output
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


# bwa alignment
print("BWA run beginning")
run_bwa =   ('run_bwa-mem_hg19.sh ' +
            bam_name + '_1.fastq ' + bam_name + '_2.fastq ' +
            '\\@RG\\\\t' + ID + '\\\\t' + SM + '\\\\t' + LB + '\\\\t' + 'PL:Illumina ' +
            bam_name + '.bwa-mem.hg19.bam')
print(run_bwa)
subprocess.check_call(run_bwa, shell=True)
print("All done!")
