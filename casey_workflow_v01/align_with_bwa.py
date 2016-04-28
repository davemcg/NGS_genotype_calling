#!/usr/local/Anaconda/envs/py3.4.3/bin/python

import argparse
import subprocess

parser = argparse.ArgumentParser(description= \
	"Takes in fastq.gz files from Casey panels, \
	aligns with bwa-mem to 1000g phaseII hs37d5 with \
	read group info  \
	Requires samtools/1.3 and bwa to be loaded \
	\
	Invoke script with sbatch --mem=30G --cpus-per-task 10 \
	\
	Example (exome): \
   		sbatch --mem=30G --cpus-per-task 10 this_script.py A_BAM_FILE_001.bam")
parser.add_argument('forward', help = 'Forward fastq reads')
parser.add_argument('reverse', help = 'Reverse fastq reads')
parser.add_argument('CaseyID', help = 'Casey sample ID, e.g. 15-00414')
parser.add_argument('panel', help = 'Casey Panel')
args = parser.parse_args()
forward = args.forward
reverse = args.reverse
caseyID = args.CaseyID
panel = args.panel

# ID is basically the NISC file name for the fastq (minus the .fq at the end)
ID = 'ID:' + forward.split('.fastq')[0]
# SM is the bam file name or sample name. Needs to the same for each sample!
SM = 'SM:' + caseyID + '_' + panel
# LB is the library
LB = 'LB:' + forward.split('.fastq')[0]
PL = 'PL:Illumina\\" \\'

Output = SM + '.bwa-mem.b37.bam'
# Joins all together
RG_core = '\\\\t'.join(['\\"\@RG',ID, SM, LB, PL])


# bwa alignment
print("BWA run beginning")
run_bwa =   ('/home/mcgaugheyd/bin/exome_workflow_v02/run_bwa-mem_hg37d5.sh ' +
            forward + ' ' + reverse + ' ' + 
            '\\@RG\\\\t' + ID + '\\\\t' + SM + '\\\\t' + LB + '\\\\t' + 'PL:Illumina ' +
            caseyID + '/' + caseyID + '_' + panel  + '.bwa-mem.b37.bam')
print(run_bwa)
subprocess.check_call(run_bwa, shell=True)
print("All done!")
