#!/bin/bash

##################
# Takes NISC bam file, extracts fastq, aligns with bwa
# Processes bam with various Picard tools to sort, mark dups, and index
# Calls GVCF for the bam

# The resulting GVCF can be combined with other GVCFs with run_GATK_GVCFtoFilteredVCF.sh 
##################

# modules needed for python j1 job
module load samtools/1.3
module load bwa/0.7.12
module load picard/2.1.1

NISC_laneBam_matrix=$1 
swarm_job_name=$2
exome_bait_bed=$3 #Give full path

############
#RE-ALIGNMENT
# pulls bwa-formatted info (read group info, bam file, etc) then hands over to sbatch 
j1=$(sbatch --job-name bwa.$1 --mem=30g --cpus-per-task 10 ~/bin/exome_workflow_v02/realign_NISC_laneBams_with_bwa.py $1 -S $2 -B $3)
############

############
# BAM processing and GVCF calling
# sort, mark duplicates, and index bam
j2=$(swarm -f $2 --dependency=afterok:$j1 --time 24:00:00 -g 20 --module picard/2.1.1,GATK/3.5-0)
############

