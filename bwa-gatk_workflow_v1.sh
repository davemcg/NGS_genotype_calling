#!/bin/bash

##################
# Takes NISC bam file, extracts fastq, aligns with bwa
# Processes bam with various Picard tools to sort, mark dups, and index
# Calls GVCF for the bam

# The resulting GVCF can be combined with other GVCFs with run_GATK_GVCFtoFilteredVCF.sh 
##################

module load samtools/1.2
module load bam2fastq/1.1.0

nisc_input_bam=$1 

############
#RE-ALIGNMENT
# pulls bwa-formatted info (read group info, bam file, etc) then hands over to sbatch 
bwa_command=$(sbatch --dependency=afterok:$j1 ~/bin/sbatch_command_withRG_align_fastq_bare.py  $1 -o $1.temp)
j2=$(sbatch --dependency=afterok:$bwa_command --mem=50g --cpus-per-task 10 $bwa_command)
############

############
# BAM processing
# sort, mark duplicates, and index bam
j3=$(sbatch --dependency=afterok:$j2 --mem=20G ~/bin/run_Picard_SortMarkIndexBam.sh ${1%.bam}.bwa-mem.hg19.bam)
############

############
# Call GVCF
# Individual GVCF file, which in conjuction with other GVCFs and filtering (hard or VQSR) are used
# to make a VCF file with called genotypes
j4=$(sbatch --dependency=afterok:$j3 --mem=16g --time=6:00:00 ~/bin/run_GATK_callGVCFs.sh ${1%.bam}.bwa-mem.hg19.sorted.markDup.bam)
############
