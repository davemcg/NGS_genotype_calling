#!/bin/bash

# outside (i.e. not NISC) trio of exomes

# first create fastq from the bam
for i in *bam; do sbatch ~/bin/run_bam2fastq $i; done

# align
sbatch --cpus-per-task 10 ~/bin/exome_workflow_v02/run_bwa-mem_hg37d5.sh 16_114290_1.fastq 16_114290_2.fastq \\@RG\\\\tID:ambry_16_114290\\\\tSM:16_114290\\\\tPL:ILLUMINA 16_114290.realigned.g1k_decoy.bwa.bam
sbatch --cpus-per-task 10 ~/bin/exome_workflow_v02/run_bwa-mem_hg37d5.sh 16_114291_1.fastq 16_114291_2.fastq \\@RG\\\\tID:ambry_16_114291\\\\tSM:16_114291\\\\tPL:ILLUMINA 16_114291.realigned.g1k_decoy.bwa.bam
sbatch --cpus-per-task 10 ~/bin/exome_workflow_v02/run_bwa-mem_hg37d5.sh 16_391899_1.fastq 16_391899_2.fastq \\@RG\\\\tID:ambry_16_391899\\\\tSM:16_391899\\\\tPL:ILLUMINA 16_391899.realigned.g1k_decoy.bwa.bam

# create GVCF
sbatch --time=36:00:00 --mem=20G ~/bin/exome_workflow_v02/process_and_callGVCF_noB.sh 16_114290.realigned.g1k_decoy.bwa.bam 
sbatch --time=36:00:00 --mem=20G ~/bin/exome_workflow_v02/process_and_callGVCF_noB.sh 16_114291.realigned.g1k_decoy.bwa.bam
sbatch --time=36:00:00 --mem=20G ~/bin/exome_workflow_v02/process_and_callGVCF_noB.sh 16_391899.realigned.g1k_decoy.bwa.bam

# make ped and gvcfs.list file

# create merged and filterd VCF
sbatch --time=8:00:00 GVCF_to_hardFilteredVCF_noB.sh gvcfs.list out.vcf the.ped

