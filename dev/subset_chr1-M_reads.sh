#!/bin/bash
#SBATCH -c18
#SBATCH --mem=18g
#SBATCH --gres=lscratch:200
#SBATCH --time=4:0:0

input=$1
output=$2
module load samtools/1.21
samtools view --threads $SLURM_CPUS_PER_TASK -b $input --output $output chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
samtools index -@ $SLURM_CPUS_PER_TASK $output