#!/bin/bash
#SBATCH -c18
#SBATCH --mem=18g
#SBATCH --gres=lscratch:200
#SBATCH --time=4:0:0

input=$1
output=$2
module load samtools/1.21
samtools view --threads $SLURM_CPUS_PER_TASK -b -F 0x4 $input --output $output
samtools index -@ $SLURM_CPUS_PER_TASK $output