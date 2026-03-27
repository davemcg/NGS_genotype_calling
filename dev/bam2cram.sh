#!/bin/bash
#SBATCH --gres=lscratch:300
#SBATCH --cpus-per-task=16
#SBATCH --mem=16g
#SBATCH --partition=quick
#SBATCH --time=2:0:0

#took < 15 min for a genome

module load samtools/1.21
ref_genome=$1
bam_file=$2
sample=$(basename $2)
mkdir -p cram
samtools view -T $ref_genome --threads $SLURM_CPUS_PER_TASK --output-fmt cram,store_md=1,store_nm=1 -o cram/${sample%.bam}.cram $bam_file
samtools index -@ $SLURM_CPUS_PER_TASK cram/${sample%.bam}.cram cram/${sample%.bam}.crai