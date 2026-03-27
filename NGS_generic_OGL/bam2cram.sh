#!/bin/bash
#SBATCH -c8
#SBATCH --mem=4g
#SBATCH --gres=lscratch:100
#SBATCH --time=1:0:0
set -e

module load samtools/1.13
BAMFILE=$1
filename=$(basename $1)
samtools view -T /data/OGL/resources/genomes/GRCh38/GRCh38Decoy2.fa --threads 8 --output-fmt cram,store_md=1,store_nm=1 -o /lscratch/$SLURM_JOB_ID/${filename%.markDup.bam}.cram $BAMFILE

mkdir -p cram

rm bam/${filename%.markDup.bam}*
rm sample_bam/${filename%.markDup.bam}*

mv /lscratch/$SLURM_JOB_ID/${filename%.markDup.bam}.cram cram/${filename%.markDup.bam}.cram 

samtools index -@ 8 cram/${filename%.markDup.bam}.cram cram/${filename%.markDup.bam}.crai

