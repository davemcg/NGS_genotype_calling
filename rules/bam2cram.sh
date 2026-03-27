#!/bin/bash
#SBATCH -c16
#SBATCH --mem=16g
#SBATCH --gres=lscratch:100
#SBATCH --time=1:0:0
#SBATCH --partition=quick

#for bamfile in bam/*.bam; do sbatch ~/git/NGS_genotype_calling/rules/bam2cram.sh config_generic.yaml $bamfile; done
#Took 7-20 min/2GB mem for GS

set -e

config_file=$1 # exome/es/wes or genome/gs/wgs
bam_input=$2

ref_genome=$(grep "^ref_genome:" $config_file | head -n 1 | cut -d"'" -f 2)
module load $(grep "^samtools_version:" $config_file | head -n 1 | cut -d"'" -f 2) 

mkdir -p cram
sample=$(basename $bam_input | cut -d. -f 1)
samtools view -T $ref_genome --threads 16 --output-fmt cram,store_md=1,store_nm=1 -o cram/$sample.cram $bam_input
samtools index -@ 16 cram/$sample.cram cram/$sample.cram.crai
