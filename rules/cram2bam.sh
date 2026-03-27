#!/bin/bash
#SBATCH -c56
#SBATCH --mem=16g
#SBATCH --gres=lscratch:100
#SBATCH --time=1:0:0
#SBATCH --partition=quick

#for bamfile in bam/*.bam; do sbatch ~/git/NGS_genotype_calling/rules/bam2cram.sh config_generic.yaml $bamfile; done
#Took 7-20 min/2GB mem for GS

set -e

config_file=$1 # exome/es/wes or genome/gs/wgs
sample=$2

ref_genome=$(grep "^cram_ref:" $config_file | head -n 1 | cut -d"'" -f 2)
module load $(grep "^samtools_version:" $config_file | head -n 1 | cut -d"'" -f 2) 

mkdir -p sample_bam
if [ -e cram/$sample.crai ] || [ -e cram/$sample.cram.crai ] ; then
 echo "index present"
 else samtools index -@ $SLURM_CPUS_PER_TASK cram/$sample.cram
fi
samtools view -T $ref_genome -@ $SLURM_CPUS_PER_TASK \
 -b -o sample_bam/$sample.bam cram/$sample.cram
samtools index -@ $SLURM_CPUS_PER_TASK sample_bam/$sample.bam