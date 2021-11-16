#!/bin/bash
#SBATCH --gres=lscratch:800
#SBATCH --cpus-per-task=56
#SBATCH --mem=1000g
#SBATCH --time=12:0:0
#SBATCH --partition=largemem

#lscratch 200g failed - used /scratch below instead.
outputVcf=$1 # i.e. novogen133.glnexus.vcf.gz
module load glnexus/1.2.7 samtools/1.13
mkdir -p glnexus
time glnexus --dir /lscratch/$SLURM_JOB_ID/glnexus --config DeepVariantWGS --threads 56 --mem-gbytes 800 deepvariant/gvcf/*.g.vcf.gz | bcftools view - | bgzip -c > glnexus/$outputVcf
#--threads 56 - -Oz -o glnexus1/$outputVcf
tabix -f -p vcf glnexus/$outputVcf

#Took 	904 GB	8:47:12; 
#Exceed 1T memory if piping to bcftools view --threads 56 - -Oz -o glnexus/$outputVcf
#After switching to hg38, make bcf files instead. Then split to ~100 files in the variant_prioritization step.