#!/bin/bash
#SBATCH --gres=lscratch:50
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --partition=quick
#SBATCH --time=2:0:0

ml samtools
ls *.vcf.gz | grep -v "phased" > temp_vcf.file.list
while read -r vcf; do bcftools norm --multiallelics -any --output-type u $vcf | bcftools norm -d exact --output-type u - | bcftools filter --threads 8 --include 'FILTER="PASS" & FORMAT/AD[0:1]>2' --output-type z --output ${vcf%.vcf.gz}.dv.hardfiltered.vcf.gz; tabix -p vcf ${vcf%.vcf.gz}.dv.hardfiltered.vcf.gz; done < temp_vcf.file.list


#normalize (left-align) for dv done at the merge_dv_hf step.