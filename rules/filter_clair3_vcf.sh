#!/bin/bash
#SBATCH -c8
#SBATCH --mem=8g
#SBATCH --gres=lscratch:100
#SBATCH --time=4:0:0
set -e
config_file=$1
sample=$2
ref_genome=$(grep "^ref_genome:" $config_file | head -n 1 | cut -d"'" -f 2)

module load samtools

bcftools norm --threads 6 --multiallelics -any --output-type u clair3/vcf/$sample.clr3.phased.vcf.gz \
 | bcftools norm --threads 6 -d exact --check-ref s --fasta-ref $ref_genome --output-type u - \
 | bcftools filter --threads 6 --write-index=tbi --include 'QUAL>10 & FORMAT/AD[0:1]>3' --output-type z --output clair3/vcf/$sample.clr3.filtered.vcf.gz
