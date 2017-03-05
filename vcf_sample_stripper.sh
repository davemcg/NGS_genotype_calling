#!/bin/bash

# removes INFO field from vcf and outputs a user given sample
# comma separate mlutiple samples

module load samtools/1.3.1

input_vcf=$1
samples=$2
output_vcf=$3

bcftools annotate -x INFO $input_vcf | bcftools view - -s $samples | bgzip > $output_vcf
tabix -p vcf $output_vcf
