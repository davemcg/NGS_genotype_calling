#!/bin/bash

# Merges vcf files together with ls and xargs.
# Need to be able to specify exact list with ls and wildcards
# Expect tabix and gzipped files
# Usage:
# sbatch run_vcf-merge.sh \\*mpg.snv.trimmed.vcf.gz CCGO_800067-68-69-70-71.mpg.snv.vcf
# Need to double escape "*"

module load vcftools


search_term=$1
output=$2

ls $1 | xargs vcf-merge > $output
bgzip $output
tabix -p vcf $output.gz
