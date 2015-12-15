#!/bin/bash

# Takes the huge NISC snv and div vcfs and only keeps non-reference variants
# retains header and recompresses with bgzip and creates tabix index

input=$1

zcat $input | head -n 1000 | grep ^# > ${input%.vcf.gz}.trimmed.vcf
zcat $input | awk '$5!="." {print $0}' >> ${input%.vcf.gz}.trimmed.vcf
bgzip ${input%.vcf.gz}.trimmed.vcf 
tabix -p vcf ${input%.vcf.gz}.trimmed.vcf.gz
