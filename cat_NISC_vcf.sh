#!/bin/bash

# Intented to be run after trim_NISC_vcf.sh and run_vcf-merge.sh
# Just glues together the snv and div files, then recompresses and creates tabix index
snv=$1
div=$2

zcat $snv | head -n 1000 | grep ^#  > /tmp/$snv.header
zcat $snv | grep -v ^# > /tmp/$snv.temp
zcat $div | grep -v ^# >> /tmp/$snv.temp

sort -k1,1d -k2,2n /tmp/$snv.temp > /tmp/$snv.ready

cat /tmp/$snv.header /tmp/$snv.ready > ${snv%.snv.vcf.gz}.snv-div.vcf

bgzip ${snv%.snv.vcf.gz}.snv-div.vcf
tabix -p vcf ${snv%.snv.vcf.gz}.snv-div.vcf.gz


