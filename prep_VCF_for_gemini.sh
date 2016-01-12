#!/bin/bash

# Assumes a GATK processed VCF after GATK-recommended filtering
# Hard coded against hg19/grch37
# 

VCF=$1
cores=$2
mkdir tmp
cat $VCF \
	| sed 's/ID=AD,Number=./ID=AD,Number=R/' \
	| ~/git/vt/./vt decompose -s - \
	| ~/git/vt/./vt normalize -r /fdb/GATK_resource_bundle/hg19-2.8/ucsc.hg19.fasta - \
	> tmp/$VCF 

/home/mcgaugheyd/bin/run_VEP82.sh tmp/$VCF GRCh37 $cores

bgzip ${VCF%.vcf}.VEP82everything.GRCh37.vcf

tabix -p vcf ${VCF%.vcf}.VEP82everything.GRCh37.vcf

rm -rf tmp/*
