#!/bin/bash

# Assumes a GATK processed VCF after GATK-recommended filtering
# Hard coded against hg19/grch37

VCF=$1

cat $VCF \
	| sed 's/ID=AD,Number=./ID=AD,Number=R/' \
	| ~/git/vt/./vt decompose -s - \
	| ~/git/vt/./vt normalize -r /fdb/GATK_resource_bundle/hg19-2.8/ucsc.hg19.fasta - \
	> test.vcf

/home/mcgaugheyd/bin/run_VEP82.sh test.vcf GRCh37

bgzip test.VEP82everything.GRCh37.vcf

tabix -p vcf test.VEP82everything.GRCh37.vcf.gz
