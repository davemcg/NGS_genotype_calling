#!/bin/bash

# Assumes a GATK processed VCF after GATK-recommended filtering
# Hard coded against hg19/grch37
# 
module load vcftools

VCF=$1
cores=$2
mkdir tmp

#only keep AF > 0.25
vcffilter -f "AF > 0.25" $VCF > tmp/$VCF.AFfiltered

# vt to "left align and trim alternative variants"
cat tmp/$VCF.AFfiltered \
	| sed 's/ID=AD,Number=./ID=AD,Number=R/' \
	| ~/git/vt/./vt decompose -s - \
	| ~/git/vt/./vt normalize -r /fdb/GATK_resource_bundle/hg19-2.8/ucsc.hg19.fasta - \
	> tmp/$VCF 

# annotate with VEP82
/home/mcgaugheyd/bin/run_VEP.sh tmp/$VCF GRCh37 $cores

# move out of tmp folder
mv tmp/${VCF%.vcf}.VEPRCh37.vcf* . 

# compress and index
bgzip ${VCF%.vcf}.VEP.GRCh37.vcf
tabix -p vcf ${VCF%.vcf}.VEP.GRCh37.vcf.gz

rm -rf tmp
