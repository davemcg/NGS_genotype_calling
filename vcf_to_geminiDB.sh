#!/bin/bash

# Assumes a GATK processed VCF after GATK-recommended filtering
# Hard coded against grch37
module load vcftools
module load gemini/0.19.0

VCF=$1
PED=$2
DBNAME=$3

mkdir tmp

#only keep AF > 0.25
#vcffilter -f "AF > 0.25" $VCF > tmp/$VCF.AFfiltered

#vt to "left align and trim alternative variants"
cat $VCF \
	| sed 's/ID=AD,Number=./ID=AD,Number=R/' \
	| ~/git/vt/./vt decompose -s - \
	| ~/git/vt/./vt normalize -r /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta - \
	> tmp/${VCF%.gz}

# annotate with VEP
/home/mcgaugheyd/bin/run_VEP.sh tmp/${VCF%.gz} GRCh37 $SLURM_CPUS_PER_TASK

# move out of tmp folder
mv tmp/${VCF%.vcf.gz}.VEP.GRCh37.vcf* . 

# compress and index
bgzip ${VCF%.vcf.gz}.VEP.GRCh37.vcf
tabix -p vcf ${VCF%.vcf.gz}.VEP.GRCh37.vcf.gz

rm -rf tmp


######## Now create gemini DB #############
gemini load --cores $SLURM_CPUS_PER_TASK -t VEP -v ${VCF%.vcf.gz}.VEP.GRCh37.vcf.gz -p $PED $DBNAME
