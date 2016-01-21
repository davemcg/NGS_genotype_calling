#!/bin/bash

module load GATK/3.5-0

gvcfs_list=$1
output_vcf_name=$2

# Merges all GVCFs into a VCF
GATK -m 8g GenotypeGVCFs \
	-R /fdb/GATK_resource_bundle/hg19-2.8/ucsc.hg19.fasta \
	-o $2 \
	-V $gvcfs_list

# Extracts all SNPs
GATK -m 8g SelectVariants \
	-R /fdb/GATK_resource_bundle/hg19-2.8/ucsc.hg19.fasta \
	-V $2 \
	-L /data/mcgaugheyd/genomes/hg19/SeqCap_EZ_Exome_v3_primary.bed \
    --interval_padding 100 \
	-selectType SNP \
	-o ${2%.vcf}.rawSNP.vcf

# Extracts all INDELS
GATK -m 8g SelectVariants \
	-R /fdb/GATK_resource_bundle/hg19-2.8/ucsc.hg19.fasta \
	-V $2 \
	-L /data/mcgaugheyd/genomes/hg19/SeqCap_EZ_Exome_v3_primary.bed \
	--interval_padding 100 \
	-selectType INDEL \
	-o ${2%.vcf}.rawINDEL.vcf

# Hard filters on GATK best practices for SNPs with MQ mod from bcbio
GATK -m 8g VariantFiltration \
	-R /fdb/GATK_resource_bundle/hg19-2.8/ucsc.hg19.fasta \
	-V ${2%.vcf}.rawSNP.vcf \
	--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
	--filterName "FAIL_McGaughey_SNP_filter_v01" \
	-o ${2%.vcf}.filterSNP.vcf

# Hard filters on GATK bests for INDELs
GATK -m 8g VariantFiltration \
	-R /fdb/GATK_resource_bundle/hg19-2.8/ucsc.hg19.fasta \
	-V ${2%.vcf}.rawINDEL.vcf \
	--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
	--filterName "FAIL_McGaughey_INDEL_filter_v01" \
	-o ${2%.vcf}.filterINDEL.vcf
	
GATK -m 8g CombineVariants \
	-R /fdb/GATK_resource_bundle/hg19-2.8/ucsc.hg19.fasta \
	--variant ${2%.vcf}.filterSNP.vcf \
	--variant ${2%.vcf}.filterINDEL.vcf \
	-o ${2%.vcf}.filterSNP-INDEL.vcf \
	--genotypemergeoption UNSORTED
