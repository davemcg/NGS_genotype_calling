#!/bin/bash

module load GATK/3.4-46

input_bam=$1

# Takes ~ 90 minutes
GATK -m 8g RealignerTargetCreator \
	-R /fdb/GATK_resource_bundle/hg19-2.8/ucsc.hg19.fasta \
	-I $input_bam \
	--known /fdb/GATK_resource_bundle/hg19-2.8/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz \
	-o ${input_bam%.bam}.forIndexRealigner.intervals \
	-L /data/mcgaugheyd/genomes/hg19/SeqCap_EZ_Exome_v3_primary.bed \
	--interval_padding 100

# Takes ~ 100 minutes
GATK -m 8g IndelRealigner \
	-R /fdb/GATK_resource_bundle/hg19-2.8/ucsc.hg19.fasta \
	-I $input_bam \
	-known /fdb/GATK_resource_bundle/hg19-2.8/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz \
	-targetIntervals ${input_bam%.bam}.forIndexRealigner.intervals \
	-o ${input_bam%.bam}.realigned.bam

 
GATK -m 8g BaseRecalibrator \
	-R /fdb/GATK_resource_bundle/hg19-2.8/ucsc.hg19.fasta \
	-I ${input_bam%.bam}.realigned.bam \
	-knownSites /fdb/GATK_resource_bundle/hg19-2.8/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz \
	-o ${input_bam%.bam}.recal_data.table

GATK -m 8g AnalyzeCovariates \
	-R /fdb/GATK_resource_bundle/hg19-2.8/ucsc.hg19.fasta \
	-BQSR ${input_bam%.bam}.recal_data.table \
	-plots ${input_bam%.bam}.BQSRplots.pdf

# Takes ~ 180 minutes
GATK -m 8g HaplotypeCaller \
	-R /fdb/GATK_resource_bundle/hg19-2.8/ucsc.hg19.fasta \
	-I ${input_bam%.bam}.realigned.bam \
	-L /data/mcgaugheyd/genomes/hg19/SeqCap_EZ_Exome_v3_primary.bed \
	--interval_padding 100 \
	--emitRefConfidence GVCF \
	-BQSR ${input_bam%.bam}.recal_data.table \
	-o ${input_bam%.bam}.realigned.raw.g.vcf
