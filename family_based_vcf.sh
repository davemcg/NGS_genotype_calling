#!/bin/bash

# takes a large, multi-sample and creates vcf files for each family
# requires a ped file for the family relationships and the multi-sample
# vcf


module load samtools

vcf_file=$1
ped_file=$2

samples=$(bcftools query -l $vcf_file)

for i in $samples
	do 	echo $i
		fam=`grep $i $ped_file | cut -f1 -d' '`
		echo $fam | grep -f - $ped_file | awk '{print $2}' > ${fam}.samples
	done

for i in *.samples
	do	echo $i
		bcftools view --force-samples -S $i $vcf_file > ${i%.samples}.b37.bwa-mem.hardFilterSNP-INDEL.vcf
		bgzip ${i%.samples}.b37.bwa-mem.hardFilterSNP-INDEL.vcf
		tabix -p vcf ${i%.samples}.b37.bwa-mem.hardFilterSNP-INDEL.vcf.gz
	done


