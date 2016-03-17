#!/bin/bash

# Move files to keep to folder
# Then remove the rest (except other folders)
# Want: 
#	original (NISC) bam
#	final bam
#	g.vcf
#	vcf


mkdir core_files
mkdir scripts

# will find the original bam files (which don't have "bwa" in their name) and move to core_files
find -mindepth 1 -maxdepth 1 | grep -v bwa.*bam | grep bam$ | xargs  -I{} mv {} core_files/
# moves last bam
mv *hg19.sorted.markDup.realigned.bam* core_files/
# moves g.vcf file. May need for re-doing genotype calls
mv *.g.vcf* core_files/
# moves vcf file ready for use (clinician, VEP, gemini)
mv *GATK.filterSNP-INDEL.vcf* core_files/

# other useful files
mv *.db core_files/	#gemini database
mv *.ped core_files/	#pedigree
mv *.list core_files/	#gvcf list files. Hand-made
mv *.sh scripts/ 	#wouldn't want to toss any scripts automatically
