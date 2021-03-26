#!/bin/bash

set -e

module load GATK/3.8-0
module load vt/0.57721

#sbatch --cpus-per-task=8 --mem=16g --time=6:00:00 /data/OGVFB/OGL_NGS/NGS_genotype_calling/gvcf_to_vcf_070419.sh OGLv1_062619.vcf /data/OGVFB/OGL_NGS/bed/OGL731_v1.sorted.bed /data/guanb/MiSeq/190626_0061_CGWJ7/OGLv1_062619.ped
#submit at the project folder.
#this version worked well with multi-allelic locus, ped, bed needs full path

output_vcf_name=$1
bed=$2
ped=$3

mkdir -p GATKprioritization

cd gvcfs

ls *.g.vcf.gz > gvcfs.list
#this file's extension has to be .list 

# Merges all GVCFs into a VCF
GATK -m 16g GenotypeGVCFs \
	-R /data/guanb/resource/GATK_resource_bundle/b37-2.8/human_g1k_v37_decoy.fasta \
	-o $1 \
	-V gvcfs.list \
	--pedigree $ped

cat $1 | vt decompose -s - > ${1%.vcf}.decomposed.vcf
	
# Extracts all SNPs
GATK -m 16g SelectVariants \
	-R /data/guanb/resource/GATK_resource_bundle/b37-2.8/human_g1k_v37_decoy.fasta \
	-V ${1%.vcf}.decomposed.vcf \
	-L $bed \
	--interval_padding 100 \
	-selectType SNP \
	-o ${1%.vcf}.rawSNP.vcf

# Extracts all INDELS
GATK -m 16g SelectVariants \
	-R /data/guanb/resource/GATK_resource_bundle/b37-2.8/human_g1k_v37_decoy.fasta \
	-V ${1%.vcf}.decomposed.vcf \
	-L $bed \
	--interval_padding 100 \
	-selectType INDEL \
	-o ${1%.vcf}.rawINDEL.vcf

# Hard filters on GATK best practices for SNPs with MQ mod from bcbio
GATK -m 16g VariantFiltration \
	-R /data/guanb/resource/GATK_resource_bundle/b37-2.8/human_g1k_v37_decoy.fasta \
	-V ${1%.vcf}.rawSNP.vcf \
	--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
	--filterName "FAIL_McGaughey_SNP_filter_v01" \
	-o ${1%.vcf}.filterSNP.vcf

# Hard filters on GATK bests for INDELs
GATK -m 16g VariantFiltration \
	-R /data/guanb/resource/GATK_resource_bundle/b37-2.8/human_g1k_v37_decoy.fasta \
	-V ${1%.vcf}.rawINDEL.vcf \
	--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
	--filterName "FAIL_McGaughey_INDEL_filter_v01" \
	-o ${1%.vcf}.filterINDEL.vcf
	
GATK -m 16g CombineVariants \
	-R /data/guanb/resource/GATK_resource_bundle/b37-2.8/human_g1k_v37_decoy.fasta \
	--variant ${1%.vcf}.filterSNP.vcf \
	--variant ${1%.vcf}.filterINDEL.vcf \
	-o ${1%.vcf}.filterSNP-INDEL.vcf \
	--genotypemergeoption UNSORTED

rm *.decomposed.vcf*
rm *.rawSNP.vcf*
rm *.rawINDEL.vcf*
rm *.filterSNP.vcf*
rm *.filterINDEL.vcf*

# Assumes a GATK processed VCF after GATK-recommended filtering
# Hard coded against grch37

#sbatch --mem=16g --time=12:0:0 ~/git/variant_prioritization/GATK_vcf_intervar_0626.sh VCF

module load samtools/1.9

bgzip -f ${1%.vcf}.filterSNP-INDEL.vcf
tabix -p vcf ${1%.vcf}.filterSNP-INDEL.vcf.gz

mv vcf ${1%.vcf}.filterSNP-INDEL.vcf.gz* ../GATKprioritization/.
