#!/bin/sh

#GIAB files downloaded based on build_GIAB_VEP_instructions.sh in the same folder.
#The VCF files in gz format need to be processed by VT and then indexed for use with rtg eval.
#$1 is the vcf files to be tested.
#The output folder name is $2
#rtg eval is memory intensive, failed at 16g mem. Worked at 32 and 64g.
#Should be very quick to finish.
#sbatch --partition=quick --mem=64g ~/git/NGS_genotype_calling/NGS_generic_OGL/rtg_eval_vcf_GIAB.sh vcf.gz output_folder
#$3 could be /data/OGL/GIAB/NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.bed
#512g memory is sufficient for rtg of a WGS dataset. 

set -e
module load  vt/0.57721
module load samtools/1.9

vcfinput=$1
filename=$(basename $1)
output=$2

#if bed given, then use it
if [ ! -z "$3" ]; then
	bed="$3"
#otherwise use the default of bed from intercept of GIAB_hign-confidence and OGLv1.
else
	bed="/data/OGL/GIAB/OGLv1_GIAB_highconf.bed"
fi

zcat $1 | vt decompose -s - | vt normalize -n -r /data/OGVFB/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta - | bgzip -c > ${filename%.vcf.gz}.vt.vcf.gz
tabix -f -p vcf ${filename%.vcf.gz}.vt.vcf.gz

zcat ${filename%.vcf.gz}.vt.vcf.gz | head -500 | grep "#" > header

module load GATK/3.8-0

GATK -m 16g SelectVariants \
	-R /data/guanb/resource/GATK_resource_bundle/b37-2.8/human_g1k_v37_decoy.fasta \
	-V ${filename%.vcf.gz}.vt.vcf.gz \
	-L $bed \
	--interval_padding 50 \
	-selectType SNP \
	-o ${filename%.vcf.gz}.SNP1.vcf.gz

zgrep -v "#" ${filename%.vcf.gz}.SNP1.vcf.gz | cat header - | bgzip -c >  ${filename%.vcf.gz}.SNP.vcf.gz
tabix -f -p vcf ${filename%.vcf.gz}.SNP.vcf.gz

# Extracts all INDELS
GATK -m 16g SelectVariants \
	-R /data/guanb/resource/GATK_resource_bundle/b37-2.8/human_g1k_v37_decoy.fasta \
	-V ${filename%.vcf.gz}.vt.vcf.gz \
	-L $bed \
	--interval_padding 50 \
	-selectType INDEL \
	-o ${filename%.vcf.gz}.INDEL1.vcf.gz

zgrep -v "#" ${filename%.vcf.gz}.INDEL1.vcf.gz | cat header - | bgzip -c >  ${filename%.vcf.gz}.INDEL.vcf.gz

tabix -f -p vcf ${filename%.vcf.gz}.INDEL.vcf.gz

module load rtg/3.8.4

rtg vcfeval -b /data/OGL/GIAB/NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.SNP.vcf.gz -c ${filename%.vcf.gz}.SNP.vcf.gz -t /data/OGL/resources/genomes/hg19/human_g1k_v37decoy_sdf --ref-overlap --evaluation-regions $bed -o $2-SNP

rtg vcfeval -b /data/OGL/GIAB/NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.INDEL.vcf.gz -c ${filename%.vcf.gz}.INDEL.vcf.gz -t /data/OGL/resources/genomes/hg19/human_g1k_v37decoy_sdf --ref-overlap --evaluation-regions $bed -o $2-INDEL

#
#--vcf-score-field QUAL 
#Default vcf-score-field is GQ.

# rm ${filename%.vcf.gz}.vt.vcf.gz*
# rm ${filename%.vcf.gz}.SNP.vcf.gz*
# rm ${filename%.vcf.gz}.INDEL.vcf.gz*
# rm header


rtg rocplot $2-SNP/snp_roc.tsv.gz --png=$2-SNP/$2.snp.png
rtg rocplot $2-SNP/non_snp_roc.tsv.gz --png=$2-SNP/$2.indel.png
rtg rocplot $2-SNP/weighted_roc.tsv.gz --png=$2-SNP/$2.roc.png
rtg rocplot -P $2-SNP/weighted_roc.tsv.gz --png=$2-SNP/$2.precission_recall.png

rtg rocplot $2-INDEL/snp_roc.tsv.gz --png=$2-INDEL/$2.snp.png
rtg rocplot $2-INDEL/non_snp_roc.tsv.gz --png=$2-INDEL/$2.indel.png
rtg rocplot $2-INDEL/weighted_roc.tsv.gz --png=$2-INDEL/$2.roc.png
rtg rocplot -P $2-INDEL/weighted_roc.tsv.gz --png=$2-INDEL/$2.precission_recall.png

#If multisample vcf files, --sample daughter; Use <baseline_sample>,<calls_sample> to select different sample names for baseline and calls.
#Interpretation of the results: Note that vcfeval reports true positives both counted using the baseline variant representation as well as counted using the call variant representation. When these numbers differ greatly, it indicates a general difference in representational conventions used between the two call sets. Since false negatives can only be measured in terms of the baseline representation, sensitivity is defined as: Sensitivity = TPbaseline/(TPbaseline + FN). Conversely since false positives can only be measured in terms of the call representation, precision is defined as: Precision = TPcall/(TPcall + FP).
