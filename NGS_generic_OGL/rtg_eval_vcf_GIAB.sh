#!/bin/sh

#GIAB files downloaded based on build_GIAB_VEP_instructions.sh in the same folder.
#The VCF files in gz format need to be processed by VT and then indexed for use with rtg eval.
#$1 is the vcf files to be tested.
#The output is in the folder named as the basename of the vcf file
#rtg eval is memory intensive, failed at 16g mem. Tried at 64g and worked. Other memory size was not tested. 
#Should be very quick to finish.
#sbatch --partition=quick --mem=64g ~/git/NGS_genotype_calling/NGS_generic_OGL/rtg_eval_vcf_GIAB.sh vcf.gz output_folder

module load  vt/0.577
module load samtools/1.9

vcfinput=$1
filename=$(basename $1)
output=$2

zcat $1 | vt decompose -s - | vt normalize -r /data/OGVFB/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta - | bgzip -c > ${filename%.vcf.gz}.vt.vcf.gz

tabix -f -p vcf ${filename%.vcf.gz}.vt.vcf.gz

module load rtg/3.8.4

rtg vcfeval -b /data/OGL/GIAB/NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.vt.vcf.gz -c ${filename%.vcf.gz}.vt.vcf.gz -t /data/OGL/resources/genomes/hg19/human_g1k_v37decoy_sdf -e /data/OGL/GIAB/OGLv1_GIAB_highconf.bed --vcf-score-field QUAL -o $2

rm ${filename%.vcf.gz}.vt.vcf.gz*

#If multisample vcf files, --sample daughter; Use <baseline_sample>,<calls_sample> to select different sample names for baseline and calls.
#Interpretation of the results: Note that vcfeval reports true positives both counted using the baseline variant representation as well as counted using the call variant representation. When these numbers differ greatly, it indicates a general difference in representational conventions used between the two call sets. Since false negatives can only be measured in terms of the baseline representation, sensitivity is defined as: Sensitivity = TPbaseline/(TPbaseline + FN). Conversely since false positives can only be measured in terms of the call representation, precision is defined as: Precision = TPcall/(TPcall + FP).
