#!/bin/sh

#GIAB files downloaded based on build_GIAB_VEP_instructions.sh in the same folder.
#The VCF files in gz format need to be processed by VT and then indexed for use with rtg eval.
#$1 is the vcf files to be tested.
#The output folder name is $2
#rtg eval is memory intensive, failed at 16g mem. Worked at 32 and 64g.
#Should be very quick to finish.
#sbatch --partition=quick --mem=64g ~/git/NGS_genotype_calling/NGS_generic_OGL/rtg_eval_vcf_GIAB.sh vcf.gz output_folder
#$3 could be /data/OGL/GIAB/NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.bed
#GIAB true set did not contain any block_substitution. 

module load  vt/0.57721
module load samtools/1.9

vcfinput=$1
filename=$(basename $1)
output=$2

# if bed given, then use it
if [ ! -z "$3" ]; then
	bed="$3"
# otherwise use the default of bed from intercept of GIAB_hign-confidence and OGLv1.
else
	bed="/data/OGL/GIAB/OGLv1_GIAB_highconf.bed"
fi


zcat $1 | vt decompose -s - | vt decompose_blocksub -a - | vt normalize -n -r /data/OGL/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta - | bgzip -c > ${filename%.vcf.gz}.vt.vcf.gz

tabix -f -p vcf ${filename%.vcf.gz}.vt.vcf.gz

module load rtg/3.8.4

rtg vcfeval -b /data/OGL/GIAB/NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.vt.vcf.gz -c ${filename%.vcf.gz}.vt.vcf.gz -t /data/OGL/resources/genomes/hg19/human_g1k_v37decoy_sdf --ref-overlap --evaluation-regions $bed -o $2

#
#--vcf-score-field QUAL 
#Default vcf-score-field is GQ.

rm ${filename%.vcf.gz}.vt.vcf.gz*

rtg rocplot $2/snp_roc.tsv.gz --png=$2/$2.snp.png
rtg rocplot $2/non_snp_roc.tsv.gz --png=$2/$2.indel.png
rtg rocplot $2/weighted_roc.tsv.gz --png=$2/$2.roc.png
rtg rocplot -P $2/weighted_roc.tsv.gz --png=$2/$2.precission_recall.png


#If multisample vcf files, --sample daughter; Use <baseline_sample>,<calls_sample> to select different sample names for baseline and calls.
#Interpretation of the results: Note that vcfeval reports true positives both counted using the baseline variant representation as well as counted using the call variant representation. When these numbers differ greatly, it indicates a general difference in representational conventions used between the two call sets. Since false negatives can only be measured in terms of the baseline representation, sensitivity is defined as: Sensitivity = TPbaseline/(TPbaseline + FN). Conversely since false positives can only be measured in terms of the call representation, precision is defined as: Precision = TPcall/(TPcall + FP).
