#!/bin/bash
#SBATCH --gres=lscratch:200
#SBATCH --cpus-per-task=56
#SBATCH --mem=128g
#SBATCH --partition=norm
#SBATCH --time=16:0:0

config=$1
sample=$2

set -e

module load $(grep "^clair3_version:" $config | head -n 1 | cut -d"'" -f 2)
module load $(grep "^whatshap_version:" $config | head -n 1 | cut -d"'" -f 2)
module load $(grep "^samtools_version:" $config | head -n 1 | cut -d"'" -f 2)

LWORK_DIR=/lscratch/$SLURM_JOB_ID

clair3 --bam_fn sample_bam/$sample.markDup.bam \
 --ref_fn=$(grep "^ref_genome:" $config | head -n 1 | cut -d"'" -f 2) \
 --threads=$SLURM_CPUS_PER_TASK --platform=ilmn --gvcf \
 --model_path=/data/OGL/resources/clair3/ilmn \
 --use_whatshap_for_final_output_phasing \
 --sample_name=$sample \
 --output=/lscratch/$SLURM_JOB_ID/
 
cp $LWORK_DIR/merge_output.gvcf.gz clair3/gvcf/$sample.clr3.gvcf.gz
cp $LWORK_DIR/merge_output.gvcf.gz.tbi clair3/gvcf/$sample.clr3.gvcf.gz.tbi
cp $LWORK_DIR/phased_merge_output.vcf.gz clair3/vcf/$sample.clr3.phased.vcf.gz
cp $LWORK_DIR/phased_merge_output.vcf.gz.tbi clair3/vcf/$sample.clr3.phased.vcf.gz.tbi

bcftools norm --multiallelics -any --output-type u clair3/vcf/$sample.clr3.phased.vcf.gz \
	| bcftools norm -d exact --output-type u - \
	| bcftools filter --threads $(($SLURM_CPUS_PER_TASK-2)) --include 'filter=="PASS" & FORMAT/AD[0:1]>2' --output-type z --output clair3/vcf/$sample.clr3.filtered.vcf.gz
sleep 2
tabix -f -p vcf clair3/vcf/$sample.clr3.filtered.vcf.gz

module load $(grep "^bedtools_version:" $config | head -n 1 | cut -d"'" -f 2)
module load $(grep "^R_version:" $config | head -n 1 | cut -d"'" -f 2)

mkdir -p /lscratch/$SLURM_JOB_ID/AutoMap
zcat clair3/vcf/$sample.clr3.filtered.vcf.gz > /lscratch/$SLURM_JOB_ID/AutoMap/$sample.vcf
bash /data/OGL/resources/git/AutoMap/AutoMap_v1.2.sh \
	--vcf /lscratch/$SLURM_JOB_ID/AutoMap/$sample.vcf \
	--out AutoMap --genome hg38 --chrX
echo "AutoMap1.2 done"

if [[ $(grep -v ^# AutoMap/$sample/$sample.HomRegions.tsv | wc -l) -eq 0 ]]; then
	touch AutoMap/$sample/$sample.HomRegions.bed
	touch AutoMap/$sample/$sample.HomRegions.annot.tsv
	echo "no ROH region detected."
else
	grep -v ^# AutoMap/$sample/$sample.HomRegions.tsv | cut -f 1-3 > AutoMap/$sample/$sample.HomRegions.bed
	module load $(grep "^annotsv_version:" $config | head -n 1 | cut -d"'" -f 2)
	AnnotSV -SVinputFile AutoMap/$sample/$sample.HomRegions.bed \
		-SVinputInfo 1 -genomeBuild GRCh38 \
		-outputDir AutoMap/$sample \
		-outputFile AutoMap/$sample/$sample.annotated.tsv
	Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/automap.R AutoMap/$sample/$sample.HomRegions.tsv AutoMap/$sample/$sample.annotated.tsv $(grep "^OGL_Dx_research_genes:" $config | head -n 1 | cut -d"'" -f 2) AutoMap/$sample/$sample.HomRegions.annot.tsv
	rm AutoMap/$sample/$sample.annotated.tsv
fi

rm AutoMap/$sample/$sample.HomRegions.tsv AutoMap/$sample/$sample.HomRegions.bed
#merge_output.gvcf.gz phase?
#merge_output.vcf.gz
#phased_merge_output.vcf.gz

# Used 10G lscrach, 40G mem, 63 max CPU, 40 min (including phasing and haplotagged bam)
# module load clair3/20250303
# AddPairEndAlleleDepth.py --bam_fn sample_bam/D2097_01_BP524859.markDup.bam \
	# --clair3_vcf_input /lscratch/$SLURM_JOB_ID/D2097_01_BP524859/phased_merge_output.vcf.gz \
	# --vcf_output /lscratch/$SLURM_JOB_ID/D2097_01_BP524859/phased_merge_output.PEAD.vcf.gz \
	# --threads $SLURM_CPUS_PER_TASK

#
#PEAD is time consuming. Started at 10:28am;

