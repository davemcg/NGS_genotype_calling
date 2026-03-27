#!/bin/bash
#SBATCH --gres=lscratch:50
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --partition=quick
#SBATCH --time=00:30:00

config=$1
sample=$2
echo $sample
set -e

module load $(grep "^samtools_version:" $config | head -n 1 | cut -d"'" -f 2)
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
