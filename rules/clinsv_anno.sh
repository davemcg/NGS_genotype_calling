#!/bin/bash
#SBATCH --gres=lscratch:50
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --partition=quick
#SBATCH --time=1:0:0

config=$1
sample=$2
#module load parallel
#cut -d',' -f1 CLIA25-10_metadata_file.csv | sort -u | parallel --no-run-if-empty -j 8 --bar 'bash ~/git/NGS_genotype_calling/rules/clinsv_anno.sh config_generic.yaml {}'
#OR using "while" command
#cut -d, -f1 CLIA25-10_metadata_file.csv | sort -u > sample.list && while read -r sample; do sbatch ~/git/NGS_genotype_calling/rules/clinsv_anno.sh config_generic.yaml $sample; done < sample.list

set -e
if [[ $(module list 2>&1 | grep "annotsv" | wc -l) -lt 1 ]]; then module load $(grep "^annotsv_version:" $config | head -n 1 | cut -d"'" -f 2); fi
project_folder=/lscratch/$SLURM_JOB_ID/clinSV/$sample
mkdir -p $project_folder
AnnotSV -SVinputFile clinSV/$sample/$sample.clinsv.SV-CNV.PASS.vcf.gz \
	-SVinputInfo 1 -genomeBuild GRCh38 \
	-outputDir $project_folder \
	-outputFile $sample.clinSV.PASS.annotated.tsv
gzip -c $project_folder/$sample.clinSV.PASS.annotated.tsv > clinSV/$sample/$sample.clinSV.PASS.annotated.tsv.gz
AnnotSV -SVinputFile clinSV/$sample/$sample.clinsv.SV-CNV.RARE_PASS_GENE.vcf.gz \
	-SVinputInfo 1 -genomeBuild GRCh38 \
	-outputDir $project_folder \
	-outputFile $sample.clinSV.RARE_PASS_GENE.annotated.tsv
gzip -c $project_folder/$sample.clinSV.RARE_PASS_GENE.annotated.tsv > clinSV/$sample/$sample.clinSV.RARE_PASS_GENE.annotated.tsv.gz

if [[ $(module list 2>&1 | grep " R/" | wc -l) -lt 1 ]]; then module load $(grep "^R_version:" $config | head -n 1 | cut -d"'" -f 2); fi
Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/clinsv_author_excel.R \
	/data/OGL/resources/OGLpanelGeneDxORcandidate.xlsx \
	clinSV/$sample/$sample.clinsv.SV-CNV.RARE_PASS_GENE.xlsx \
	clinSV/$sample/$sample.clinSV.RARE_PASS_GENE.eG.tsv.gz \
	clinSV/$sample/$sample.clinSV.RARE_PASS_GENE.eG.filtered.xlsx
