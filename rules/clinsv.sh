#!/bin/bash
#SBATCH --gres=lscratch:200
#SBATCH --cpus-per-task=28
#SBATCH --mem=64g
#SBATCH --partition=norm
#SBATCH --time=32:0:0

config=$1
sample=$2

#selecting mapped reads using subset_mapped_reads.sh reduced sample processing time to 44 hours.
#selecting reads using subset_chr1-M_reads.sh ???
#The bam file has to be saved in the directory sample_bam and named as *.markDup.bam, unless changing the code below.
set -e
if (( $(module list 2>&1 | grep " R/" | wc -l) == 1 )); then module unload R; fi
module load clinsv/1.1
refdata_path=/usr/local/apps/clinsv/ref_hg38/clinsv/refdata-b38
input_path=$PWD
project_folder=/lscratch/$SLURM_JOB_ID/clinSV/$sample
mkdir -p $project_folder
singularity run --bind $refdata_path:/app/ref-data/refdata-b38 \
	--bind $project_folder:/app/project_folder \
	--bind $input_path:/app/input \
	/usr/local/apps/clinsv/1.1/libexec/1.1.sif \
	/app/clinsv/bin/clinsv \
	-r all -f \
	-p /app/project_folder/ \
	-i "/app/input/sample_bam/$sample.markDup.bam" \
	-ref /app/ref-data/refdata-b38
mv $project_folder/SVs/joined/SV-CNV.vcf.gz clinSV/$sample/$sample.clinsv.SV-CNV.vcf.gz
mv $project_folder/SVs/joined/SV-CNV.vcf.gz.tbi clinSV/$sample/$sample.clinsv.SV-CNV.vcf.gz.tbi
mv $project_folder/SVs/joined/SV-CNV.PASS.vcf clinSV/$sample/$sample.clinsv.SV-CNV.PASS.vcf
mv $project_folder/SVs/joined/SV-CNV.RARE_PASS_GENE.vcf clinSV/$sample/$sample.clinsv.SV-CNV.RARE_PASS_GENE.vcf
mv $project_folder/results/$sample.QC_report.pdf clinSV/$sample/$sample.QC_report.pdf

module load $(grep "^samtools_version:" $config | head -n 1 | cut -d"'" -f 2)
bgzip clinSV/$sample/$sample.clinsv.SV-CNV.PASS.vcf
bgzip clinSV/$sample/$sample.clinsv.SV-CNV.RARE_PASS_GENE.vcf
tabix -f -p vcf clinSV/$sample/$sample.clinsv.SV-CNV.PASS.vcf.gz
tabix -f -p vcf clinSV/$sample/$sample.clinsv.SV-CNV.RARE_PASS_GENE.vcf.gz

module load $(grep "^annotsv_version:" $config | head -n 1 | cut -d"'" -f 2)
AnnotSV -SVinputFile clinSV/$sample/$sample.clinsv.SV-CNV.PASS.vcf.gz \
	-SVinputInfo 1 -genomeBuild $(grep "^genomeBuild:" $config | head -n 1 | cut -d"'" -f 2) \
	-outputDir $project_folder \
	-outputFile clinSV/$sample/$sample.clinSV.PASS.annotated.tsv
mv clinSV/$sample/$sample.clinSV.PASS.annotated.tsv clinSV/$sample.clinSV.PASS.annotated.tsv
AnnotSV -SVinputFile clinSV/$sample.clinsv.SV-CNV.RARE_PASS_GENE.vcf.gz \
	-SVinputInfo 1 -genomeBuild $(grep "^genomeBuild:" $config | head -n 1 | cut -d"'" -f 2)  \
	-outputDir $project_folder \
	-outputFile clinSV/$sample/$sample.clinSV.RARE_PASS_GENE.annotated.tsv

module load $(grep "^R_version:" $config | head -n 1 | cut -d"'" -f 2)
Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/clinsv_author_excel.R \
	$(grep "^OGL_Dx_research_genes:" $config | head -n 1 | cut -d"'" -f 2) \
	$project_folder/SVs/joined/SV-CNV.RARE_PASS_GENE.xlsx \
	clinSV/$sample/$sample.clinSV.RARE_PASS_GENE.eG.tsv \
	clinSV/$sample/$sample.clinSV.RARE_PASS_GENE.eG.filtered.xlsx

if [ ! -e clinSV/result_description.docx ];	then
	cp $project_folder/results/result_description.docx clinSV/result_description.docx
fi