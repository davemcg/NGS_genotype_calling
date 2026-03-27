#!/bin/bash
#SBATCH --gres=lscratch:200
#SBATCH --cpus-per-task=28
#SBATCH --mem=64g
#SBATCH --partition=norm
#SBATCH --time=48:0:0

config=$1
sample=$2

#selecting mapped reads using subset_mapped_reads.sh reduced sample processing time to 44 hours.
#selecting reads using subset_chr1-M_reads.sh: 40 hours (without P flag. Adding P flag makes the run > 48h).

set -e
mkdir -p sample_bam/subset
module load $(grep "^samtools_version:" $config | head -n 1 | cut -d"'" -f 2)
samtools view --threads $SLURM_CPUS_PER_TASK -b sample_bam/$sample.markDup.bam --output sample_bam/subset/$sample.markDup.bam \
	chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
samtools index -@ $SLURM_CPUS_PER_TASK sample_bam/subset/$sample.markDup.bam

if (( $(module list 2>&1 | grep " R/" | wc -l) == 1 )); then module unload R; fi
module load clinsv/1.1
refdata_path=/usr/local/apps/clinsv/ref_hg38/clinsv/refdata-b38
input_path=$PWD
project_folder=$PWD/clinSV/$sample
rm -rf $project_folder
mkdir -p $project_folder
singularity run --bind $refdata_path:/app/ref-data/refdata-b38 \
	--bind $project_folder:/app/project_folder \
	--bind $input_path:/app/input \
	/usr/local/apps/clinsv/1.1/libexec/1.1.sif \
	/app/clinsv/bin/clinsv \
	-r all -f \
	-p /app/project_folder/ \
	-i "/app/input/sample_bam/subset/$sample.markDup.bam" \
	-ref /app/ref-data/refdata-b38
cp clinSV/$sample/SVs/joined/SV-CNV.vcf.gz clinSV/$sample.clinsv.SV-CNV.vcf.gz
cp clinSV/$sample/SVs/joined/SV-CNV.vcf.gz.tbi clinSV/$sample.clinsv.SV-CNV.vcf.gz.tbi
cp clinSV/$sample/SVs/joined/SV-CNV.PASS.vcf clinSV/$sample.clinsv.SV-CNV.PASS.vcf
cp clinSV/$sample/SVs/joined/SV-CNV.RARE_PASS_GENE.vcf clinSV/$sample.clinsv.SV-CNV.RARE_PASS_GENE.vcf
cp clinSV/$sample/results/$sample.QC_report.pdf clinSV/$sample.QC_report.pdf

bgzip clinSV/$sample.clinsv.SV-CNV.PASS.vcf
bgzip clinSV/$sample.clinsv.SV-CNV.RARE_PASS_GENE.vcf
tabix -f -p vcf clinSV/$sample.clinsv.SV-CNV.PASS.vcf.gz
tabix -f -p vcf clinSV/$sample.clinsv.SV-CNV.RARE_PASS_GENE.vcf.gz

module load $(grep "^annotsv_version:" $config | head -n 1 | cut -d"'" -f 2)
AnnotSV -SVinputFile clinSV/$sample.clinsv.SV-CNV.PASS.vcf.gz \
	-SVinputInfo 1 -genomeBuild $(grep "^genomeBuild:" $config | head -n 1 | cut -d"'" -f 2) \
	-outputDir clinSV/$sample \
	-outputFile clinSV/$sample/$sample.clinSV.PASS.annotated.tsv
mv clinSV/$sample/$sample.clinSV.PASS.annotated.tsv clinSV/$sample.clinSV.PASS.annotated.tsv
AnnotSV -SVinputFile clinSV/$sample.clinsv.SV-CNV.RARE_PASS_GENE.vcf.gz \
	-SVinputInfo 1 -genomeBuild $(grep "^genomeBuild:" $config | head -n 1 | cut -d"'" -f 2)  \
	-outputDir clinSV/$sample \
	-outputFile clinSV/$sample/$sample.clinSV.RARE_PASS_GENE.annotated.tsv
mv clinSV/$sample/$sample.clinSV.RARE_PASS_GENE.annotated.tsv clinSV/$sample.clinSV.RARE_PASS_GENE.annotated.tsv

module load $(grep "^R_version:" $config | head -n 1 | cut -d"'" -f 2)
Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/clinsv_author_excel.R \
	$(grep "^OGL_Dx_research_genes:" $config | head -n 1 | cut -d"'" -f 2) \
	clinSV/$sample/SVs/joined/SV-CNV.RARE_PASS_GENE.xlsx \
	clinSV/$sample.clinSV.RARE_PASS_GENE.eG.tsv \
	clinSV/$sample.clinSV.RARE_PASS_GENE.eG.filtered.xlsx

if [ ! -e clinSV/result_description.docx ];	then
	cp clinSV/$sample/results/result_description.docx clinSV/result_description.docx
fi