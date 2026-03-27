#!/bin/bash
#SBATCH --gres=lscratch:200
#SBATCH --cpus-per-task=2
#SBATCH --mem=72g
#SBATCH --partition=norm
#SBATCH --time=10:0:0

set -e
config=$1
sample=$2

module load JAX-CNV/20240208
RefGenome=$(grep "^ref_genome:" $config | head -n 1 | cut -d"'" -f 2)
JAX-CNV GetCnvSignal -f $RefGenome \
-k ${RefGenome%/*}/genome.kmer \
-b sample_bam/$sample.markDup.bam -o jax-cnv/$sample
JaxCNVMerge.R -i jax-cnv/$sample

module load $(grep "^annotsv_version:" $config | head -n 1 | cut -d"'" -f 2)
AnnotSV -SVinputFile jax-cnv/$sample.merge.bed \
	-SVinputInfo 1 -genomeBuild $(grep "^genomeBuild:" $config | head -n 1 | cut -d"'" -f 2) \
	-outputDir jax-cnv \
	-outputFile jax-cnv/$sample.jaxcnv.annotated.0.tsv
awk -F "\t" 'BEGIN{OFS="\t"} {print $0, "$sample"}' jax-cnv/$sample.merge.bed \
	> jax-cnv/$sample.merge.bed.temp && mv jax-cnv/$sample.merge.bed.temp jax-cnv/$sample.merge.bed

module load $(grep "^R_version:" $config | head -n 1 | cut -d"'" -f 2)
Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/jaxCNV.R \
	jax-cnv/$sample.jaxcnv.annotated.0.tsv \
	$sample jax-cnv/$sample.jaxcnv.annotated.tsv
rm jax-cnv/$sample jax-cnv/$sample.jaxcnv.annotated.0.tsv

