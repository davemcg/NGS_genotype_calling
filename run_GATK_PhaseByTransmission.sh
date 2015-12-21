#!/bin/bash

module load GATK/3.5-0

input_vcf=$1
ped=$2
output_vcf_name=$3
mendelian_violation_file=$4

# No need for sbatch. Will finish in a few minutes.
GATK -m 2g PhaseByTransmission \
	-R /fdb/GATK_resource_bundle/hg19-2.8/ucsc.hg19.fasta \
	-V $1 \
	-ped $2 \
	-o $3 \
	--MendelianViolationsFile $4 \
	--DeNovoPrior 1e-7 \
	--pedigreeValidationType SILENT 
