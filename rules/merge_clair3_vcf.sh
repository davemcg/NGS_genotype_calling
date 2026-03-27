#!/bin/bash
#SBATCH -c8
#SBATCH --mem=16g
#SBATCH --gres=lscratch:100
#SBATCH --time=4:0:0

config_file=$1
hf_vcf=$2
glnexus_vcf=$3
ref_genome=$(grep "^ref_genome:" $config_file | head -n 1 | cut -d"'" -f 2)
analysis_batch_name=$(grep "^analysis_batch_name:" $config_file | head -n 1 | cut -d"'" -f 2)
gt_call_version=$(grep "^gt_call_version:" $config_file | head -n 1 | cut -d"'" -f 2)

set -e
module load samtools
LWORK_DIR=/lscratch/$SLURM_JOB_ID
mkdir -p prioritization-clair3g
bcftools isec --threads 8 -p $LWORK_DIR -n=2 -w1 --collapse none --output-type z \
  $hf_vcf $glnexus_vcf

bcftools annotate --threads 8 --write-index=tbi --set-id 'clr3g_%CHROM\:%POS%REF\>%ALT' \
  -a $glnexus_vcf -c INFO \
  $LWORK_DIR/0000.vcf.gz -Oz -o prioritization-clair3g/$analysis_batch_name.$gt_call_version.clr3.vcf.gz
			