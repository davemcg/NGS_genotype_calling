#!/bin/bash

# to run snakemake as batch job 
# to go from GVCFs to vcf
# run in the data folder for this project
# 
# Example config.yaml file:
# ~/git/NGS_genotype_calling/src/Snakefile_gvcf_to_vcf_example_config.yaml

module load snakemake || exit 1

mkdir -p 00log 

sbcmd="sbatch --cpus-per-task={threads} \
--mem={cluster.mem} \
--time={cluster.time} \
--partition={cluster.partition} \
--output={cluster.output} \
--error={cluster.error} \
{cluster.extra}"

config_yaml=$1

if [ -z "$1" ]; then
	echo "config.yaml file not given"
	exit 1
fi


snakemake -s ~/git/NGS_genotype_calling/src/Snakefile_gvcf_to_vcf \
-pr --local-cores 2 --jobs 1999 \
--cluster-config ~/git/NGS_genotype_calling/src/cluster_gvcf_to_vcf.json \
--configfile $config_yaml \
--cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
-k \
--resources parallel=2

# --restart-times 4 \
