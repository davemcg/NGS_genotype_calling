#!/bin/bash

# to run snakemake as batch job
# run in the data folder for this project
# on biowulf2:
# /data/mcgaugheyd/projects/nei/mcgaughey/EGA_EGAD00001002656
# uses config.yaml in git directory

mkdir -p 00log
module load snakemake || exit 1

sbcmd="sbatch --cpus-per-task={threads} \
--mem={cluster.mem} \
--time={cluster.time} \
--partition={cluster.partition} \
--output={cluster.output} \
--error={cluster.error} \
{cluster.extra}"


snakemake -s /home/mcgaugheyd/git/NGS_genotype_calling/NGS_generic/Snakefile \
-pr --local-cores 2 --jobs 1999 \
--cluster-config /home/mcgaugheyd/git/NGS_genotype_calling/NGS_generic/cluster.json \
--cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
-k --restart-times 4 \
--resources parallel=4 \
--configfile $1
