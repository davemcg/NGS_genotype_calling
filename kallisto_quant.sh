#!/bin/bash

module load kallisto/0.42.4

kallisto quant -i /data/mcgaugheyd/genomes/GRCh38/0.42.4_Homo_sapiens.GRCh38.rel84.cdna.all.idx \
			   -o ${1%_*}_kallisto \
               -b 100 \
               -t $SLURM_CPUS_PER_TASK \
			   $1 $2
