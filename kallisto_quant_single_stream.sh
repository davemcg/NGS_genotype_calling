#!/bin/bash

module load kallisto/0.42.4

folder=$1
fastq_url=$2
frag_length=$3
frag_sd=$4

kallisto quant -i /data/mcgaugheyd/genomes/GRCh38/0.42.4_Homo_sapiens.GRCh38.rel84.cdna.all.idx \
			   -o $folder \
               -b 100 \
               -t $SLURM_CPUS_PER_TASK \
			   --single -l $frag_length -s $frag_sd  \
			   <(curl -s $fastq_url | zcat )
