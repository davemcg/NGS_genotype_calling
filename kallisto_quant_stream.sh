#!/bin/bash

module load kallisto/0.42.4
# example input: full ftp url
# will stream to kallisto, not saving fastq
# ftp.sra.ebi.ac.uk/vol1/fastq/ERR030/ERR030887/ERR030887_1.fastq.gz
output_folder=$1
fastq_url=$2
kallisto quant -i /data/mcgaugheyd/genomes/GRCh38/0.42.4_Homo_sapiens.GRCh38.rel84.cdna.all.idx \
			   -o $1 \
               -b 100 \
               -t $SLURM_CPUS_PER_TASK \
			   <(curl -s $fastq_url | zcat ) \
			   <(curl -s ${fastq_url%_1.fastq.gz}_2.fastq.gz | zcat )
