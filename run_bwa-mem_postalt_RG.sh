#!/bin/bash

module load samtools

fastq1=$1
fastq2=$2
rg=$3
output_bam_name=$4

# parse out Read Group 
RG=`samtools view -h $ | head -n 1000 | grep ^@RG | cut -f5 | cut -f2 -d":"`

~/bin/bwa.kit/bwa mem -t 10 -B 4 -O 6 -E 1 -M -R $rg /data/mcgaugheyd/genomes/GRCh38/hs38DH.fa $1 $2 | \
	~/bin/bwa.kit/k8 ~/bin/bwa.kit/bwa-postalt.js /data/mcgaugheyd/genomes/GRCh38/hs38DH.fa.alt | \
	samtools view -1 - > $output_bam_name

