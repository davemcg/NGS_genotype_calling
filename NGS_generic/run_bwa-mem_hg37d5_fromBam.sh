#!/bin/bash

module load bwa/0.7.12
module load samtools/1.4

input_bam=$1
rg=$2
output_name=$3

samtools collate -uOn 128 $input_bam /tmp/TMP-$input_bam | \
	samtools fastq - | \
	bwa mem -M -t $SLURM_CPUS_PER_TASK -B 4 -O 6 -E 1 -M -p -R $rg /data/mcgaugheyd/genomes/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta - | \
	samtools view -1 - > $output_name

