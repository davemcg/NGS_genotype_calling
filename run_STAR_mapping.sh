#!/bin/bash

module load STAR/2.4.2a

GENOME=$1
forward=$2
reverse=$3

mkdir STAR_bams

STAR \
	--runThreadN $SLURM_CPUS_PER_TASK \
	--genomeDir $GENOME \
	--readFilesIn $2 $3 \
	--outFilterType BySJout \
	--outFilterMultimapNmax 20 \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax 999 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix STAR_bams/${2%_1.fastq}.  
