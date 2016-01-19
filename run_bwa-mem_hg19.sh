#!/bin/bash

module load bwa/0.7.12
module load samtools/1.2

fastq1=$1
fastq2=$2
rg=$3
output_name=$4

bwa mem -M -t 10 -B 4 -O 6 -E 1 -M -R $rg \
	 /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa \
	$1 $2 | samtools view -1 - > $4

rm $1
rm $2
