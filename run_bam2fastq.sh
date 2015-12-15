#!/bin/bash

module load bam2fastq/1.1.0


input_bam=$1
bam2fastq -o ${1%.bam}#.fastq $1
