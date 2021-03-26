#!/bin/bash

input=$1

mkdir -p fastq

java -Xmx12g -Dsamjdk.reference_fasta=$2 -jar ~/git/bazam/build/libs/bazam.jar -r1 fastq/${input%.cram}_R1_001.fastq.gz -r2 fastq/${input%.cram}_R2_001.fastq.gz -bam $input
