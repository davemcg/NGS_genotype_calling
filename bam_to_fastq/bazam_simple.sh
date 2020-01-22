#!/bin/bash
sm=$1
bam=$2
java -Xmx12g -jar ~/git/bazam/build/libs/bazam.jar -r1 fastq/"$sm"_R1_001.fastq.gz -r2 fastq/"$sm"_R2_001.fastq.gz -bam RAWDATA/$bam