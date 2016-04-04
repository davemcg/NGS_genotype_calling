#!/bin/bash

module load picard

mkdir ${1%.bam}_fastq
java -Xmx20g -jar -jar $PICARDJARPATH/picard.jar SamToFastq \
	I=$1 \
	OUTPUT_PER_RG="true" \
	OUTPUT_DIR=${1%.bam}_fastq
