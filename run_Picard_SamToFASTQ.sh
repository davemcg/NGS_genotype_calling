#!/bin/bash

module load picard

java -Xmx8g -jar $PICARDJARPATH/picard.jar SamToFastq \
	I=$1 \
	FASTQ=${1%.bam}_1.fastq.gz \
	SECOND_END_FASTQ=${1%.bam}_2.fastq.gz
