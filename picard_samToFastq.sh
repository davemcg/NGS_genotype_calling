#!/bin/bash

module load picard

java -Xmx20g -jar -jar $PICARDJARPATH/picard.jar SamToFastq \
	I=$1 \
	FASTQ=${1%.bam}_1.fastq \
	SECOND_END_FASTQ=${1%.bam}_2.fastq
