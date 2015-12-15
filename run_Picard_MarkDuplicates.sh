#!/bin/bash

module load picard

bam_name=$1


java -Xmx40g -jar $PICARDJARPATH/picard.jar MarkDuplicates\
	INPUT=$bam_name \
	OUTPUT=${bam_name%.bam}.dupsMarked.bam \
	METRICS_FILE=${bam_name%.bam}.dupsMarked.metrics
