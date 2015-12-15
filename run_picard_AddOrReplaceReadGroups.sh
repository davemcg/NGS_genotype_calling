#!/bin/bash

module load picard

bam_name=$1


java -Xmx4g -jar $PICARDJARPATH/picard.jar AddOrReplaceReadGroups \
	INPUT=$bam_name \
	OUTPUT=${bam_name%.bam}.RG.bam \
	$2
#	RGID=$3 \
#	RGLB=$4 \
#	RGPL=illumina \
#	RGSM=$5 \
#	RGPU=$6 
