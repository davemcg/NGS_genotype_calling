#!/bin/bash

module load QoRTs/1.1.6


mkdir qorts

bam=$1
gtf=$2

java -Xmx20G -jar $QORTS_JARFILE QC \
	$bam \
	$gtf \
	qorts/${1%.bam}	
