#!/bin/bash

module load picard

# Sort bam, then marks duplicates, then creates index

echo $1

java -Xmx20g -jar $PICARDJARPATH/picard.jar SortSam \
	INPUT=$1 OUTPUT=${1%.bam}.sorted.bam SORT_ORDER=coordinate

sorted_bam=${1%.bam}.sorted.bam

java -Xmx20g -jar $PICARDJARPATH/picard.jar MarkDuplicates \
	INPUT=$sorted_bam OUTPUT=${sorted_bam%.bam}.markDup.bam METRICS_FILE=${sorted_bam%.bam}.markDup.metrics

sorted_markDup_bam=${sorted_bam%.bam}.markDup.bam

java -Xmx20g -jar $PICARDJARPATH/picard.jar BuildBamIndex \
	INPUT=$sorted_markDup_bam OUTPUT=$sorted_markDup_bam.bai

# Delete intermediate files
rm $1
rm ${1%.bam}.sorted.bam
#rm ${sorted_bam%.bam}.markDup.bam
