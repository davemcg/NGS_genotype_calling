#!/bin/bash

# Converts from numeric chromosomes (e.g. 1,2,3) to named (e.g. chr1, chr2, chr3)
# REMOVES NON-CANONICAL chromosomes (e.g. Contig, JG, Un)
# Author: David McGaughey
# Modified slightly from http://seqanswers.com/forums/showthread.php?t=22504
# usage: convert_chrNames.sh INPUT.bam OUTPUT.bam 

module load samtools
INPUT=$1
OUTPUT=$2

samtools view -h $INPUT | awk 'BEGIN{FS=OFS="\t"} (/^@/ && !/@SQ/){print $0} $2~/^SN:[1-9]|^SN:X|^SN:Y|^SN:MT/{print $0}  $3~/^[1-9]|X|Y|MT/{$3="chr"$3; print $0} ' | sed 's/SN:/SN:chr/g' | sed 's/chrMT/chrM/g' | samtools view -bS - > $OUTPUT
