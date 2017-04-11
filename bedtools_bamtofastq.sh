#!/bin/bash

module load bedtools/2.25.0

bam=$1
setup=$2 #paired or single


if [ "$setup" == "single" ]
then
	bedtools bamtofastq -i $1 -fq ${1%.bam}.fq
elif [ "$setup" == "paired" ]
then
	samtools sort -n $1 | bedtools bamtofastq -i - -fq ${1%.bam}_1.fq -fq2 ${1%.bam}_2.fq
else
	echo 'Specify paired or single bam file'
fi

