#!/bin/bash

module load lumpy/0.2.11

bam=$1

lumpyexpress -B $1 -o ${1%.bam}.lumpy0.2.11.vcf
