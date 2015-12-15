#!/bin/bash

module load bedtools

bam=$1
output=$2


genomeCoverageBed -bg -ibam $bam -g /data/mcgaugheyd/genomes/hg19/hg19.genome > /tmp/$bam.temp

sort -k1,1 -k2,2n /tmp/$bam.temp > /tmp/$bam.sorted.temp

rm /tmp/$bam.temp

bedGraphToBigWig /tmp/$bam.sorted.temp /data/mcgaugheyd/genomes/hg19/hg19.noheader.genome $output

rm /tmp/$bam.sorted.temp
