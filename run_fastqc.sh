#!/bin/bash

module load fastqc


input=$1
echo $PWD
mkdir $PWD/fastqc
fastqc -o $PWD/fastqc -f bam $input
