#!/bin/bash

module load fastqc
input=$1
mkdir fastqc

fastqc -o fastqc -t 2 -f bam $input



