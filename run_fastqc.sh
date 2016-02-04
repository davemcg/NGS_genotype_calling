#!/bin/bash

module load fastqc
input=$1
mkdir fastqc

fastqc -o fastqc -f bam $input



