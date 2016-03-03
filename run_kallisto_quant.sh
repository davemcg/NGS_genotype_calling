#!/bin/bash

module load kallisto/0.42.4

index=$1
output_folder=$2
fastq1=$3
fastq2=$4


kallisto quant -b 100 -i $index -o $output_folder -t $SLURM_CPUS_PER_TASK $3 $4

