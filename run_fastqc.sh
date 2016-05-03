#!/bin/bash

module load fastqc
input=$1
mkdir fastqc

fastqc -o fastqc -t $SLURM_CPUS_PER_TASK  $input



