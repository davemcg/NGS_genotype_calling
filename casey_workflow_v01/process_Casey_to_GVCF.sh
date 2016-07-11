#!/bin/bash

# Run in folder containing subfolders of fastq

folder=$1
panel=$2
bed=$3

fastqs=$(ls $folder/*fastq.gz)
#splits fastqs into $1 and $2 by space
set $fastqs


############
# Alignment
j1=$( sbatch --cpus-per-task 10 ~/bin/casey_workflow_v01/align_with_bwa.py $1 $2 $folder $panel )
############

output_bam=$folder'/'$folder'_'$panel'.bwa-mem.b37.bam'

############
# BAM processing and GVCF calling
sbatch --mem=20G --time=16:00:00 --dependency=afterok:$j1 ~/bin/casey_workflow_v01/process_and_callGVCF.sh $output_bam $bed
############


