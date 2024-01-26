#!/bin/bash

FASTQ=$1
sample=$2

module load minimap2 sambamba
minimap2 -a -x map-ont -Y -t 3 -R "@RG\tLB:plasmidExpress\tID:plEx-$sample\tPL:ONT\tSM:$sample" /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.mmi $FASTQ \
	| sambamba sort -u --compression-level 6 --tmpdir=/lscratch/$SLURM_JOB_ID -t 3 -o $sample.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin)
