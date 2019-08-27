#!/bin/bash

#May need more mem than 8gb. Took a few minutes, Coloboma panel needs 32gb with samtools.

bam=$1
no_of_regions=$2
path_filename=$3 # such as /data/OGVFB/OGL_NGS/bed/freebayes.OGLv1.500.region

module load freebayes
module load samtools

samtools depth $bam | coverage_to_regions.py /data/OGVFB/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta.fai $no_of_regions > $3


##Alternative

#module load bamtools

#bamtools coverage -in $bam | coverage_to_regions.py /data/OGVFB/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta.fai $no_of_regions > $3
#
