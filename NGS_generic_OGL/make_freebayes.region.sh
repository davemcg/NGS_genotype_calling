#!/bin/bash

#May need more mem than 8gb. Took a few minutes

bam=$1
no_of_regions=$2

module load freebayes
module load samtools

samtools depth $bam | coverage_to_regions.py /data/OGVFB/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta.fai $no_of_regions > /data/OGVFB/OGL_NGS/bed/freebayes.OGLv1.$2.region