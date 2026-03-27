#!/bin/bash

#May need more mem than 8gb. Took a few minutes, Coloboma panel needs 32gb with samtools.

bam=$1
ref=$2
no_of_regions=$3
path_filename=$4 # such as /data/OGVFB/OGL_NGS/bed/freebayes.OGLv1.500.region

module load freebayes
module load samtools

samtools depth $bam | coverage_to_regions.py $2.fai $no_of_regions > $3

#samtools depth bam/D565_001.bam | coverage_to_regions.py /data/OGL/resources/genomes/GRCh38/GRCh38Decoy2.fa.fai 10000 > /data/OGL/resources/freebayesRegion/xGenV1.10000.region

#samtools depth bam/2226.bam | coverage_to_regions.py /data/OGL/resources/genomes/GRCh38/GRCh38Decoy2.fa.fai 500 > /data/OGL/resources/freebayesRegion/freebayes.OGLv1.GRCh38.500.region


##Alternative

#module load bamtools

#bamtools coverage -in $bam | coverage_to_regions.py /data/OGVFB/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta.fai $no_of_regions > $3
#
#wgs

# mkdir -p GRCh38
# fasta_generate_regions.py /data/OGL/resources/genomes/GRCh38/GRCh38Decoy2.fa.fai 100000 > GRCh38/GRCh38.100kb.region
# for chr in chr{1..22} chrX chrY; do grep "$chr:" GRCh38.100kb.region > region.$chr; done

# #line no. of the last chrY, after check tail region.chrY
# sed -n '/chrY:57200000-57227415/=' GRCh38.100kb.region
# #showed 30894
# tail -n +30895 GRCh38.100kb.region > region.MT_contigs
# #Contigs have many small regions without chr thus
# tail -n +30895 GRCh38.100kb.region | grep '^chr' > region.MT_contigs