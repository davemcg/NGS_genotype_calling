#!/bin/bash

#https://github.com/ekg/freebayes/blob/master/scripts/freebayes-parallel

module load freebayes/1.3.1
fasta_generate_regions.py /data/OGVFB/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta.fai 100000 > v37d5.100kb.region

for chr in {1..22} X Y; do grep ^$chr: v37d5.100kb.region > region.$chr; done
grep -e "^MT" -e "^GL" -e "^NC" v37d5.100kb.region > region.MT_contigs

module load bamtools

bamtools coverage -in aln.bam | coverage_to_regions.py ref.fa 500 >ref.fa.500.regions