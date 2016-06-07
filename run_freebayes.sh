#!/bin/bash

module load freebayes/1.0.2
vcf_file=$1
genome=$2
bam_files=$3
region=$4
# chr:start-stop (0 base)

#usage: sbatch --mem=20g run_freebayes.sh output.vcf 1000G_hg19 '-b bam1.bam -b bam2.bam -b bam3.bam'

if [ "$genome" == "GRCh38" ]; then
	freebayes -v $vcf_file -f /data/mcgaugheyd/genomes/GRCh38/hs38DH.fa $bam_files
elif [ "$genome" == "1000G_hg19" ]; then
	freebayes  -f /fdb/GATK_resource_bundle/b37-2.8/human_g1k_v37_decoy.fasta --bam-list $bam_files > $vcf_file
else
	echo "Pick either GRCh38 or 1000G_hg19 genomes"
fi
