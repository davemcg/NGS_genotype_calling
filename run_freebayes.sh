#!/bin/bash

module load freebayes/0.9.21

vcf_file=$1
genome=$2
bam_files=$4
region=$3
# chr:start-stop (0 base)

#usage: sbatch --mem=20g run_freebayes.sh output.vcf hg19 '-b bam1.bam -b bam2.bam -b bam3.bam'

if [ "$genome" == "GRCh38" ]; then
	freebayes -v $vcf_file -f /data/mcgaugheyd/genomes/GRCh38/hs38DH.fa $bam_files
elif [ "$genome" == "hg19" ]; then
	freebayes -v $vcf_file -f /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa $bam_files
else
	echo "Pick either GRCh38 or hg19 genomes"
fi
