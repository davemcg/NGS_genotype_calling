#!/bin/bash
mkdir -p fastq
while IFS= read -r line; do
	sm=$(echo $line | cut -d ' ' -f 1);
	bam=$(echo $line | cut -d ' ' -f 3);
	sbatch --mem=12g ~/git/NGS_genotype_calling/bam_to_fastq/bazam_simple.sh $sm $bam
done < RAWDATA/map_file.txt
