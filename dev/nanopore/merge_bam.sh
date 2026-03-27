#!/bin/bash
#SBATCH -c48
#SBATCH --mem=64g
#SBATCH --gres=lscratch:200
#SBATCH --time=0:30:0

# merge upto 500 files.
#ls temp_bam/*.bam > bam.list.txt
#split -l 500 bam.list.txt
#for file in xa{a..j}; do sbatch merge_s1.sh $file; done

module load sambamba
mkdir -p bam
bam_file=""
while read -r line; do bam_file+=" $line"; done < $1

sambamba merge -t 48 bam/$2.bam $bam_file

#get the reads in specific regions only
#samtools view -b --threads 6 rico.hom.bam chr6:124700000-124800000 chr13:74700000-81130000 -o Eno2.chr13.bam --write-index

