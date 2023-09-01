#!/bin/bash
#SBATCH -c48
#SBATCH --mem=64g
#SBATCH --gres=lscratch:200
#SBATCH --time=0:30:0

# merge upto 500 files.
#ls bam_temp/*.bam > bam.list.txt
#split -l 500 bam.list.txt
#for file in xa{a..j}; do sbatch merge_s1.sh $file; done

module load sambamba
mkdir -p bam
bam_file=""
while read -r line; do bam_file+=" $line"; done < $1

sambamba merge -t 48 bam/D15.$1.bam $bam_file
