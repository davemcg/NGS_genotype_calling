#!/bin/bash
#SBATCH -c48
#SBATCH --mem=64g
#SBATCH --gres=lscratch:200
#SBATCH --time=0:30:0

# merge upto 500 files.
#ls fastq/*.bam > bam.list.txt
#split -l 500 bam.list.txt
#for file in xa{b..i}; do sbatch merge_test.sh $file; done

module load sambamba

bam_file=""
while read -r line; do bam_file+=" $line"; done < $1

sambamba merge -t 48 D15.$1.bam $bam_file
