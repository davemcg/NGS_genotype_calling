#!/bin/bash
#SBATCH --gres=lscratch:200
#SBATCH --cpus-per-task=28
#SBATCH --mem=64g
#SBATCH --partition=norm
#SBATCH --time=6:0:0

module load samblaster/0.1.26 samtools/1.23 
samtools collate -@ 28 -O --output-fmt SAM sample_bam/D1427x01.markDup.bam \
| samblaster --acceptDupMarks --addMateTags --quiet \
| samtools view -@ 28 -Sb - \
| samtools sort -@ 28 -O BAM -o sample_bam/D1427x01mtag.markDup.bam -

samtools index -@ 28 sample_bam/D1427x01mtag.markDup.bam 
#--output sample_bam/D1427x01mtag.markDup.bam.bai

#34048608: 33 cpu, 46g mem, 52 min
#clinSV: D1427x01 34055815
#		D1427x01mtag 34055823
#mtag reduced the lumpy step from 4 to 3 hours. Total time: 9 vs 8 hours, this step also reduced bam file sizes from 