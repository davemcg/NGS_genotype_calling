#!/bin/bash
#SBATCH --gres=lscratch:150
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH --partition=quick
#SBATCH --time=4:0:0

export TMPDIR=/lscratch/$SLURM_JOB_ID
module load samtools/1.13 bazam/1.0.1 bwa/0.7.17 samblaster/0.1.25 sambamba/0.8.1

RG=$(samtools view -H sample_bam/G05089.markDup.bam | grep "^@RG" | head -n 1 | sed 's/\t/\\\\t/g')
java -Xmx32g -jar $BAZAMPATH/bazam.jar -n 16 -bam sample_bam/G05089.markDup.bam --regions chrX:153929000-154373500 | bwa mem -t 16 -K 100000000 -M -Y -B 4 -O 6 -E 1 -p -R $RG /data/OGL/resources/genomes/GRCh38/GRCh38Decoy2.fa - | samblaster -M --addMateTags --quiet | sambamba sort -u --tmpdir=/lscratch/$SLURM_JOB_ID -t 16 -o bam/bcmlocus/G05089.bcm.bam <(sambamba view -S -f bam -l 0 -t 16 /dev/stdin)