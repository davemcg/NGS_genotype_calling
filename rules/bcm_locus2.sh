#!/bin/bash
#SBATCH --gres=lscratch:100
#SBATCH --cpus-per-task=32
#SBATCH --mem=64g
#SBATCH --partition=norm
#SBATCH --time=16:0:0

config=$1
sample=$2

#selecting mapped reads using subset_mapped_reads.sh reduced sample processing time to 44 hours.
#selecting reads using subset_chr1-M_reads.sh: 40 hours (without P flag. Adding P flag makes the run > 48h).

set -e
mkdir -p sample_bam/subset
module load $(grep "^samtools_version:" $config | head -n 1 | cut -d"'" -f 2)

export TMPDIR=/lscratch/$SLURM_JOB_ID
module load bazam/1.0.1 bwa/0.7.17 samblaster/0.1.26 sambamba/1.0.0
RG=$(samtools view -H sample_bam/$sample.markDup.bam | grep "^@RG" | head -n 1 | sed 's/\t/\\t/g')
java -Xmx32g -jar $BAZAMPATH/bazam.jar -bam sample_bam/$sample.markDup.bam --regions chrX:153929000-154373500 | bwa mem -t 16 -K 100000000 -M -Y -B 4 -O 6 -E 1 -p -R $RG /data/OGL/resources/genomes/GRCh38/GRCh38Decoy2.shortread.bcm.fa - | samblaster -M --addMateTags --quiet | sambamba sort -u --tmpdir=/lscratch/$SLURM_JOB_ID -t 16 -o bcmlocus/bam/$sample.2.bcm.bam <(sambamba view -S -f bam -l 0 -t 16 /dev/stdin)
