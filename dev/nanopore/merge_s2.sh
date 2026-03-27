#!/bin/bash
#SBATCH -c48
#SBATCH --mem=64g
#SBATCH --gres=lscratch:200
#SBATCH --time=6:0:0

set -e
date
module load sambamba samtools/1.13

mkdir -p bam_s2 bam

find bam_s1 -name "*.bam" | parallel samtools view -h {} | awk 'length($10) > 60000 || $1 ~ /^@/' | samtools view -bS - > bam_s2/{/} 

for file in bam_s1/D15.xa{a..j}.bam; do 
filename=$(basename $file)
echo $filename
samtools view -h --threads 4 $file | awk 'length($10) > 60000 || $1 ~ /^@/' | samtools view -bS - > bam_s2/${filename%.bam}.60k.bam
sambamba index -t 4 bam_s2/${filename%.bam}.60k.bam bam_s2/${filename%.bam}.60k.bam.bai
date
done

bam_file=""
for file in bam_s2/*.bam; do bam_file+=" $file"; done
sambamba merge -t 48 bam/D15.60k.bam $bam_file
date
#rm -r bam_s2
bam_file=""
for file in bam_s1/*.bam; do bam_file+=" $file"; done
sambamba merge -t 48 bam/D15.bam $bam_file
date
#rm -r bam_s1
module load mosdepth
mkdir -p coverage
cd coverage
mosdepth -t 10 --no-per-base --by /data/OGL/resources/bed/OGLv1.GRCh38.sorted.bed --use-median --mapq 0 --fast-mode --thresholds 10,20,30 D15.md ../bam/D15.bam
mosdepth -t 10 --no-per-base --by /data/OGL/resources/bed/OGLv1.GRCh38.sorted.bed --use-median --mapq 0 --fast-mode --thresholds 10,20,30 D15.60k.md ../bam/D15.60k.bam
cd ..
samtools view --min-MQ 5 -h bam/D15.bam chrX:153974198-154559350 | awk 'length($10) > 10000 || $1 ~ /^@/' | samtools view -bS - > bam/MQ5.10k.bam
sambamba index bam/MQ5.10k.bam bam/MQ5.10k.bam.bai
samtools view -h bam/MQ5.10k.bam | awk 'length($10) > 60000 || $1 ~ /^@/' | samtools view -bS - > bam/MQ5.60k.bam
sambamba index bam/MQ5.60k.bam bam/MQ5.60k.bam.bai
