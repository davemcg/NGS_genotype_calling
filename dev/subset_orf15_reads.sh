#!/bin/bash
#SBATCH -c8
#SBATCH --mem=18g
#SBATCH --gres=lscratch:100
#SBATCH --time=4:0:0

input=$1
output=$2
module load samtools/1.21
samtools view --threads $SLURM_CPUS_PER_TASK -b --target /data/OGL/resources/bed/RPGR_ORF15.bed --output $output $input
samtools index -@ $SLURM_CPUS_PER_TASK $output


for bamfile in bam/*.bam; do 
filename=$(basename $bamfile)
echo $filename
samtools view --threads $SLURM_CPUS_PER_TASK -bP --target /data/OGL/resources/bed/RPGR_ORF15.bed --output /data/OGL/genome/orf15_mini_bam/positive/$filename $bamfile
samtools index -@ $SLURM_CPUS_PER_TASK /data/OGL/genome/orf15_mini_bam/positive/$filename
done

for sample in G01877 G02164 G03208 G03377 G04033 G05456 G05462; do 
echo $sample
samtools view --threads $SLURM_CPUS_PER_TASK -bP --target /data/OGL/resources/bed/RPGR_ORF15.bed --output /data/OGL/genome/orf15_mini_bam/negative/$sample.bam $sample.bam
samtools index -@ $SLURM_CPUS_PER_TASK /data/OGL/genome/orf15_mini_bam/negative/$sample.bam
done

for sample in D1163_001 D1280_01 D1287_03 D1323_01 D1620_01 D1631_01 D2030_02; do 
echo $sample
samtools view --threads $SLURM_CPUS_PER_TASK -bP --target /data/OGL/resources/bed/RPGR_ORF15.bed --output /data/OGL/genome/orf15_mini_bam/negative/$sample.bam bam/$sample.bam
samtools index -@ $SLURM_CPUS_PER_TASK /data/OGL/genome/orf15_mini_bam/negative/$sample.bam
done

for sample in RY1P RY2P RY2; do 
echo $sample
samtools view --threads $SLURM_CPUS_PER_TASK -bP --target /data/OGL/resources/bed/RPGR_ORF15.bed --output /data/OGL/genome/orf15_mini_bam/negative/$sample.bam bam/$sample.bam
samtools index -@ $SLURM_CPUS_PER_TASK /data/OGL/genome/orf15_mini_bam/negative/$sample.bam
done

for sample in RY1P RY2P RY2; do 
echo $sample
samtools view --threads $SLURM_CPUS_PER_TASK -bP --target /data/OGL/resources/bed/RPGR_ORF15.bed --output /data/OGL/genome/orf15_mini_bam/negative/$sample.bam bam/$sample.bam
samtools index -@ $SLURM_CPUS_PER_TASK /data/OGL/genome/orf15_mini_bam/negative/$sample.bam
done

for sample in D2173x03_REV24A0010319; do 
echo $sample
samtools view --threads $SLURM_CPUS_PER_TASK -bP --target /data/OGL/resources/bed/RPGR_ORF15.bed --output /data/OGL/genome/orf15_mini_bam/negative/$sample.bam bam/$sample.bam
samtools index -@ $SLURM_CPUS_PER_TASK /data/OGL/genome/orf15_mini_bam/negative/$sample.bam
done