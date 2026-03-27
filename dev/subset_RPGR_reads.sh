#!/bin/bash
#SBATCH -c8
#SBATCH --mem=18g
#SBATCH --gres=lscratch:100
#SBATCH --time=4:0:0

input_bam_dir=$1 #do not end with "/"
output_bam_dir=$2 #do not end with "/"

module load samtools/1.21
#RPGR.bed  also in git/dev, slightly larger than RPGR gene, 60kb, chrX:38,268,201-38,328,200
mkdir -p $output_bam_dir

for bamfile in $input_bam_dir/*.bam; do 
filename=$(basename $bamfile)
echo "subsetting" $filename
samtools view --threads $SLURM_CPUS_PER_TASK -bP --target /data/OGL/resources/bed/RPGR.bed --output $output_bam_dir/${filename%.bam}.RPGR.bam $bamfile
samtools index -@ $SLURM_CPUS_PER_TASK $output_bam_dir/${filename%.bam}.RPGR.bam
done

# for sample in G01877 G02164 G03208 G03377 G04033 G05456 G05462; do 
# echo $sample
# samtools view --threads $SLURM_CPUS_PER_TASK -bP --target /data/OGL/resources/bed/RPGR_ORF15.bed --output /data/OGL/genome/orf15_mini_bam/negative/$sample.bam $sample.bam
# samtools index -@ $SLURM_CPUS_PER_TASK /data/OGL/genome/orf15_mini_bam/negative/$sample.bam
# done

# for sample in D1163_001 D1280_01 D1287_03 D1323_01 D1620_01 D1631_01 D2030_02; do 
# echo $sample
# samtools view --threads $SLURM_CPUS_PER_TASK -bP --target /data/OGL/resources/bed/RPGR_ORF15.bed --output /data/OGL/genome/orf15_mini_bam/negative/$sample.bam bam/$sample.bam
# samtools index -@ $SLURM_CPUS_PER_TASK /data/OGL/genome/orf15_mini_bam/negative/$sample.bam
# done

# for sample in RY1P RY2P RY2; do 
# echo $sample
# samtools view --threads $SLURM_CPUS_PER_TASK -bP --target /data/OGL/resources/bed/RPGR_ORF15.bed --output /data/OGL/genome/orf15_mini_bam/negative/$sample.bam bam/$sample.bam
# samtools index -@ $SLURM_CPUS_PER_TASK /data/OGL/genome/orf15_mini_bam/negative/$sample.bam
# done

# for sample in RY1P RY2P RY2; do 
# echo $sample
# samtools view --threads $SLURM_CPUS_PER_TASK -bP --target /data/OGL/resources/bed/RPGR_ORF15.bed --output /data/OGL/genome/orf15_mini_bam/negative/$sample.bam bam/$sample.bam
# samtools index -@ $SLURM_CPUS_PER_TASK /data/OGL/genome/orf15_mini_bam/negative/$sample.bam
# done

# for sample in D2173x03_REV24A0010319; do 
# echo $sample
# samtools view --threads $SLURM_CPUS_PER_TASK -bP --target /data/OGL/resources/bed/RPGR_ORF15.bed --output /data/OGL/genome/orf15_mini_bam/negative/$sample.bam bam/$sample.bam
# samtools index -@ $SLURM_CPUS_PER_TASK /data/OGL/genome/orf15_mini_bam/negative/$sample.bam
# done