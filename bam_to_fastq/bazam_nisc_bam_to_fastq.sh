#!/bin/bash

#for file in RAWDATA/*.bam; do sbatch --mem=12g ~/git/NGS_genotype_calling/bam_to_fastq/bazam_nisc_bam_to_fastq.sh $file; done

original_bam=$1

mkdir -p fastq

module load samtools/1.9

samtools index $original_bam 

mv $original_bam.bai ${original_bam%.bam}.bai

oldrg=$(samtools view -H $original_bam | grep '^@RG' )
echo $oldrg
id=$(echo $oldrg | cut -d\  -f 2 | cut -d: -f 2)
echo $id
sm=$(echo $oldrg | cut -d\  -f 3 | cut -d: -f 2)
echo $sm

#sm=$(echo $1 | awk -F/ '{print $NF}' | cut -d\. -f 1)
#echo $sm

java -Xmx12g -jar ~/git/bazam/build/libs/bazam.jar -r1 fastq/"$sm"_"$id"_R1_001.fastq.gz -r2 fastq/"$sm"_"$id"_R2_001.fastq.gz -bam $original_bam

