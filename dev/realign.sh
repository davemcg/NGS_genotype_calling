#!/bin/bash
#SBATCH --gres=lscratch:500
#SBATCH --cpus-per-task=56
#SBATCH --mem=128g
#SBATCH --partition=norm
#SBATCH --time=12:0:0


#Previously: 6h:16min

export TMPDIR=/lscratch/$SLURM_JOB_ID
module load bazam/1.0.1
module load bwa/0.7.17 samblaster/0.1.25 sambamba/0.8.1
BAMFILE=old_bam/G04747.bam

java -Xmx112g -jar $BAZAMPATH/bazam.jar -n 54 -bam old_bam/G04747.bam | bwa mem -t $((56-2)) -K 100000000 -M -Y -B 4 -O 6 -E 1 -p -R @RG\\tID:HK3GHDSX3.2.19274584\\tLB:19274584\\tPL:ILLUMINA\\tSM:G04747\\tPU:HK3GHDSX3.2.19274584 /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa - | samblaster -M --addMateTags --quiet | sambamba sort -u --tmpdir=/lscratch/$SLURM_JOB_ID -t $((56-2)) -o sample_bam/G04747.markDup.bam <(sambamba view -S -f bam -l 0 -t $((56-2)) /dev/stdin)
mv sample_bam/G04747.markDup.bam.bai sample_bam/G04747.markDup.bai