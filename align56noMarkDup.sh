#!/bin/bash
#SBATCH --gres=lscratch:400
#SBATCH --cpus-per-task=56
#SBATCH --mem=100g
#SBATCH --time=8:0:0

export TMPDIR=/lscratch/$SLURM_JOB_ID
echo @RG\\tID:HTVF5DSXY.4_D939_001\\tSM:D939_001\\tLB:novogene20_D939_001\\tPL:ILLUMINA
mkdir -p lane_bam
THREADS=$SLURM_CPUS_PER_TASK
module load bwa-mem2/2-2.1 samblaster/0.1.25 sambamba/0.8.0
bwa-mem2 mem -t $((THREADS - 6)) -K 100000000 -M -Y -B 4 -O 6 -E 1 -R @RG\\tID:HTVF5DSXY.4_D939_001\\tSM:D939_001\\tLB:novogene20_D939_001\\tPL:ILLUMINA /data/OGL/resources/1000G_phase2_GRCh37/bwa-mem2/human_g1k_v37_decoy.fasta fastq/D939-001_FDSW202654548-1r_HTVF5DSXY_L4_1.fq.gz fastq/D939-001_FDSW202654548-1r_HTVF5DSXY_L4_2.fq.gz | samblaster -M --acceptDupMarks --addMateTags --quiet | sambamba sort -u --tmpdir=/lscratch/$SLURM_JOB_ID -t 50 -o lane_bam/D939-001_FDSW202654548-1r_HTVF5DSXY_L4.bam <(sambamba view -S -f bam -l 0 -t 50 /dev/stdin)