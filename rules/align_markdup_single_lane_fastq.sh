#!/bin/bash
#SBATCH --gres=lscratch:800
#SBATCH --cpus-per-task=56
#SBATCH --mem=100g
#SBATCH --partition=norm,quick
#SBATCH --time=4:0:0

config=$1
sample=$2

set -e
mkdir -p sample_bam

module load $(grep "^bwa-mem2_version:" $config | head -n 1 | cut -d"'" -f 2)
module load $(grep "^samblaster_version:" $config | head -n 1 | cut -d"'" -f 2)
module load $(grep "^sambamba_version:" $config | head -n 1 | cut -d"'" -f 2)

METADATA_FILE=$(grep "^metadata_file:" $config | head -n 1 | cut -d"'" -f 2)
RG=$(grep "$sample" $METADATA_FILE | head -n 1 | cut -d, -f 3)
echo $RG
FASTQ=$(grep "$sample" $METADATA_FILE | cut -d, -f 2 | sed 's#^#fastq/#'| tr '\n' ' ' )

bwa-mem2 mem -t $(($SLURM_CPUS_PER_TASK-2)) -K 100000000 -M -Y -B 4 -O 6 -E 1 -R $RG \
	$(grep "^bwa-mem2_ref:" $config | head -n 1 | cut -d"'" -f 2) $FASTQ \
	| samblaster -M --addMateTags --quiet \
	| sambamba sort -u --tmpdir=/lscratch/$SLURM_JOB_ID -t $(($SLURM_CPUS_PER_TASK-2)) -o sample_bam/$sample.markDup.bam \
	<(sambamba view -S -f bam -l 0 -t $(($SLURM_CPUS_PER_TASK-2)) /dev/stdin)
