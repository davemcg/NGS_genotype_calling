#!/bin/bash
#SBATCH --gres=lscratch:200
#SBATCH --cpus-per-task=56
#SBATCH --mem=128g
#SBATCH --partition=norm
#SBATCH --time=16:0:0

set -e
module load samtools/1.21 clair3/1.1.2


#config=$1
REF=$1 #/data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa
sample=$2
LWORK_DIR=/lscratch/$SLURM_JOB_ID
# D2206x02.cram 49G; D2238x01.cram 36G
BAMFILE="cram/$sample.cram"
if [ -e $BAMFILE.cram.crai ] || [ -e ${BAMFILE%.cram}.crai ] ; then
 echo "index present"
  else
   samtools index -@ $SLURM_CPUS_PER_TASK $BAMFILE
fi
samtools view -T $REF -@ $SLURM_CPUS_PER_TASK \
 -b -o $LWORK_DIR/$sample.bam $BAMFILE
samtools index -@ $SLURM_CPUS_PER_TASK $LWORK_DIR/$sample.bam

# module load $(grep "^clair3_version:" $config | head -n 1 | cut -d"'" -f 2)
# module load $(grep "^whatshap_version:" $config | head -n 1 | cut -d"'" -f 2)
# module load $(grep "^samtools_version:" $config | head -n 1 | cut -d"'" -f 2)

clair3 --bam_fn $LWORK_DIR/$sample.bam \
 --ref_fn=$REF \
 --threads=$SLURM_CPUS_PER_TASK --platform=ilmn --gvcf \
 --model_path=/data/OGL/resources/clair3/ilmn \
 --sample_name=$sample \
 --output=/lscratch/$SLURM_JOB_ID/
 
cp $LWORK_DIR/merge_output.gvcf.gz clair3/gvcf/$sample.gvcf.gz
cp $LWORK_DIR/merge_output.gvcf.gz.tbi clair3/gvcf/$sample.gvcf.gz.tbi
cp $LWORK_DIR/phased_merge_output.vcf.gz clair3/vcf/$sample.vcf.gz
cp $LWORK_DIR/phased_merge_output.vcf.gz.tbi clair3/vcf/$sample.vcf.gz.tbi
