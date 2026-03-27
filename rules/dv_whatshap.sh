#!/bin/bash
#SBATCH --gres=lscratch:200
#SBATCH --cpus-per-task=26
#SBATCH --mem=72g
#SBATCH --partition=norm,quick
#SBATCH --time=4:0:0

config=$1
sample=$2

module load $(grep "^samtools_version:" $config | head -n 1 | cut -d"'" -f 2)
module load $(grep "^whatshap_version:" $config | head -n 1 | cut -d"'" -f 2)
module load parallel
WORK_DIR=/lscratch/$SLURM_JOB_ID
FILTEREDVCF=/lscratch/$SLURM_JOB_ID/$sample.vcf.gz
rm -rf $WORK_DIR/*
cp deepvariant/$(grep "^analysis_batch_name:" $config | head -n 1 | cut -d"'" -f 2).glnexus.vcf.gz* /lscratch/$SLURM_JOB_ID
bcftools view --threads $SLURM_CPUS_PER_TASK -Oz --samples $sample $WORK_DIR/$(grep "^analysis_batch_name:" $config | head -n 1 | cut -d"'" -f 2).glnexus.vcf.gz \
	-o $WORK_DIR/$sample.vcf.gz && rm $WORK_DIR/$(grep "^analysis_batch_name:" $config | head -n 1 | cut -d"'" -f 2).glnexus.vcf.gz*
tabix -f -p vcf $WORK_DIR/$sample.vcf.gz
CONTIGFILE="/data/OGL/resources/whatshap/vcf.contig.filename.$(grep "^genomeBuild:" $config | head -n 1 | cut -d"'" -f 2).txt"
cp $(grep "^ref_genome:" $config | head -n 1 | cut -d"'" -f 2) $(grep "^ref_genome:" $config | head -n 1 | cut -d"'" -f 2).fai $WORK_DIR
REF_GENOME=$WORK_DIR/$(basename $(grep "^ref_genome:" $config | head -n 1 | cut -d"'" -f 2))
cp sample_bam/$sample.markDup.bam sample_bam/$sample.markDup.bam.bai $WORK_DIR
mkdir -p /lscratch/$SLURM_JOB_ID/filtered
mkdir -p /lscratch/$SLURM_JOB_ID/phased
( cat $CONTIGFILE | parallel -C "\t" -j 21 "bcftools filter --threads $(($SLURM_CPUS_PER_TASK-6)) -r {1} --output-type z $FILTEREDVCF -o $WORK_DIR/filtered/{2}.filtered.vcf.gz" ) && echo "Filtered vcf split to chr" || exit 5
 cat $CONTIGFILE | parallel -C "\t" -j 21 "tabix -f -p vcf $WORK_DIR/filtered/{2}.filtered.vcf.gz" ) && echo "Chr vcf index created" || exit 6
( cat $CONTIGFILE | parallel -C "\t" -j 21 --tmpdir $WORK_DIR --eta --halt 2 --line-buffer \
 --tag "whatshap phase --reference $REF_GENOME \
--indels --ignore-read-groups $WORK_DIR/filtered/{2}.filtered.vcf.gz $WORK_DIR/$sample.markDup.bam \
| bgzip -f --threads $(($SLURM_CPUS_PER_TASK-10)) > $WORK_DIR/phased/{2}.phased.vcf.gz" \
) && echo "whatshap on chr completed" || exit 7
( cat $CONTIGFILE | parallel -C "\t" -j 21 "tabix -f -p vcf $WORK_DIR/phased/{2}.phased.vcf.gz" ) && echo "Phased-chr vcf index created" || exit 8
PHASEDCHRFILE=""
cut -f 2 $CONTIGFILE > $WORK_DIR/temp.chr.txt
while read line; do PHASEDCHRFILE+=" /lscratch/$SLURM_JOB_ID/phased/$line.phased.vcf.gz"; done < $WORK_DIR/temp.chr.txt
echo "chr files are $PHASEDCHRFILE"
bcftools concat --threads $SLURM_CPUS_PER_TASK --output-type z $PHASEDCHRFILE > deepvariant/vcf/$sample.dv.glnexus.phased.vcf.gz || exit 8
tabix -f -p vcf deepvariant/vcf/$sample.dv.glnexus.phased.vcf.gz
