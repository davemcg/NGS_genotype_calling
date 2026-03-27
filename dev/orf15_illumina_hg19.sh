#!/bin/bash
#SBATCH -c8
#SBATCH --mem=64g
#SBATCH --gres=lscratch:100
#SBATCH --time=2:0:0
set -e
module load clair3/1.0.10 annovar/2020-06-08

mkdir -p clair3 clairAnnotation

for file in old_bam/*.bam; do
sample=$(basename $file | sed "s/.bam//")
clair3 \
  --bam_fn=old_bam/$sample.bam \
  --ref_fn=/data/OGL/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta \
  --threads=8 \
  --platform=ilmn \
  --bed_fn=/data/OGL/resources/bed/RPGR_ORF15_hg19.bed \
  --model_path=/data/OGL/resources/clair3/ilmn \
  --ref_pct_full=1.0 --var_pct_full=0.5 \
  --chunk_size=-1 \
  --output=/lscratch/$SLURM_JOB_ID
cp /lscratch/$SLURM_JOB_ID/merge_output.vcf.gz clair3/$sample.vcf.gz
cp /lscratch/$SLURM_JOB_ID/merge_output.vcf.gz.tbi clair3/$sample.vcf.gz.tbi
convert2annovar.pl -format vcf4old clair3/$sample.vcf.gz -includeinfo --outfile clairAnnotation/$sample.avinput
ver=hg19
table_annovar.pl clairAnnotation/$sample.avinput \
$ANNOVAR_DATA/$ver \
-buildver $ver \
-remove \
-out clairAnnotation/$sample.avinput \
--protocol refGeneWithVer \
-operation g \
--argument '-hgvs' \
--polish -nastring . \
--thread 1 \
--otherinfo
sed "1 s/Otherinfo1\tOtherinfo2\tOtherinfo3\tOtherinfo4\tOtherinfo5\tOtherinfo6\tOtherinfo7\tOtherinfo8\tOtherinfo9\tOtherinfo10/CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGT_FIELDS/" clairAnnotation/$sample.avinput."$ver"_multianno.txt
cut -f 6- clairAnnotation/$sample.avinput."$ver"_multianno.txt > clairAnnotation/$sample.annovar.tsv
rm clairAnnotation/$sample.avinput clairAnnotation/$sample.avinput."$ver"_multianno.txt
done