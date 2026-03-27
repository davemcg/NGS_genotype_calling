#!/bin/bash
#SBATCH -c8
#SBATCH --mem=64g
#SBATCH --gres=lscratch:100
#SBATCH --time=2:0:0
set -e
bed=$1
module load clair3/20250303 annovar/2020-06-08 samtools/1.21

mkdir -p clair3 clairAnnotation

for file in bam/*.bam; do
sample=$(basename $file | sed "s/.bam//")
clair3 \
  --bam_fn=bam/$sample.bam \
  --ref_fn=/data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
  --threads=8 \
  --platform=ilmn \
  --bed_fn=$bed \
  --model_path=/data/OGL/resources/clair3/ilmn \
  --snp_min_af=0.05 \
  --print_ref_calls \
  --ref_pct_full=1.0 --var_pct_full=0.5 \
  --chunk_size=-1 \
  --output=/lscratch/$SLURM_JOB_ID/$sample
AddPairEndAlleleDepth.py --bam_fn bam/$sample.bam \
	--clair3_vcf_input /lscratch/$SLURM_JOB_ID/$sample/merge_output.vcf.gz \
	--vcf_output /lscratch/$SLURM_JOB_ID/$sample/merge_output.pead.vcf.gz \
	--threads 1
bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u \
	/lscratch/$SLURM_JOB_ID/$sample/merge_output.pead.vcf.gz \
	| bcftools annotate --set-id '%CHROM-%POS-%REF-%ALT' --output-type v --no-version \
	| bcftools norm -d exact --output-type z -o clair3/$sample.vcf.gz
tabix -f -p vcf clair3/$sample.vcf.gz
convert2annovar.pl -format vcf4old clair3/$sample.vcf.gz -includeinfo --outfile clairAnnotation/$sample.avinput
ver=hg38
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
sed -i "1 s/Otherinfo1\tOtherinfo2\tOtherinfo3\tOtherinfo4\tOtherinfo5\tOtherinfo6\tOtherinfo7\tOtherinfo8\tOtherinfo9\tOtherinfo10/CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGT_FIELDS/" clairAnnotation/$sample.avinput."$ver"_multianno.txt
awk -F "\t" -v NAME="$sample" 'BEGIN{OFS="\t"} NR==1 {print "Sample", $11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$6,$7,$8,$9,$10, "Note"; next} NR>1 {print NAME, $11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$6,$7,$8,$9,$10, ""}' clairAnnotation/$sample.avinput."$ver"_multianno.txt > clairAnnotation/$sample.annovar.tsv
rm clairAnnotation/$sample.avinput clairAnnotation/$sample.avinput."$ver"_multianno.txt
done