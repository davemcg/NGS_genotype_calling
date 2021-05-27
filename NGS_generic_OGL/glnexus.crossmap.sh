#!/bin/bash
#SBATCH --gres=lscratch:200
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH --time=6:0:0


module load samtools/1.12 crossmap/0.5.4

WORK_DIR=/lscratch/$SLURM_JOB_ID

hg38ref=/data/OGL/resources/genomes/GRCh38/GRCh38Decoy2.fa
hg19refM=/data/OGL/resources/1000G_phase2_GRCh37/human_g1k_v37_decoyM.fasta
hg19ref=/data/OGL/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta

bcftools norm --multiallelics -any --output-type u --no-version $1 \
| bcftools norm --check-ref s --fasta-ref $hg38ref --output-type u --no-version - \
| bcftools +fill-tags - -Ou -- -t AC,AC_Hemi,AC_Hom,AC_Het,AN,AF \
| bcftools annotate --threads 16 --set-id 'dv_%CHROM\:%POS%REF\>%ALT' -x FORMAT/PL,FORMAT/RNC \
--no-version - -Oz -o $WORK_DIR/dv.vcf.gz

#AF from glnexus seems to be strange. Asked at GitHub: https://github.com/dnanexus-rnd/GLnexus/issues/247
#bcftools AF above is AC/AN. AC_Hom is the no. of Allele count, while gnomAD AC_Hom could be the no. of individuals (at least gnomAD website has the no. of individuals)

crossmap vcf /data/OGL/resources/ucsc/hg38ToHg19.over.chain.gz \
$WORK_DIR/dv.vcf.gz \
$hg19refM \
$WORK_DIR/GRCh37.vcf && rm $WORK_DIR/dv.vcf.gz

time sed -e 's/^chrM/MT/' -e 's/<ID=chrM/<ID=MT/' $WORK_DIR/GRCh37.vcf \
| sed -e 's/^chr//' -e 's/<ID=chr/<ID=/' - \
| bcftools norm --check-ref s --fasta-ref $hg19ref --output-type u - \
| bcftools sort -m 48G -T $WORK_DIR -Ou - \
| bcftools norm --threads 16 -d exact --output-type z - -o ${1%.vcf.gz}.GRCh37.vcf.gz 

#sort before -d exact, otherwise when two sites from hg38 to a single hg19 site will not be removed. -d exact only removes neighboring line.

tabix -f -p vcf ${1%.vcf.gz}.GRCh37.vcf.gz
