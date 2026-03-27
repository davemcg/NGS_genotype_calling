#!/bin/bash
#SBATCH --gres=lscratch:400
#SBATCH --cpus-per-task=72
#SBATCH --mem=750g
#SBATCH --partition=largemem
#SBATCH --time=6:0:0

#Required 51G mem, took 21 min for the following 3 sample WGS.
#131 wgs - 512g/480 not sufficient. When using 512g/25l, 200g lscratch is not sufficient.
#current setting is faster than 56 cpus and 1.5TB memory.

module load glnexus/1.2.7 samtools/1.21

# glnexus --dir /lscratch/$SLURM_JOB_ID/glnexus --config DeepVariant \
			# --threads 64 --mem-gbytes 512 \
			# deepvariant/gvcf/*.g.vcf.gz \
			# | bcftools norm --multiallelics -any --output-type u --no-version \
			# | bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/GRCh38/bwa-mem2/GRCh38Decoy2.fa --output-type u --no-version - \
			# | bcftools +fill-tags - -Ou -- -t AC,AC_Hom,AC_Het,AN,AF \
			# | bcftools annotate --threads 64 --set-id 'dv_%CHROM\:%POS%REF\>%ALT' --no-version - -Oz -o  glnexus.512g.vcf.gz

# tabix -f -p vcf glnexus.512g.vcf.gz
# rm -r /lscratch/$SLURM_JOB_ID/glnexus

module load crossmap
hg19refM=/data/OGL/resources/1000G_phase2_GRCh37/human_g1k_v37_decoyM.fasta
hg19ref=/data/OGL/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta
WORK_DIR=/lscratch/$SLURM_JOB_ID
crossmap vcf /data/OGL/resources/ucsc/hg38ToHg19.over.chain.gz \
glnexus.vcf.gz \
$hg19refM \
$WORK_DIR/GRCh37.vcf
sed -e 's/^chrM/MT/' -e 's/<ID=chrM/<ID=MT/' $WORK_DIR/GRCh37.vcf \
| sed -e 's/^chr//' -e 's/<ID=chr/<ID=/' - \
| bcftools norm --check-ref s --fasta-ref $hg19ref --output-type u - \
| bcftools sort -m 512G -T $WORK_DIR/ -Ou - \
| bcftools norm --threads 64 -d exact --output-type z - -o glnexus.512g.GRCh37.vcf.gz
tabix -f -p vcf glnexus.512g.GRCh37.vcf.gz


