#!/bin/bash
#SBATCH -c56
#SBATCH --mem=128g
#SBATCH --gres=lscratch:100
#SBATCH --time=4:0:0

#Took 12 min, 52G mem, 8G lscratch for a GS trio.
#Took 50 min, 355G mem for 50 GS samples. Can use norm for this now.
#DeepVariant gvcf files need to be in the deepvariant/gvcf folder
set -e

config_file=$1 # exome/es/wes or genome/gs/wgs
ngstype=$(grep "^ngstype:" $config_file | head -n 1 | cut -d"'" -f 2)
ref_genome=$(grep "^ref_genome:" $config_file | head -n 1 | cut -d"'" -f 2)
analysis_batch_name=$(grep "^analysis_batch_name:" $config_file | head -n 1 | cut -d"'" -f 2)
gt_call_version=$(grep "^gt_call_version:" $config_file | head -n 1 | cut -d"'" -f 2)
module load $(grep "^glnexus_version:" $config_file | head -n 1 | cut -d"'" -f 2)
module load $(grep "^samtools_version:" $config_file | head -n 1 | cut -d"'" -f 2) 

WORK_DIR=/lscratch/$SLURM_JOB_ID

case "${ngstype^^}" in
	"EXOME"|"WES"|"ES")
	glnexus_cli --dir /lscratch/$SLURM_JOB_ID/glnexus \
		--config DeepVariantWES --bed $(grep "^padded_bed:" $config_file | head -n 1 | cut -d"'" -f 2) \
		--threads 24  --mem-gbytes 96 \
		deepvariant/gvcf/*.g.vcf.gz \
		| bcftools norm --multiallelics -any --output-type u --no-version \
		| bcftools norm -d exact --check-ref s --fasta-ref $ref_genome --output-type u --no-version - \
		| bcftools +fill-tags - -Ou -- -t AC,AC_Hom,AC_Het,AN,AF \
		| bcftools annotate --threads 24 --set-id 'dvg_%CHROM\:%POS%REF\>%ALT' --no-version - -Oz -o deepvariant/$analysis_batch_name.$gt_call_version.vcf.gz
	tabix -f -p vcf deepvariant/$analysis_batch_name.$gt_call_version.vcf.gz
	;;
	"GENOME"|"WGS"|"GS")
	glnexus_cli --dir /lscratch/$SLURM_JOB_ID/glnexus --config DeepVariant \
		--threads 54 --mem-gbytes 128 \
		deepvariant/gvcf/*.g.vcf.gz \
		| bcftools view --threads 54 -Ob -o $WORK_DIR/$analysis_batch_name.glnexus.bcf
	df -h $WORK_DIR
	bcftools norm --multiallelics -any --output-type u --no-version $WORK_DIR/$analysis_batch_name.glnexus.bcf \
		| bcftools norm -d exact --check-ref s --fasta-ref $ref_genome --output-type u --no-version - \
		| bcftools +fill-tags - -Ou -- -t AC,AC_Hom,AC_Het,AN,AF \
		| bcftools annotate --threads 54 --set-id 'dvg_%CHROM\:%POS%REF\>%ALT' --no-version - -Oz -o deepvariant/$analysis_batch_name.$gt_call_version.vcf.gz
	tabix -f -p vcf deepvariant/$analysis_batch_name.$gt_call_version.vcf.gz
	;;
	"GENOME-CLAIR3")
	glnexus_cli --dir /lscratch/$SLURM_JOB_ID/glnexus --config /data/OGL/resources/clair3/clair3.yml \
		--threads 54 --mem-gbytes 128 \
		clair3/gvcf/*.gvcf.gz \
		| bcftools view --threads 54 -Ob -o $WORK_DIR/$analysis_batch_name.glnexus.bcf
	df -h $WORK_DIR
	bcftools norm --multiallelics -any --output-type u --no-version $WORK_DIR/$analysis_batch_name.glnexus.bcf \
		| bcftools norm -d exact --check-ref s --fasta-ref $ref_genome --output-type u --no-version - \
		| bcftools +fill-tags - -Ou -- -t AC,AC_Hom,AC_Het,AN,AF \
		| bcftools annotate --threads 54 --set-id 'clr3g_%CHROM\:%POS%REF\>%ALT' --no-version - -Oz -o clair3/$analysis_batch_name.$gt_call_version.vcf.gz
	tabix -f -p vcf clair3/$analysis_batch_name.$gt_call_version.vcf.gz
	;;
esac
