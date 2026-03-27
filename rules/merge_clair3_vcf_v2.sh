#!/bin/bash
#SBATCH -c8
#SBATCH --mem=16g
#SBATCH --gres=lscratch:100
#SBATCH --time=4:0:0

config_file=$1
glnexus_vcf=$2
ref_genome=$(grep "^ref_genome:" $config_file | head -n 1 | cut -d"'" -f 2)
analysis_batch_name=$(grep "^analysis_batch_name:" $config_file | head -n 1 | cut -d"'" -f 2)
gt_call_version=$(grep "^gt_call_version:" $config_file | head -n 1 | cut -d"'" -f 2)

set -e
module load samtools
LWORK_DIR=/lscratch/$SLURM_JOB_ID

bcftools merge --merge none --missing-to-ref --output-type u --threads 8 clair3/vcf/*.clr3.filtered.vcf.gz \
	| bcftools annotate --threads 8 --set-id 'clr3_%CHROM\:%POS%REF\>%ALT' --no-version - -Ou \
	| bcftools +fill-tags - -Ov -- -t AC,AC_Hom,AC_Het,AN,AF \
	| sed 's#0/0:\.:\.:\.#0/0:10:10:10,0#g' - \
	| bgzip -f > clair3/$analysis_batch_name.clr3.vcf.gz
tabix -f -p vcf clair3/$analysis_batch_name.clr3.vcf.gz
mkdir -p prioritization-clair3g
bcftools isec --threads 8 -p $LWORK_DIR -n=2 -w1 --collapse none --output-type z \
  clair3/$analysis_batch_name.clr3.vcf.gz $glnexus_vcf

bcftools annotate --threads 8 --write-index=tbi --set-id 'clr3g_%CHROM\:%POS%REF\>%ALT' \
  -a $glnexus_vcf -c INFO \
  $LWORK_DIR/0000.vcf.gz -Oz -o prioritization-clair3g/$analysis_batch_name.$gt_call_version.clr3.vcf.gz
