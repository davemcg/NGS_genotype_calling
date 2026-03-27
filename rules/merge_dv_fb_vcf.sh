#!/bin/bash
#SBATCH --gres=lscratch:200
#SBATCH --cpus-per-task=16
#SBATCH --mem=16g
#SBATCH --partition=norm
#SBATCH --time=6:0:0

set -e
config=$1
new_batch_name=$2
WORK_DIR=/lscratch/$SLURM_JOB_ID
#to merge the three vcf types from 2 or more Step 1 runs (NGS_genotype_calling).

module load $(grep "^samtools_version:" $config | head -n 1 | cut -d"'" -f 2)
ngstype=$(grep "^ngstype:" $config | head -n 1 | cut -d"'" -f 2)

mkdir -p $WORK_DIR/input
vcf_inputs=""
for vcf in freebayes/*.freebayes.vcf.gz; do
	vcf_inputs+="$vcf "; done
bcftools merge --merge none --missing-to-ref --output-type z --threads 8 $vcf_inputs -o $WORK_DIR/input/freebayes.vcf.gz && tabix -p vcf /$WORK_DIR/input/freebayes.vcf.gz &

vcf_inputs=""
case "${ngstype^^}" in
 "EXOME"|"WES"|"ES"|"PANEL")
  for vcf in deepvariant/*.dv.phased.vcf.gz; do
	vcf_inputs+="$vcf "; done
  ;;
 *)
  for vcf in deepvariant/*.dv.hf.vcf.gz; do
	vcf_inputs+="$vcf "; done
  ;;
esac
bcftools merge --merge none --missing-to-ref --output-type z --threads 8 $vcf_inputs -o $WORK_DIR/input/dv.hf.vcf.gz
tabix -p vcf $WORK_DIR/input/dv.hf.vcf.gz

vcf_inputs=""
for vcf in deepvariant/*.dv.glnexus.phased.vcf.gz; do
	vcf_inputs+="$vcf "; done
bcftools merge --merge none --missing-to-ref --output-type v --threads 8 $vcf_inputs \
	| sed 's#0/0:\.:\.:\.#0/0:10:10:10,0#g' - \
	| bcftools +fill-tags -Oz --threads 8 -o $WORK_DIR/input/dv.glnexus.phased.vcf.gz -- -t AC,AC_Hom,AC_Het,AN,AF
tabix -p vcf $WORK_DIR/input/dv.glnexus.phased.vcf.gz


bcftools isec -p $WORK_DIR/dv -w 2 --collapse none --output-type u --threads 16 \
	$WORK_DIR/input/dv.glnexus.phased.vcf.gz \
	$WORK_DIR/input/dv.hf.vcf.gz
bcftools +fill-tags $WORK_DIR/dv/0001.bcf -Ov -- -t AC,AC_Hom,AC_Het,AN,AF \
	| sed 's#0/0:\.:\.:\.#0/0:10:10:10,0#g' - \
	| bcftools annotate --threads 16 --set-id 'dv_%CHROM\:%POS%REF\>%ALT' --no-version - -Oz -o $WORK_DIR/dv/dv.hf.vcf.gz
tabix -f -p vcf $WORK_DIR/dv/dv.hf.vcf.gz
bcftools concat --threads 16 -a --rm-dups none --no-version \
	$WORK_DIR/input/dv.glnexus.phased.vcf.gz $WORK_DIR/dv/dv.hf.vcf.gz -Oz \
	-o $WORK_DIR/dv.glnexus.hf.vcf.gz
tabix -f -p vcf $WORK_DIR/dv.glnexus.hf.vcf.gz
rm -r -f $WORK_DIR/dv $WORK_DIR/input/dv.glnexus.phased.vcf.gz*
bcftools isec --threads 16 -p $WORK_DIR --collapse none -Oz \
	$WORK_DIR/dv.glnexus.hf.vcf.gz \
	$WORK_DIR/input/freebayes.vcf.gz
rm $WORK_DIR/0003.vcf.gz &
bcftools annotate --threads 16 --set-id 'fb_%CHROM\:%POS%REF\>%ALT' -x ^INFO/QA,FORMAT/RO,FORMAT/QR,FORMAT/AO,FORMAT/QA,FORMAT/GL \
	--no-version $WORK_DIR/0001.vcf.gz -Ov - \
	| sed 's#0/0:\.:\.:\.#0/0:10:10:10,0#g' - \
	| bcftools +fill-tags - -Oz -o $WORK_DIR/fb.vcf.gz -- -t AC,AC_Hom,AC_Het,AN,AF
rm $WORK_DIR/0001.vcf.gz &
zcat $WORK_DIR/0002.vcf.gz | sed 's/\tdv/\tfbDv/' | bgzip -@ 16 > $WORK_DIR/dvFb.vcf.gz
rm $WORK_DIR/0002.vcf.gz &
tabix -f -p vcf $WORK_DIR/0000.vcf.gz
tabix -f -p vcf $WORK_DIR/fb.vcf.gz
tabix -f -p vcf $WORK_DIR/dvFb.vcf.gz
bcftools concat --threads 16 -a --rm-dups none --no-version \
	$WORK_DIR/dvFb.vcf.gz $WORK_DIR/0000.vcf.gz $WORK_DIR/fb.vcf.gz -Oz \
	-o prioritization/$new_batch_name.vcf.gz
tabix -f -p vcf prioritization/$new_batch_name.vcf.gz
