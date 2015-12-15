#!/bin/bash

module load VEP/82


input_vcf=$1
genome=$2


if [ "$genome" == "GRCh38" ] || [ "$genome" == "GRCh37" ]; then
	variant_effect_predictor.pl -i $input_vcf --offline \
	--cache --dir_cache $VEPCACHEDIR \
	--fasta $VEPCACHEDIR/human.fa --species human --assembly $genome  \
	--output ${input_vcf%.vcf}.VEP82everything.$genome.vcf \
	--everything --vcf --force_overwrite --fork 10
else
    echo "Pick either GRCh38 or GRCh37 genomes"
fi
