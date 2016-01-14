#!/bin/bash

module load VEP/82


input_vcf=$1
genome=$2
cores=$3

if [ "$genome" == "GRCh38" ] || [ "$genome" == "GRCh37" ]; then
	variant_effect_predictor.pl -i $input_vcf --offline \
	--cache --dir_cache $VEPCACHEDIR \
	--fasta $VEPCACHEDIR/$genome.fa --species human --assembly $genome  \
	--output ${input_vcf%.vcf}.VEP82everything.$genome.vcf \
	--plugin ExAC,/data/mcgaugheyd/genomes/hg19/ExAC.r0.3.nonTCGA.sites.vep.vcf.gz,AC \
	--plugin Grantham \
	--total_length \
	--fields Consequence, Codons, Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,CANONICAL,CADD_raw,CADD_phred,Grantham
	--everything --vcf --force_overwrite --fork $cores
else
    echo "Pick either GRCh38 or GRCh37 genomes"
fi
