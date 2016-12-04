#!/bin/bash

module load VEP/86


input_vcf=$1
genome=$2
cores=$3

if [ "$genome" == "GRCh38" ] || [ "$genome" == "GRCh37" ]; then
	variant_effect_predictor.pl -i $input_vcf --offline \
	--cache --dir_cache $VEPCACHEDIR \
	--fasta $VEPCACHEDIR/$genome.fa --species human --assembly $genome  \
	--output ${input_vcf%.vcf}.VEP.$genome.vcf \
	--plugin Grantham \
	--plugin MaxEntScan,/home/mcgaugheyd/bin/MaxEntScan \
	--plugin Phenotypes,file=$VEPCACHEDIR/Phenotypes.pm_human_86_GRCh37.bed.gz,include_sources=HGMD_PUBLIC \
	--total_length \
    --hgvs \
	--sift b \
    --polyphen b \
    --symbol \
    --numbers \
    --biotype \
    --total_length \
	--pubmed \
    --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,CANONICAL,Grantham,MaxEntScan,HGVSc,HGVSp,PUBMED,Phenotypes \
	--vcf --force_overwrite --fork $cores
else
    echo "Pick either GRCh38 or GRCh37 genomes"
fi
