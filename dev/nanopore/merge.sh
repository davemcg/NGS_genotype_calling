#!/bin/bash
module load sambamba
sambamba merge -t 8 D15.bam D15.xab.bam D15.xac.bam D15.xad.bam D15.xag.bam D15.xae.bam D15.xah.bam D15.xaf.bam D15.xai.bam

rm D15.xab.bam D15.xac.bam D15.xad.bam D15.xag.bam D15.xae.bam D15.xah.bam D15.xaf.bam D15.xai.bam
rm D15.xab.bam.bai D15.xac.bam.bai D15.xad.bam.bai D15.xag.bam.bai D15.xae.bam.bai D15.xah.bam.bai D15.xaf.bam.bai D15.xai.bam.bai

module load R/3.6.3 mosdepth

mosdepth -t 10 --no-per-base --by /data/OGL/resources/bed/OGLv1.GRCh38.sorted.bed --use-median --mapq 0 --fast-mode --thresholds 10,20,30 D15.md ../D15.bam


#variant calls.

sinteractive --mem=32g -c10 --gres=lscratch:200

module load PEPPER_deepvariant

run_pepper_margin_deepvariant call_variant -b D15.bam -f /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa -o D15dv -t 10 --ont_r9_guppy5_sup --sample_name D15 --gvcf --phased_output 

