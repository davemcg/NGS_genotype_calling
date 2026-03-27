#!/bin/bash
#SBATCH --gres=lscratch:10
#SBATCH --cpus-per-task=8
#SBATCH --mem=8g
#SBATCH --time=2:00:00


rm -rf .snakemake
rm -rf 00log gemini_tsv/
rm -rf slurm*.out
rm -rf *.gt3.anno3.dvg.vcf.gz*
rm *.gemini.db

#rm -rf old_bam bam cram
echo "File deletion task done"

