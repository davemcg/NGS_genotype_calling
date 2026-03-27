#!/bin/bash
#SBATCH --gres=lscratch:20
#SBATCH --cpus-per-task=2
#SBATCH --mem=32g
#SBATCH --partition=quick
#SBATCH --time=2:0:0

config=$1

#ls cram | grep -v ".crai" | sed s/.cram// > sample.list.txt
#while IFS= read -r line; do sbatch ~/git/NGS_genotype_calling/rules/STRipy.sh config_generic.eG_STGD3.yaml $line; done < sample.list.txt

#selecting mapped reads using subset_mapped_reads.sh reduced sample processing time to 44 hours.
#selecting reads using subset_chr1-M_reads.sh ???
#The bam file has to be saved in the directory sample_bam and named as *.markDup.bam, unless changing the code below.
set -e
if (( $(module list 2>&1 | grep " R/" | wc -l) == 0 )); then module load R/4.3.0; fi
analysis_batch_name=$(grep "^analysis_batch_name:" $config | head -n 1 | cut -d"'" -f 2)

LWORK="/lscratch/$SLURM_JOB_ID"

head -n 1 $(ls STRipy/*.tsv | head -n 1) > STRipy/all.STRipy.$analysis_batch_name.tsv
awk -v IGNORECASE=1 'BEGIN{OFS="\t"} ($0 ~ /(^|[^[:alnum:]_])(intermediate|pathogenic)([^[:alnum:]_]|$)/){print $0}' STRipy/*.STRipy.tsv >> STRipy/all.STRipy.$analysis_batch_name.tsv

if [[ $(cat STRipy/all.STRipy.$analysis_batch_name.tsv | wc -l) -gt 1 ]]; then
 Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/STRipy.R STRipy/all.STRipy.$analysis_batch_name.tsv STRipy/all.STRipy.$analysis_batch_name.xlsx
fi

