#!/bin/bash
#SBATCH --gres=lscratch:10
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --time=2:00:00

# only the first @RG from the old bam/cram is used.
# on Notepad++ the RG \t shows as \t. If editing manually, keep one \t.  This is different to fastq file metadata_file where it showed as \\t. Because  different processings in snakefile, do not switch these types.

module load samtools/1.21 #orignal 1.11

for file in old_bam/*.*am; do filename=$(basename $file); sample=$(echo $filename | cut -d. -f 1); RG=$(samtools view -H $file | grep "^@RG" | head -n 1 | sed 's/\t/\\t/g'); ref=$(samtools view -H $file  | grep "^@PG" | head -n 1 | cut -d " " -f 6); echo $sample,$filename,"$RG","$ref" >> metadata_file.csv; done

#awk -F"," 'BEGIN{OFS=","} {sub("HUF_", "", $1); sub("_D1", "", $1); sub("HUF_", "", $3); sub("_D1", "", $3); print $0 }' metadata_file.csv > metadata_file.csv.bk && mv metadata_file.csv.bk metadata_file.csv

#remove accession number (ie. 2264_) from sample name
#awk -F"," 'BEGIN{OFS=","} {sub("^[0-9]+_", "", $1); sub("SM:[0-9]+_", "SM:", $3); print $0 }' metadata_file.csv > metadata_file.csv.bk && mv metadata_file.csv.bk metadata_file.csv
