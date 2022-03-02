#!/bin/bash

# only the first @RG from the old bam/cram is used.
# on Notepad++ the RG \t shows as \t. If editing manually, keep one \t.  This is different to fastq file metadata_file where it showed as \\t. Because  different processings in snakefile, do not switch these types.

module load samtools/1.11
 
for file in old_bam/*.*am; do filename=$(basename $file); sample=$(echo $filename | cut -d. -f 1); RG=$(samtools view -H $file | grep "^@RG" | head -n 1 | sed 's/\t/\\t/g'); echo $sample,$filename,"$RG" >> metadata_file.csv; done

#awk -F"," 'BEGIN{OFS=","} {sub("HUF_", "", $1); sub("_D1", "", $1); sub("HUF_", "", $3); sub("_D1", "", $3); print $0 }' metadata_file.csv > metadata_file.csv.bk && mv metadata_file.csv.bk metadata_file.csv