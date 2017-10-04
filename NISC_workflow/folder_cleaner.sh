#!/bin/bash

# dangerous script designed to keep core files from a NISC processed exome set
# will keep:
# 1. Commands run
# 2. the final realigned, recalibrated, mark-duped bam and index
# 3. The gvcf and index
# EVERYTHING ELSE IS DELETED!!!!!!!!!!

# Expected Usage
# ~/git/NGS_genotype_calling/NISC_workflow/folder_cleaner.sh .
# Or if you a risk taker and ARE CERTAIN that all sub-folders are just for sequencing data
# for i in */; do ~/git/NGS_genotype_calling/NISC_workflow/folder_cleaner.sh $i; done

function my_date { date "+%Y-%m-%d_%H:%M:%S"; }
cur_date=$(my_date)
path=$1

if [[ -z $path ]]; then
	echo "Path not given"
	exit 1
fi

mkdir $path/to_keep
mv $path/*.realignmentCommands.txt $path/to_keep/
mv $path/*bwa-mem.b37.merged.sorted.markDup.realigned.recalibrated.bam* $path/to_keep/
mv $path/*bwa-mem.b37.merged.sorted.markDup.realigned.raw.g.vcf.gz* $path/to_keep/

find $path -maxdepth 1 -type f > $path/to_keep/deleted_files_$cur_date.txt
find $path -maxdepth 1 -type f -delete

mv $path/to_keep/* $path 
rmdir $path/to_keep/
