#!/bin/bash

echo -e "sample\tlength\tmean\ton_target_rate" > coverage/$1.mean.coverage.summary.tsv

for file in coverage/mosdepth/*.md.mosdepth.summary.txt; do
	filename=$(basename $file)
	sm=$(echo $filename | sed 's/.md.mosdepth.summary.txt//')
	on_target=$(tail -n 2 $file | awk 'NR>1{print $3/p} {p=$3}')
	tail -n 1 $file | awk -v on_target="$on_target" -v sm="$sm" -F"\t" 'BEGIN{OFS="\t"}{print sm,$2,$4,on_target}' - >> coverage/$1.mean.coverage.summary.tsv
done