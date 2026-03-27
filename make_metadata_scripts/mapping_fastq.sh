#!/bin/bash

calling_configure_file=$1
mapping_file=$2

set -e

sed -i 's/\r$//' $mapping_file

analysis_batch_name=$(grep "^analysis_batch_name:" $calling_configure_file | head -n 1 | cut -d"'" -f 2)
lib=$(grep "^lib:" $calling_configure_file | head -n 1 | cut -d"'" -f 2)
metadata_file=$(grep "^metadata_file:" $calling_configure_file | head -n 1 | cut -d"'" -f 2)

#awk -F'\t' '{print $1 "\t" $2}' $mapping_file | 
while IFS=$'\t' read -r sm filename; do
  if [[ ! -f fastq/$filename ]]; then
    echo "First line is header"
    continue
  fi
  header=$(zcat fastq/$filename | head -n 1)
  id=$(echo $header | cut -d: -f 3,4 | sed 's/\:/\./g')
  echo "$sm,$filename,@RG\\\tID:$id"_"$sm\\\tSM:$sm\\\tLB:$lib"_"$sm\\\tPL:ILLUMINA" >> $metadata_file
done < $mapping_file
