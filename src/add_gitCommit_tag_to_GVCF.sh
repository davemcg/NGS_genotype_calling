#!/bin/bash

gvcf_file=$1
git_repo_url=$2
git_commit=$3

tabix -H $gvcf_file > /scratch/mcgaugheyd/$1.header

sed -i "$ i\\##git_repo_url:$git_repo_url" /scratch/mcgaugheyd/$1.header
sed -i "$ i\\##git_commit:$git_commit" /scratch/mcgaugheyd/$1.header
tabix -r /scratch/mcgaugheyd/$1.header $1 > /scratch/mcgaugheyd/$1.newvcf
mv /scratch/mcgaugheyd/$1.newvcf $1
rm $1.tbi
tabix -p vcf $1


