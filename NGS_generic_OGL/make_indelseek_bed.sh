#!/bin/bash

bed=$1

awk -F"\t" 'NR == 1 {print "#region"} NR > 1 {print $1 ":" $2 "-" $3}' $bed > ${bed%.bed}.indelseek.bed 

#$1, $2 will have space between the fields or FS if specifying FS