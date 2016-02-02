#!/bin/bash

module load gemini/0.18

VCF=$1
PED=$2
DBNAME=$3

gemini load --cores 8 -t VEP -v $VCF -p $PED $DBNAME
