#!/bin/bash
cd fastq
ls > ../rename.sh
cd ..
cut -d "_" -f 2- rename.sh | paste -d " " rename.sh - | sed 's/^/mv /' > fastq/rename.sh
rm rename.sh
cd fastq
bash rename.sh
