#!/bin/bash

module load kallisto/0.42.4

# indexes
# /data/mcgaugheyd/genomes/GRCm38/kallisto_index/gencode.vM8.idx 
# /data/mcgaugheyd/genomes/GRCh38/0.42.4_Homo_sapiens.GRCh38.rel84.cdna.all.idx

index=$1
index=$(echo $index | tr 'a-z' 'A-Z')
read_1=$2
read_2=$3

if [ $index == 'HUMAN' ]; then
		kallisto quant -i /data/mcgaugheyd/genomes/GRCh38/0.42.4_Homo_sapiens.GRCh38.rel84.cdna.all.idx \
		   -o ${2%_*}_kallisto \
           -b 100 \
           -t $SLURM_CPUS_PER_TASK \
		   $read_1 $read_2 
elif [ $index == 'MOUSE' ]; then
		kallisto quant -i /data/mcgaugheyd/genomes/GRCm38/kallisto_index/gencode.vM8.idx \
    	   -o ${2%_*}_kallisto \
           -b 100 \
           -t $SLURM_CPUS_PER_TASK \
           $read_1 $read_2 
else
	echo "Give either Mouse or Human"
fi

