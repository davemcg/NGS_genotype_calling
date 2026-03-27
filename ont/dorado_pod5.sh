#!/bin/bash
#SBATCH -c8
#SBATCH --mem=16g
#SBATCH --partition=gpu
#SBATCH --gres=lscratch:50,gpu:1 --constraint='gpuv100|gpua100|gpuv100x'
#SBATCH --time=6:0:0

#dorado works on fast5 files as well.
#copy data exactly as pod5_pass/barcode01 etc.
#save rename.sh using info from the excel template to rename the files. The rename.sh should be under the project folder.
#took 4-10 minutes for batch3 rerun (5 pod5 files), similarly to batch 4 (4h25min) which run longer and had 38 pod5 files. Thus, increase time to 5 hours. 

set -e
module load dorado/0.8.1
mkdir -p fastq

rm -rf pod5_pass/mixed pod5_pass/unclassified

#multiple pod5 files for each sample
for barcodeDir in pod5_pass/barcode*/; do 
	sample=$(echo $barcodeDir | cut -d "/" -f 2)
	dorado basecaller  ${DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 $barcodeDir --trim all --emit-fastq | gzip > fastq/$sample.fastq.gz
	done

cp rename.sh fastq/
cd fastq
bash rename.sh
rm rename.sh
cd ..
sstat -j $SLURM_JOB_ID
sacct -j $SLURM_JOB_ID
