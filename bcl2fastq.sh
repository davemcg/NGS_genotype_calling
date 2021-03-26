#!/bin/bash
#SBATCH --cpus-per-task=20
#SBATCH --mem=128g
#SBATCH -p norm,quick
#SBATCH --time=2:00:00

#sinteractive --cpus-per-task=20 --mem=128g
#finished one midi-run in 5-10 min.
#samplesheet should be named as SampleSheet.csv if no sample sheet is specified in the command below
#$1 Default = 1 mismatches for dual index; Use 0 mismatch for single 6nt indexes ( 1 = 1,1; 0 = 0,0)

module load bcl2fastq/2.20.0

bcl2fastq --runfolder-dir . --output-dir ./bcl2fastq -r 4 -w 4 -p 12 --barcode-mismatches $1
#read 4, write 4, processing 12

mkdir -p bcl2fastq/undertmined
mv bcl2fastq/Undetermined* bcl2fastq/undertmined/.