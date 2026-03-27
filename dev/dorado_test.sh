

#dorado works on fast5 files as well.
#module load pod5
#find fast5_pass/ -name "*.fast5" -exec pod5 convert fast5 {} --output pod5/ --one-to-one fast5_pass \;
#Or put all fast5 in the same folder, then
#pod5 convert fast5 fast5_pass/*.fast5 --output pod5/ --one-to-one fast5_pass/
module load dorado
while read -r filename; do
	sample=$(basename $filename | cut -d "." -f 1 | sed "s/-/_/")
	dorado basecaller  ${DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 pod5/$filename --trim all --emit-fastq | gzip > supAccu/fastq/$sample.fastq.gz
	done < pod5files.txt
	
while read -r filename; do
	sample=$(basename $filename | cut -d "_" -f 3 | sed "s/-/_/")
	dorado basecaller  ${DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 pod5/$filename --trim all --emit-fastq | gzip > fastq/$sample.fastq.gz
	done < pod5.txt
#directory renamed as bam_mm2

#multiple pod5 files for each sample
for barcodeDir in pod5_pass/barcode*/; do 
	sample=$(echo $barcodeDir | cut -d "/" -f 2)
	dorado basecaller  ${DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 $barcodeDir | gzip > fastq/$sample.fastq.gz
	done


while read -r filename; do
	sample=$(basename $filename | cut -d "_" -f 3 | sed "s/-/_/")
	dorado basecaller ${DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_hac@v5.0.0 pod5/$filename --trim all --emit-fastq | gzip > fastq/$sample.fastq.gz
	done < pod5.txt
#Unfortunately it led to very low filtered fast reads when using the hac.
	
#does not use Chopper here:
while read -r filename; do
	sample=$(basename $filename | cut -d "_" -f 3 | sed "s/-/_/");
	dorado basecaller --trim all -Y --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.mmi ${DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 pod5/$filename \
	| sambamba sort -u --compression-level 6 --tmpdir=$WORK_DIR -t 3 -o bam/$sample.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin);
	done < pod5.txt

