#Run WISExome for CNV detection in whole exome data.
#WISExome github page https://github.com/VUmcCGP/wisexome

#Allocate an interactive session and run the program. Sample session:

sinteractive  --mem=32g  --gres=gpu:p100:1,lscratch:100 -c4
module load WISExome
git clone https://github.com/VUmcCGP/wisexome.git
cd wisexome
sed -i 's/sys.stderr.write/# sys.stderr.write/g' test.py
cp $WEXOME_BIN/*.sh .

#Make necessary directories
mkdir input
mkdir leno
mkdir convert
mkdir refdata
mkdir refout
mkdir occout

#Get sample data (examples)

ln -s /fdb/WISExome/WISExome_SRR1273288.bam
ln -s /fdb/WISExome/ucscrefseq.bed
ln -s /fdb/WISExome/SeqCap_EZ_Exome_v3_hg19_capture_targets.bed

#Sort and index the BAM file(s)

samtools sort -T /tmp/WISExome_SRR1273288.sorted -o input/WISExome_SRR1273288.sorted.bam input/WISExome_SRR1273288.bam
samtools index input/WISExome_SRR1273288.sorted.bam

#Produce a file suitable for analysis or training from the BAM file

wexome python consam.py  input/WISExome_SRR1273288.sorted.bam SeqCap_EZ_Exome_v3_hg19_capture_targets.bed convert/WISExome_SRR1273288.sorted.hits
python lennormalize.py convert SeqCap_EZ_Exome_v3_hg19_capture_targets.bed leno

#Create a reference
./create_ref.sh

#Combine the results to make a final selection of reference regions for every region on every target chromosome (will take ~ 1 hour)

python takeref.py refdata/ refout/WISExome_SRR1273288

#Determine unreliable target regions
ln -s convert refsamples
./determine_unreliable_target_regions.sh WISExome_SRR1273288

#Run tests, this will create files: testout/*.flameview.pdf, testout/*.overview.pdf, and testout/*.pickle.

./run_tests.sh WISExome_SRR1273288 SeqCap_EZ_Exome_v3_hg19_capture_targets.bed

#exit
