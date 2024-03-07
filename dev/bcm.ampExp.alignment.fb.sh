#!/bin/bash
#SBATCH -c8
#SBATCH --mem=64g
#SBATCH --gres=lscratch:50
#SBATCH --time=2:0:0

#~20 min per sample, thus adjust the time accordingly. 
#fastq files have to be SampleID-LW and SampleID-MW. If missing either LW or MW, touch an empty file.

set -e

module load minimap2/2.26 sambamba/0.8.2

mkdir -p bam

for fastq in fastq/*-LW*.fastq.gz; do
sample=$(basename $fastq | cut -d "_" -f 1 | sed "s/-LW//");
echo -e "$sample" >> manifest.csv
minimap2 -a -x map-ont -Y -t 3 -R "@RG\tLB:QBplEx\tID:$sample\tSM:$sample" /data/OGL/resources/genomes/GRCh38/OPN1LW/GRCh38Decoy_OPN1LW.mmi $fastq | sambamba sort -u --compression-level 6 --tmpdir=/lscratch/$SLURM_JOB_ID -t 3 -o bam/$sample.LW.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin);
done

for fastq in fastq/*-MW*.fastq.gz; do
sample=$(basename $fastq | cut -d "_" -f 1 | sed "s/-MW//");
minimap2 -a -x map-ont -Y -t 3 -R "@RG\tLB:QBplEx\tID:$sample\tSM:$sample" /data/OGL/resources/genomes/GRCh38/OPN1MW/GRCh38Decoy_OPN1MW.mmi $fastq | sambamba sort -u --compression-level 6 --tmpdir=/lscratch/$SLURM_JOB_ID -t 3 -o bam/$sample.MW.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin);
done

while read -r sample;
do
sambamba merge -t 3 -l 6 bam/$sample.bcm.bam bam/$sample.LW.bam bam/$sample.MW.bam;
mv bam/$sample.bcm.bam.bai bam/$sample.bcm.bai;
done < manifest.csv

rm bam/*.LW.bam bam/*.LW.bam.bai bam/*.MW.bam bam/*.MW.bam.bai

#add coverage
#test freebayes

module load samtools/1.17 vcflib/1.0.3 vt/0.57721

module load freebayes/1.3.5 R/4.3.0 annovar/2020-06-08
mkdir -p freebayes
mkdir -p annotation

#changed -C 5 to 3 to account for low coverage Plasmid Express reads.
while read -r sample;
do
freebayes -f /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --max-complex-gap 80 -p 10 -C 3 -F 0.05 \
--genotype-qualities --strict-vcf --use-mapping-quality --min-base-quality 5 \
--targets /data/OGL/resources/bed/OPN1_e2-e5.bed \
bam/$sample.bcm.bam \
| vcffilter -f "QUAL > 10" \
| bcftools +fill-tags - -Ou -- -t VAF \
| bcftools norm --multiallelics -any --output-type u - \
| bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa  --output-type u --no-version - \
| bcftools annotate --set-id 'haplo_%CHROM\:%POS%REF\>%ALT' -x ^INFO/AO,^FORMAT/GT,FORMAT/DP,FORMAT/VAF --output-type u --no-version \
| bcftools norm -d exact --output-type z -o freebayes/$sample.haplo.vcf.gz
tabix -p vcf freebayes/$sample.haplo.vcf.gz
freebayes -f /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --limit-coverage 1000 --min-alternate-fraction 0.05 -C 3 -p 10 \
--min-mapping-quality 0 --genotype-qualities --strict-vcf --use-mapping-quality --min-base-quality 5 \
--targets /data/OGL/resources/bed/OPN1_geneBody_LR-PCR.bed \
bam/$sample.bcm.bam \
| vcffilter -f "QUAL > 10 & SAF > 0 & SAR > 0 & RPR > 0 & RPL > 0 & AO > 2 & DP > 5" \
| bcftools +fill-tags - -Ou -- -t VAF \
| bcftools norm --multiallelics -any --output-type u - \
| vt decompose_blocksub -p -m -d 2 - \
| bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u - \
| bcftools annotate --set-id 'sm_%CHROM\:%POS%REF\>%ALT' -x ^INFO/AO,^FORMAT/GT,FORMAT/DP,FORMAT/VAF --output-type v --no-version \
| bcftools norm -d exact --output-type z -o freebayes/$sample.small.vcf.gz
tabix -p vcf freebayes/$sample.small.vcf.gz
bcftools concat -a --rm-dups none --no-version \
freebayes/$sample.haplo.vcf.gz freebayes/$sample.small.vcf.gz \
-Oz -o freebayes/$sample.vcf.gz
tabix -p vcf freebayes/$sample.vcf.gz
convert2annovar.pl -format vcf4old freebayes/$sample.vcf.gz -includeinfo --outfile annotation/$sample.avinput
ver=hg38
table_annovar.pl annotation/$sample.avinput \
$ANNOVAR_DATA/$ver \
-buildver $ver \
-remove \
-out annotation/$sample.avinput \
--protocol refGeneWithVer \
-operation g \
--argument '-hgvs' \
--polish -nastring . \
--thread 1 \
--otherinfo
sed -i "1 s/Otherinfo1\tOtherinfo2\tOtherinfo3\tOtherinfo4\tOtherinfo5\tOtherinfo6\tOtherinfo7\tOtherinfo8\tOtherinfo9\tOtherinfo10/CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGT_FIELDS/" annotation/$sample.avinput."$ver"_multianno.txt
Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/bcmlocus.R \
/data/OGL/resources/bcmlocus.xlsx \
$sample annotation/$sample.avinput."$ver"_multianno.txt annotation/$sample.bcm.tsv
rm annotation/$sample.avinput annotation/$sample.avinput."$ver"_multianno.txt
done < manifest.csv

#perhaps it's better to use the freebayes raw output to represent the variants, as bcftools norm may shorten the variants.

#rm annotation/$sample.avinput."$ver"_multianno.txt


#find LW/fastq -name "*fastq*" | parallel -j 10 -I% --max-args 1 --tag "minimap2 -a -Y -t 3 /data/OGL/resources/genomes/GRCh38/OPN1LW/GRCh38Decoy_OPN1LW.mmi % | sambamba sort -u --compression-level 6 --tmpdir=/lscratch/$SLURM_JOB_ID -t 3 -o %.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin)"

#find MW/fastq -name "*fastq*" | parallel -j 10 -I% --max-args 1 --tag "minimap2 -a -Y -t 3 /data/OGL/resources/genomes/GRCh38/OPN1MW/GRCh38Decoy_OPN1MW.mmi % | sambamba sort -u --compression-level 6 --tmpdir=/lscratch/$SLURM_JOB_ID -t 3 -o MW/bam/%.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin)"

#The minimap2 index were generated by: minimap2 -d genome.mmi  /fdb/igenomes/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa 
#The index is now moved to /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.mmi


