#!/bin/bash
#SBATCH -c8
#SBATCH --mem=64g
#SBATCH --gres=lscratch:50
#SBATCH --time=2:0:0

#Could submit with dependency then start automatically after the previous step completes. 
#sbatch --dependency=afterok:11254323 ~/git/NGS_genotype_calling/dev/bcm_RPGR.porechop.minimap2.alignment.sh OGL
#~20 min per sample, thus adjust the time accordingly. 
#fastq files have to be SampleID-LW and SampleID-MW. If missing either LW or MW, touch an empty file.

set -e
library=$1 #QBplEx or OGL; will be part of @RG
metadata=$2 #excel file with sample name, target name, Fwd primer and reverse primer.
batchname=$(date '+%Y%m%d_%H%M') 
WORK_DIR=/lscratch/$SLURM_JOB_ID

##Have not started editing
module load R/4.3.0
Rscript ~/git/NGS_genotype_calling/ont/phasing_metadata.R $metadata $batchname.metadata.tsv

module load ucsc bedtools
cut -f 2- $batchname.metadata.tsv | isPcr -maxSize=30000 -out=bed /data/OGL/resources/ucsc/hg38.2bit stdin $batchname.metadata.bed
paste $batchname.metadata.tsv $batchname.metadata.bed > temp.$batchname.metadata.tsv && mv temp.$batchname.metadata.tsv $batchname.metadata.tsv
#cut -f 2- $batchname.metadata.tsv | isPcr -maxSize=30000 -out=bed /data/OGL/resources/ucsc/hg38.2bit stdin stdout | bedtools merge -i - > $batchname.metadata.bed

module load porechop/0.2.4 chopper/0.7.0 bbtools/39.06 minimap2/2.26 sambamba/0.8.2 mosdepth 

mkdir -p bam porechop mosdepth
export PORECHOP_ADAPTERS=porechop.adapters.py
while read -r line;
do
SAMPLE=$(echo $line | cut -d" " -f 1)
TARGET=$(echo $line | cut -d" " -f 2 | sed 's/_//' )
FWD_PRIMER=$(echo $line | cut -d" " -f 3 | tr 'atcg' 'ATCG')
RevCom=$(echo $FWD_PRIMER | rev | tr 'ATCG' 'TAGC')
sed -e "s/TARGETNAME/$TARGET/" -e "s/PRIMER_SEQ/$FWD_PRIMER/" -e "s/PRIMER_REVCOMPL/$RevCom/" ~/git/NGS_genotype_calling/ont/porechop.adapters.py > porechop.adapters.py
porechop --discard_unassigned --discard_middle --extra_end_trim 0 -i fastq/$SAMPLE.fastq.gz -b porechop/fwd/$SAMPLE
REV_PRIMER=$(echo $line | cut -d" " -f 4 | tr 'atcg' 'ATCG')
RevCom=$(echo $FWD_PRIMER | rev | tr 'ATCG' 'TAGC')
sed -e "s/TARGETNAME/$TARGET/" -e "s/PRIMER_SEQ/$REV_PRIMER/" -e "s/PRIMER_REVCOMPL/$RevCom/" ~/git/NGS_genotype_calling/ont/porechop.adapters.py > porechop.adapters.py
porechop --discard_unassigned --discard_middle --extra_end_trim 0 -i porechop/fwd/$SAMPLE/BC$TARGET.fastq.gz -b porechop/trimmed/$SAMPLE
bbtools readlength in=fastq/$SAMPLE.fastq.gz out=qc/$SAMPLE.readlengh.histogram.txt
bbtools readlength in=porechop/trimmed/$SAMPLE/BC$TARGET.fastq.gz out=qc/$SAMPLE.readlengh.histogram.both.primers.txt
maxl=$(echo $line | awk -F" " '{size = $7- $6 + 200; print size}')
minl=$(echo $line | awk -F" " '{size = $7- $6 - 200; print size}')
zcat porechop/trimmed/$SAMPLE/BC$TARGET.fastq.gz | chopper -q 10 --maxlength $maxl --minlength $minl | gzip > porechop/trimmed/$SAMPLE.length.filtered.fastq.gz
bbtools readlength in=porechop/trimmed/$SAMPLE.length.filtered.fastq.gz out=qc/$SAMPLE.readlengh.histogram.length.filtered.txt
minimap2 -a -x map-ont -Y -t 3 -R "@RG\tLB:$library\tID:$SAMPLE\tSM:$SAMPLE\tPL:ONT" /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.mmi porechop/trimmed/$SAMPLE/BC$TARGET.fastq.gz \
| sambamba sort -u --compression-level 6 --tmpdir=$WORK_DIR -t 3 -o bam/$SAMPLE.porechop.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin)
minimap2 -a -x map-ont -Y -t 3 -R "@RG\tLB:$library\tID:$SAMPLE\tSM:$SAMPLE\tPL:ONT" /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.mmi porechop/trimmed/$SAMPLE.length.filtered.fastq.gz \
| sambamba sort -u --compression-level 6 --tmpdir=$WORK_DIR -t 3 -o bam/$SAMPLE.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin)
echo $line | cut -d" " -f 5-7 | sed 's/ /\t/g' > $SAMPLE.$TARGET.bed
cd mosdepth
mosdepth -t 1 --no-per-base --by ../$SAMPLE.$TARGET.bed --use-median --mapq 1 --fast-mode \
$SAMPLE.md ../bam/$SAMPLE.bam
cd ..
done < $batchname.metadata.tsv
rm porechop.adapters.py
module load samtools/1.21 vcflib/1.0.3 vt/0.57721

module load deepvariant/1.6.0 whatshap/2.3 annovar/2020-06-08
mkdir -p deepvariant annotation

while read -r line;
do
SAMPLE=$(echo $line | cut -d" " -f 1)
TARGET=$(echo $line | cut -d" " -f 2 | sed 's/_//' )
coverage_float=$(tail -n 1 mosdepth/$SAMPLE.md.mosdepth.summary.txt | cut -f 4)
coverage=$(perl -e "print int($coverage_float)")
echo "Coverage is" $coverage
if (( $coverage > 3 )); then
run_deepvariant --model_type ONT_R104 --num_shards 1 \
 --ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
 --regions $SAMPLE.$TARGET.bed \
 --reads bam/$SAMPLE.bam \
 --output_vcf deepvariant/$SAMPLE.dv.vcf.gz \
 --sample_name $SAMPLE \
 --intermediate_results_dir /lscratch/$SLURM_JOB_ID
bcftools norm --multiallelics -any --output-type u deepvariant/$SAMPLE.dv.vcf.gz \
| bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa  --output-type u --no-version - \
| bcftools annotate --set-id 'sml_%CHROM\:%POS%REF\>%ALT' -x ^FORMAT/GT,FORMAT/DP,FORMAT/VAF --output-type u --no-version \
| bcftools filter --include 'QUAL>0' --output-type u \
| bcftools norm -d exact --output-type z -o deepvariant/$SAMPLE.dv.filtered.vcf.gz
tabix -p vcf deepvariant/$SAMPLE.dv.filtered.vcf.gz
whatshap phase --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
 deepvariant/$SAMPLE.dv.filtered.vcf.gz bam/$SAMPLE.bam \
 | bgzip -f > deepvariant/$SAMPLE.dv.phased.vcf.gz
tabix -p vcf deepvariant/$SAMPLE.dv.phased.vcf.gz
rm deepvariant/$SAMPLE.dv.filtered.vcf.gz*
whatshap haplotag --output-threads 3 -o bam/$SAMPLE.haplotag.bam --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
  --ignore-read-groups --skip-missing-contigs deepvariant/$SAMPLE.dv.phased.vcf.gz bam/$SAMPLE.bam
samtools index -@ 3 bam/$SAMPLE.haplotag.bam
rm bam/$SAMPLE.bam bam/$SAMPLE.bam.bai
mv bam/$SAMPLE.haplotag.bam bam/$SAMPLE.bam
mv bam/$SAMPLE.haplotag.bam.bai bam/$SAMPLE.bam.bai
convert2annovar.pl -format vcf4old deepvariant/$SAMPLE.dv.phased.vcf.gz -includeinfo --outfile annotation/$SAMPLE.avinput
ver=hg38
table_annovar.pl annotation/$SAMPLE.avinput \
$ANNOVAR_DATA/$ver \
-buildver $ver \
-remove \
-out annotation/$SAMPLE.avinput \
--protocol refGeneWithVer \
-operation g \
--argument '-hgvs -splicing 50' \
--polish -nastring . \
--thread 1 \
--otherinfo
sed "1 s/Otherinfo1\tOtherinfo2\tOtherinfo3\tOtherinfo4\tOtherinfo5\tOtherinfo6\tOtherinfo7\tOtherinfo8\tOtherinfo9\tOtherinfo10/CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGT_FIELDS/" annotation/$SAMPLE.avinput."$ver"_multianno.txt \
	| awk -F"\t" 'BEGIN{OFS="\t"} {print $11,$12,$13,$14,$15,$16,$17,$18,$6,$7,$8,$9,$10,$19,$20}' > annotation/$SAMPLE.avinput."$ver"_multianno.txt.temp && mv annotation/$SAMPLE.avinput."$ver"_multianno.txt.temp annotation/$SAMPLE.avinput."$ver"_multianno.txt
rm annotation/$SAMPLE.avinput $SAMPLE.$TARGET.bed
else echo -e "####################\n$SAMPLE low coverage @ $coverage\n####################";
fi
done < $batchname.metadata.tsv
rm $batchname.metadata.bed

#try VEP instead, but how to select transcript? manually should be fine for amplicon approach.

#clean up directory
#rm -r fastq pod5_pass porechop

#rm -rf fastq porechop

##clair3 didn't call some SNP, even with the last line of --var_pct_full=1 --ref_pct_full=1 --var_pct_phasing=1
# mkdir -p clair3/$sample
# clair3 \
  # --bam_fn=bam/$sample.bam \
  # --ref_fn=/data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
  # --threads=8 \
  # --platform=ont \
  # --bed_fn=/data/OGL/resources/bed/OPN1_geneBody_LR-PCR.bed \
  # --model_path=/data/OGL/resources/clair3/ont_r1041/r1041_e82_400bps_sup_v420 \
  # --sample_name=$sample \
  # --output=clair3/$sample \
  # --var_pct_full=1 --ref_pct_full=1 --var_pct_phasing=1
  
##freebayes is not as good as deepvariant
# freebayes -f /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --limit-coverage 1000 --min-alternate-fraction 0.05 -C 3 -p 10 \
# --min-mapping-quality 0 --genotype-qualities --strict-vcf --use-mapping-quality --min-base-quality 5 \
# --targets /data/OGL/resources/bed/OPN1_geneBody_LR-PCR.bed \
# bam/$sample.bam \
# | vcffilter -f "QUAL > 10 & SAF > 0 & SAR > 0 & RPR > 0 & RPL > 0 & AO > 2 & DP > 5" \
# | bcftools +fill-tags - -Ou -- -t VAF \
# | bcftools norm --multiallelics -any --output-type u - \
# | vt decompose_blocksub -p -m -d 2 - \
# | bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u - \
# | bcftools annotate --set-id 'sm_%CHROM\:%POS%REF\>%ALT' -x ^INFO/AO,^FORMAT/GT,FORMAT/DP,FORMAT/VAF --output-type v --no-version \
# | bcftools norm -d exact --output-type z -o freebayes/$sample.small.vcf.gz
# tabix -p vcf freebayes/$sample.small.vcf.gz

#perhaps it's better to use the freebayes raw output to represent the variants, as bcftools norm may shorten the variants.

#rm annotation/$sample.avinput."$ver"_multianno.txt


#find LW/fastq -name "*fastq*" | parallel -j 10 -I% --max-args 1 --tag "minimap2 -a -Y -t 3 /data/OGL/resources/genomes/GRCh38/OPN1LW/GRCh38Decoy_OPN1LW.mmi % | sambamba sort -u --compression-level 6 --tmpdir=/lscratch/$SLURM_JOB_ID -t 3 -o %.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin)"

#find MW/fastq -name "*fastq*" | parallel -j 10 -I% --max-args 1 --tag "minimap2 -a -Y -t 3 /data/OGL/resources/genomes/GRCh38/OPN1MW/GRCh38Decoy_OPN1MW.mmi % | sambamba sort -u --compression-level 6 --tmpdir=/lscratch/$SLURM_JOB_ID -t 3 -o MW/bam/%.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin)"

#The minimap2 index were generated by: minimap2 -d genome.mmi  /fdb/igenomes/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa 
#The index is now moved to /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.mmi


