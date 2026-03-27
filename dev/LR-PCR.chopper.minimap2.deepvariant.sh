#!/bin/bash
#SBATCH -c8
#SBATCH --mem=64g
#SBATCH --gres=lscratch:50
#SBATCH --time=2:0:0

#~20 min per sample, thus adjust the time accordingly. 
#fastq files have to be SampleID-LW and SampleID-MW. If missing either LW or MW, touch an empty file.

#set -e
library=$1 #QBplEx or OGL; will be part of @RG
WORK_DIR=/lscratch/$SLURM_JOB_ID

module load chopper/0.7.0 minimap2/2.26 sambamba/0.8.2

mkdir -p bam filteredFastq
case "${library^^}" in
	"QBPLEX")
	declare -i maxl=17000
	declare -i minl=200
	;;
	*)
	declare -i maxl=17000
	declare -i minl=12000
	;;
esac
echo $minl
if [ -e manifest.csv ];
then
	echo "manifest.csv provided, old bam file used for calling"
	echo "delete manifest.csv to re-align file"
else
	for fastq in fastq/*.fastq.gz; do
		sample=$(basename $fastq | sed "s/_reads.fastq.gz//" | sed "s/.fastq.gz//" | sed "s/-/_/");
		echo -e "$sample" >> manifest.csv
		zcat $fastq | chopper -q 10 --maxlength $maxl --minlength $minl | gzip > filteredFastq/$sample.filtered.fastq.gz
		minimap2 -a -x map-ont -Y -t 3 -R "@RG\tLB:$library\tID:$sample\tSM:$sample\tPL:ONT" /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.mmi filteredFastq/$sample.filtered.fastq.gz \
			| sambamba sort -u --compression-level 6 --tmpdir=$WORK_DIR -t 3 -o bam/$sample.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin);
	done
fi


#add coverage
#test freebayes

module load samtools/1.17 vcflib/1.0.3 vt/0.57721

module load freebayes/1.3.5 deepvariant/1.6.0 R/4.3.0 whatshap/1.1 annovar/2020-06-08
mkdir -p freebayes deepvariant
mkdir -p annotation

#changed -C 5 to 3 to account for low coverage Plasmid Express reads.
while read -r sample;
do
freebayes -f /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --max-complex-gap 80 -p 10 -C 3 -F 0.05 \
--genotype-qualities --strict-vcf --use-mapping-quality --min-base-quality 5 \
--targets /data/OGL/resources/bed/OPN1_e2-e5.bed \
bam/$sample.bam \
| vcffilter -f "QUAL > 10" \
| bcftools +fill-tags - -Ou -- -t VAF \
| bcftools norm --multiallelics -any --output-type u - \
| bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa  --output-type u --no-version - \
| bcftools annotate --set-id 'haplo_%CHROM\:%POS%REF\>%ALT' -x ^INFO/AO,^FORMAT/GT,FORMAT/DP,FORMAT/VAF --output-type u --no-version \
| bcftools norm -d exact --output-type z -o freebayes/$sample.haplo.vcf.gz
tabix -p vcf freebayes/$sample.haplo.vcf.gz
run_deepvariant --model_type ONT_R104 --num_shards 1 \
 --ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
 --regions /data/OGL/resources/bed/OPN1_geneBody_LR-PCR.bed \
 --reads bam/$sample.bam \
 --output_vcf deepvariant/$sample.dv.vcf.gz \
 --sample_name $sample \
 --intermediate_results_dir /lscratch/$SLURM_JOB_ID
bcftools norm --multiallelics -any --output-type u deepvariant/$sample.dv.vcf.gz \
| bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa  --output-type u --no-version - \
| bcftools annotate --set-id 'sml_%CHROM\:%POS%REF\>%ALT' -x ^FORMAT/GT,FORMAT/DP,FORMAT/VAF --output-type u --no-version \
| bcftools filter --include 'QUAL>0' --output-type u \
| bcftools norm -d exact --output-type z -o deepvariant/$sample.dv.filtered.vcf.gz
tabix -p vcf deepvariant/$sample.dv.filtered.vcf.gz
whatshap phase --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
 deepvariant/$sample.dv.filtered.vcf.gz bam/$sample.bam \
 | bgzip -f > deepvariant/$sample.dv.phased.vcf.gz
tabix -p vcf deepvariant/$sample.dv.phased.vcf.gz
rm deepvariant/$sample.dv.filtered.vcf.gz*
bcftools concat -a --rm-dups none --no-version \
deepvariant/$sample.dv.phased.vcf.gz freebayes/$sample.haplo.vcf.gz \
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
$sample annotation/$sample.avinput."$ver"_multianno.txt annotation/$sample.bcm.tsv annotation/$sample.bcm.xlsx
rm annotation/$sample.avinput annotation/$sample.avinput."$ver"_multianno.txt
done < manifest.csv

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


