#!/bin/bash
#SBATCH -c8
#SBATCH --mem=64g
#SBATCH --gres=lscratch:50
#SBATCH --time=2:0:0

#a few minutes per sample, thus adjust the time accordingly. 
#Run gel to make sure the amplicon is between 1.6-2.6 kb. If outside of this range, edit the max and min length accordingly
#set -e

library=$1 #QBplEx or OGL; will be part of @RG
WORK_DIR=/lscratch/$SLURM_JOB_ID

module load chopper/0.7.0 minimap2/2.26 sambamba/0.8.2

mkdir -p bam filteredFastq
case "${library^^}" in
	"QBPLEX")
	declare -i maxl=2600
	declare -i minl=200
	;;
	*)
	declare -i maxl=2600
	declare -i minl=1600
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
#rm -r filteredFastq
module load clair3/1.0.10 deepvariant/1.6.0 samtools/1.19 vcflib/1.0.3 annovar/2020-06-08 R/4.3.0

mkdir -p deepvariant clair3 vcf annotation

while read -r sample; do
	run_deepvariant --model_type ONT_R104 --num_shards 1 \
		--ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
		--regions /data/OGL/resources/bed/RPGR_ORF15.bed \
		--reads bam/$sample.bam \
		--output_vcf deepvariant/$sample.dv.vcf.gz \
		--sample_name $sample \
		--intermediate_results_dir $WORK_DIR/dv/$sample
	vcffilter -f "QUAL > 0" deepvariant/$sample.dv.vcf.gz \
		| bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u --no-version - \
		| bcftools annotate --set-id 'dv_%CHROM\:%POS%REF\>%ALT' -x FORMAT/PL --output-type z -o deepvariant/$sample.dv.filtered.vcf.gz
	tabix -p vcf deepvariant/$sample.dv.filtered.vcf.gz
	clair3 --bam_fn=bam/$sample.bam \
		--ref_fn=/data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
		--threads=8 \
		--platform=ont \
		--bed_fn=/data/OGL/resources/bed/RPGR_ORF15.bed \
		--model_path=/data/OGL/resources/clair3/ont_r1041/r1041_e82_400bps_sup_v420 \
		--sample_name=$sample \
		--ref_pct_full=1.0 --var_pct_full=0.5 \
		--chunk_size=-1 \
		--output=$WORK_DIR/clair3/$sample 
	###--ctg_name=chrX \ ## --bed_fn=/data/OGL/resources/bed/RPGR_ORF15.bed 
	##	--var_pct_full=1.0 \ --ref_pct_full=1.0 \
	##Will false negative increase if setting minimum read = 1 instead of 2 as in the default??? ##added var_pct_full & ref_pcf_full 8/6/24
	cp $WORK_DIR/clair3/$sample/merge_output.vcf.gz clair3/$sample.vcf.gz
	cp $WORK_DIR/clair3/$sample/merge_output.vcf.gz.tbi clair3/$sample.vcf.gz.tbi
	bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u --no-version clair3/$sample.vcf.gz \
		| bcftools annotate --set-id 'clr_%CHROM\:%POS%REF\>%ALT' --output-type z -o clair3/$sample.norm.vcf.gz
	tabix -p vcf clair3/$sample.norm.vcf.gz
	bcftools concat --threads 1 -a --no-version \
		clair3/$sample.norm.vcf.gz deepvariant/$sample.dv.filtered.vcf.gz -Oz \
		-o vcf/$sample.clair.dv.vcf.gz
	tabix -p vcf vcf/$sample.clair.dv.vcf.gz
	rm -rf $WORK_DIR/*
	convert2annovar.pl -format vcf4old vcf/$sample.clair.dv.vcf.gz -includeinfo --outfile annotation/$sample.avinput
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
	Rscript ~/git/NGS_genotype_calling/dev/annovar2excel.R annotation/$sample.avinput."$ver"_multianno.txt annotation/$sample.clair.dv.anno.xlsx
	rm annotation/$sample.avinput
	done < manifest.csv


#deepvariant higher Quality cut-off such as >0.5 or >1?? check the positive/negative cases.
#bcftools isec function #

#use bcftools normalize, VEP for proper right-alighment VEP vcf, bcftools, R 
 
#perhaps it's better to use the freebayes raw output to represent the variants, as bcftools norm may shorten the variants.

#rm annotation/$sample.avinput."$ver"_multianno.txt


#find LW/fastq -name "*fastq*" | parallel -j 10 -I% --max-args 1 --tag "minimap2 -a -Y -t 3 /data/OGL/resources/genomes/GRCh38/OPN1LW/GRCh38Decoy_OPN1LW.mmi % | sambamba sort -u --compression-level 6 --tmpdir=/lscratch/$SLURM_JOB_ID -t 3 -o %.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin)"

#find MW/fastq -name "*fastq*" | parallel -j 10 -I% --max-args 1 --tag "minimap2 -a -Y -t 3 /data/OGL/resources/genomes/GRCh38/OPN1MW/GRCh38Decoy_OPN1MW.mmi % | sambamba sort -u --compression-level 6 --tmpdir=/lscratch/$SLURM_JOB_ID -t 3 -o MW/bam/%.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin)"

#The minimap2 index were generated by: minimap2 -d genome.mmi  /fdb/igenomes/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa 
#The index is now moved to /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.mmi

