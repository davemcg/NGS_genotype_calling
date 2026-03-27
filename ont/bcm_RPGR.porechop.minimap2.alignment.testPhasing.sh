#!/bin/bash
#SBATCH -c8
#SBATCH --mem=64g
#SBATCH --gres=lscratch:50
#SBATCH --time=4:0:0

#Could submit with dependency then start automatically after the previous step completes. 
#sbatch --dependency=afterok:11254323 ~/git/NGS_genotype_calling/dev/bcm_RPGR.porechop.minimap2.alignment.sh OGL
#~20 min per sample, thus adjust the time accordingly. 
#fastq files have to be SampleID-LW and SampleID-MW. If missing either LW or MW, touch an empty file.

set -e
library=$1 #QBplEx or OGL; will be part of @RG
batchname=$(date '+%Y%m%d_%H%M') 
WORK_DIR=/lscratch/$SLURM_JOB_ID

module load porechop/0.2.4 minimap2/2.26 sambamba/0.8.2 mosdepth

mkdir -p bam porechop mosdepth
#code below not used
#module load chopper/0.7.0
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


if [ -e manifest.csv ]; then
	echo "manifest.csv provided, old bam file used for calling"
	echo "delete manifest.csv to re-align file"
else
	for fastq in fastq/*.fastq.gz; do
		sample=$(basename $fastq | sed "s/_reads.fastq.gz//" | sed "s/.fastq.gz//" | sed "s/-/_/");
		echo -e "$sample" >> manifest.csv
	done
	export PORECHOP_ADAPTERS=/data/OGL/resources/porechop/adaptersRPGRr_BCM_E6R.py
	while read -r sample; do
		mkdir -p porechop/$sample
		zcat fastq/$sample.fastq.gz | cat - ~/git/NGS_genotype_calling/ont/bcmORF15.fastq | gzip > porechop/$sample/input.fastq.gz
		porechop --discard_unassigned  --discard_middle --extra_end_trim 0 -i porechop/$sample/input.fastq.gz -b porechop/$sample
	done < manifest.csv
	export PORECHOP_ADAPTERS=/data/OGL/resources/porechop/adaptersBCMf.py
	while read -r sample; do
		porechop --discard_unassigned  --discard_middle --extra_end_trim 0 --adapter_threshold 60 -i porechop/$sample/BCE6R.fastq.gz -b porechop/$sample/bcm
	done < manifest.csv
	export PORECHOP_ADAPTERS=/data/OGL/resources/porechop/adaptersRPGRf.py
	while read -r sample; do
	porechop --discard_unassigned  --discard_middle --extra_end_trim 0  --adapter_threshold 60 -i porechop/$sample/BCRPGRr.fastq.gz -b porechop/$sample/RPGR
	if (( $(zcat porechop/$sample/RPGR/BCRPGRf.fastq.gz | wc -l) > 4 )); then zcat porechop/$sample/RPGR/BCRPGRf.fastq.gz | head -n -4 | gzip > /lscratch/$SLURM_JOB_ID/temp.$sample.fastq.gz && mv /lscratch/$SLURM_JOB_ID/temp.$sample.fastq.gz porechop/$sample/RPGR/BCRPGRf.fastq.gz; else echo -e "$sample has 0 RPGR full-length amplicon."; fi
	if (( $(zcat porechop/$sample/bcm/BCLW.fastq.gz | wc -l) > 4 )); then zcat porechop/$sample/bcm/BCLW.fastq.gz | head -n -4 | gzip > /lscratch/$SLURM_JOB_ID/temp.$sample.fastq.gz && mv /lscratch/$SLURM_JOB_ID/temp.$sample.fastq.gz porechop/$sample/bcm/BCLW.fastq.gz; else echo -e "$sample has 0 LW full-length amplicon."; fi
	if (( $(zcat porechop/$sample/bcm/BCMW.fastq.gz | wc -l) > 4 )); then zcat porechop/$sample/bcm/BCMW.fastq.gz | head -n -4 | gzip > /lscratch/$SLURM_JOB_ID/temp.$sample.fastq.gz && mv /lscratch/$SLURM_JOB_ID/temp.$sample.fastq.gz porechop/$sample/bcm/BCMW.fastq.gz; else echo -e "$sample has 0 MW full-length amplicon."; fi
	minimap2 -a -x map-ont -Y -t 3 -R "@RG\tLB:$library\tID:$sample\tSM:$sample\tPL:ONT" /data/OGL/resources/genomes/GRCh38/OPN1LW/GRCh38Decoy_OPN1LW.mmi porechop/$sample/RPGR/BCRPGRf.fastq.gz \
		| sambamba sort -u --compression-level 6 --tmpdir=$WORK_DIR -t 3 -o bam/$sample.RPGR.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin);
	minimap2 -a -x map-ont -Y -t 3 -R "@RG\tLB:$library\tID:$sample\tSM:$sample\tPL:ONT" /data/OGL/resources/genomes/GRCh38/OPN1LW/GRCh38Decoy_OPN1LW.mmi porechop/$sample/bcm/BCLW.fastq.gz \
		| sambamba sort -u --compression-level 6 --tmpdir=$WORK_DIR -t 3 -o bam/$sample.LW.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin);
	minimap2 -a -x map-ont -Y -t 3 -R "@RG\tLB:$library\tID:$sample\tSM:$sample\tPL:ONT" /data/OGL/resources/genomes/GRCh38/OPN1MW/GRCh38Decoy_OPN1MW.mmi porechop/$sample/bcm/BCMW.fastq.gz \
		| sambamba sort -u --compression-level 6 --tmpdir=$WORK_DIR -t 3 -o bam/$sample.MW.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin);
	sambamba merge -t 3 -l 6 bam/$sample.bam bam/$sample.RPGR.bam bam/$sample.LW.bam bam/$sample.MW.bam
	#sambamba merge -t 3 -l 6 bam/$sample.RPGR-LW.bam bam/$sample.RPGR.bam bam/$sample.LW.bam
	#rm bam/$sample.RPGR.bam* bam/$sample.LW.bam* bam/$sample.MW.bam*
	cd mosdepth
	mosdepth -t 1 --no-per-base --by /data/OGL/resources/bed/bcmRPGR.bed --use-median --mapq 1 --fast-mode \
	$sample.md ../bam/$sample.bam
	cd ..
	done < manifest.csv
fi

echo -e "****** Step 1: Alignment completed."

module load R/4.3.0
echo -e "Sample\tRPGR\tOPN1LW\tOPN1MW\tMax_OPN1" > $batchname.coverage.summary.tsv
while read -r sample; do
	Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/bcmRPGRcoverageWide.R $sample mosdepth/$sample.md.regions.bed.gz mosdepth/$sample.coverage.wide.tsv
	tail -n 1 mosdepth/$sample.coverage.wide.tsv >> $batchname.coverage.summary.tsv
done < manifest.csv

#zcat mosdepth/$sample.md.regions.bed.gz | sed "s/^/$sample\t/" >> mosdepth/coverage.summary.tsv
#test freebayes

module load samtools/1.21 vt/0.57721 #vcflib/1.0.3
module load freebayes/1.3.5 deepvariant/1.6.0 whatshap/2.3 annovar/2020-06-08
mkdir -p freebayes deepvariant vcf annotation

#changed -C 5 to 3 to account for low coverage Plasmid Express reads.
while read -r sample; do
	mkdir -p bam/$sample
	if (( $(grep $sample $batchname.coverage.summary.tsv | cut -f 4) > 3 )); then
	#freebayes for haplotype-tagging of MW, vcf not annotated by annovar
	freebayes -f /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa -p 6 -C 5 -F 0.1 \
		--genotype-qualities --strict-vcf --use-mapping-quality --min-base-quality 5 \
		--region chrX:154181866-154196279 \
		bam/$sample.bam \
		| bcftools view -i "QUAL > 10" -Ou - \
		| bcftools +fill-tags - -Ou -- -t VAF \
		| bcftools norm --multiallelics -any --output-type u - \
		| bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa  --output-type u --no-version - \
		| bcftools annotate --set-id 'haplo_%CHROM\:%POS%REF\>%ALT' -x ^INFO/AO,^FORMAT/GT,FORMAT/DP,FORMAT/VAF --output-type u --no-version \
		| bcftools norm -d exact --output-type z -o freebayes/$sample.MW.vcf.gz
	tabix -f -p vcf freebayes/$sample.MW.vcf.gz
	whatshap polyphase --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
		freebayes/$sample.MW.vcf.gz bam/$sample.bam --ploidy 6 \
		| bgzip -f > freebayes/$sample.MW.phased.vcf.gz
	tabix -f -p vcf freebayes/$sample.MW.phased.vcf.gz
	whatshap haplotag --output-threads 3 --ploidy 6 -o bam/$sample/$sample.MW.haplotag.bam --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
		--ignore-read-groups --skip-missing-contigs freebayes/$sample.MW.phased.vcf.gz bam/$sample.MW.bam
	samtools index -@ 3 bam/$sample/$sample.MW.haplotag.bam
	else
	echo -e "#################### $sample low MW coverage @ $coverage ####################"
	fi
done < manifest.csv

echo -e "****** Step 2: MW haplotag - polyphase 6 completed."

while read -r sample; do
coverage=$(grep $sample $batchname.coverage.summary.tsv | cut -f 5)
echo -e "Coverage of $sample is $coverage"
if (( $coverage > 3 )); then
#variant calls
freebayes -f /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --max-complex-gap 80 -p 6 -C 3 -F 0.05 \
--genotype-qualities --strict-vcf --use-mapping-quality --min-base-quality 5 \
--targets /data/OGL/resources/bed/OPN1_e2-e5.bed  \
bam/$sample.bam \
| bcftools view -i "QUAL > 10" -Ou - \
| bcftools +fill-tags - -Ou -- -t VAF \
| bcftools norm --multiallelics -any --output-type u - \
| bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa  --output-type u --no-version - \
| bcftools annotate --set-id 'haplo_%CHROM\:%POS%REF\>%ALT' -x ^INFO/AO,^FORMAT/GT,FORMAT/DP,FORMAT/VAF --output-type u --no-version \
| bcftools norm -d exact --output-type z -o freebayes/$sample.haplo.vcf.gz
tabix -f -p vcf freebayes/$sample.haplo.vcf.gz
run_deepvariant --model_type ONT_R104 --num_shards 1 \
 --ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
 --regions /data/OGL/resources/bed/OPN1_geneBody_LR-PCR.bed \
 --reads bam/$sample.bam \
 --output_vcf deepvariant/$sample.OPN1.dv.vcf.gz \
 --sample_name $sample \
 --intermediate_results_dir /lscratch/$SLURM_JOB_ID
bcftools norm --multiallelics -any --output-type u deepvariant/$sample.OPN1.dv.vcf.gz \
| bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa  --output-type u --no-version - \
| bcftools annotate --set-id 'sml_%CHROM\:%POS%REF\>%ALT' -x ^FORMAT/GT,FORMAT/DP,FORMAT/VAF --output-type u --no-version \
| bcftools filter --include 'QUAL>0' --output-type u \
| bcftools norm -d exact --output-type z -o deepvariant/$sample.OPN1.dv.filtered.vcf.gz
tabix -f -p vcf deepvariant/$sample.OPN1.dv.filtered.vcf.gz
whatshap phase --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
 deepvariant/$sample.OPN1.dv.filtered.vcf.gz bam/$sample.bam \
 | bgzip -f > deepvariant/$sample.OPN1.dv.phased.vcf.gz
tabix -f -p vcf deepvariant/$sample.OPN1.dv.phased.vcf.gz
rm deepvariant/$sample.OPN1.dv.filtered.vcf.gz*
whatshap haplotag --output-threads 3 -o bam/$sample/$sample.LW.haplotag.bam --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
  --ignore-read-groups --skip-missing-contigs deepvariant/$sample.OPN1.dv.phased.vcf.gz bam/$sample.LW.bam
samtools index -@ 3 bam/$sample/$sample.LW.haplotag.bam
bcftools concat -a --rm-dups none --no-version \
deepvariant/$sample.OPN1.dv.phased.vcf.gz freebayes/$sample.haplo.vcf.gz \
-Oz -o vcf/$sample.OPN1.vcf.gz
tabix -f -p vcf vcf/$sample.OPN1.vcf.gz
convert2annovar.pl -format vcf4old vcf/$sample.OPN1.vcf.gz -includeinfo --outfile annotation/$sample.avinput
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
rm annotation/$sample.avinput*
else echo -e "#################### $sample low LW and MW coverage @ $coverage ####################";fi
done < manifest.csv

echo -e "****** Step 3: LW/MW variant call completed."

module load clair3/1.0.10

mkdir -p deepvariant clair3 vcf annotation

while read -r sample; do
	mkdir -p bam/$sample
	coverage=$(grep $sample $batchname.coverage.summary.tsv | cut -f 2)
	echo -e "RPGR coverage of $sample is $coverage."
	if (( $coverage > 3 )); then
	run_deepvariant --model_type ONT_R104 --num_shards 1 \
	--ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
	--regions /data/OGL/resources/bed/RPGR.bed \
	--reads bam/$sample.bam \
	--output_vcf deepvariant/$sample.RPGR.dv.vcf.gz \
	--sample_name $sample \
	--intermediate_results_dir $WORK_DIR/dv/$sample
	bcftools view -i "QUAL > 0" -Ou deepvariant/$sample.RPGR.dv.vcf.gz \
	| bcftools norm --multiallelics -any --output-type u - \
	| bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u --no-version - \
	| bcftools annotate --set-id 'dv_%CHROM\:%POS%REF\>%ALT' -x FORMAT/PL --output-type z -o deepvariant/$sample.RPGR.dv.filtered.vcf.gz
	tabix -f -p vcf deepvariant/$sample.RPGR.dv.filtered.vcf.gz
	whatshap phase --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
	 deepvariant/$sample.RPGR.dv.filtered.vcf.gz bam/$sample.bam \
	 | bgzip -f > deepvariant/$sample.RPGR.dv.phased.vcf.gz
	tabix -f -p vcf deepvariant/$sample.RPGR.dv.phased.vcf.gz
	rm deepvariant/$sample.RPGR.dv.filtered.vcf.gz*
	whatshap haplotag --output-threads 3 -o bam/$sample/$sample.RPGR.haplotag.bam --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
	  --ignore-read-groups --skip-missing-contigs deepvariant/$sample.RPGR.dv.phased.vcf.gz bam/$sample.RPGR.bam
	samtools index -@ 3 bam/$sample/$sample.RPGR.haplotag.bam
	clair3 --bam_fn=bam/$sample.bam \
	--ref_fn=/data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
	--threads=8 \
	--platform=ont \
	--bed_fn=/data/OGL/resources/bed/RPGR.bed \
	--model_path=/data/OGL/resources/clair3/ont_r1041/r1041_e82_400bps_sup_v420 \
	--sample_name=$sample \
	--ref_pct_full=1.0 --var_pct_full=0.5 \
	--chunk_size=-1 \
	--output=$WORK_DIR/clair3/$sample 
	###--ctg_name=chrX \ ## --bed_fn=/data/OGL/resources/bed/RPGR_ORF15.bed 
	###--var_pct_full=1.0 \ --ref_pct_full=1.0 \
	##Will false negative increase if setting minimum read = 1 instead of 2 as in the default??? ##added var_pct_full & ref_pcf_full 8/6/24
	bcftools norm --multiallelics -any --output-type u $WORK_DIR/clair3/$sample/merge_output.vcf.gz \
	| bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u --no-version - \ 
	| bcftools filter --include 'FORMAT/AD[0:1]>2' --output-type u - \
	| bcftools annotate --set-id 'clr_%CHROM\:%POS%REF\>%ALT' --output-type z -o clair3/$sample.vcf.gz
	tabix -f -p vcf clair3/$sample.vcf.gz
	rm -rf $WORK_DIR/*
	else echo -e "####################$sample low RPGR coverage @ $coverage ####################";
	fi
done < manifest.csv
##deepvariant RefCall variants could have QUAL score of 0 or a value like 0.2, but genotype will be ./.
echo -e "****** Step 4: RPGR variant call Part A completed."

#whatshap haplog merge
while read -r sample; do
bamfile=""
if [ -e bam/$sample/$sample.RPGR.haplotag.bam ]; then 
bamfile+=" bam/$sample/$sample.RPGR.haplotag.bam"
fi
if [ -e bam/$sample/$sample.LW.haplotag.bam ]; then
bamfile+=" bam/$sample/$sample.LW.haplotag.bam"
fi
if [ -e bam/$sample/$sample.MW.haplotag.bam ]; then
bamfile+=" bam/$sample/$sample.MW.haplotag.bam"
fi
if [ $bamfile=="" ]; then
echo -e "#################### $sample low coverage for RPGR, LW, & MW ####################"
else
sambamba merge -t 3 -l 6 bam/$sample.haplotag.bam $bamfile
fi
done < manifest.csv

echo -e "****** Step 5: Deepvariant Haplotagged RPGR/LW + freebayes MW bams combined."

#Phase RPGR locus using het variants in the non-rep regions as 1st choide,
#then phase using all het variants if 0 in non-rep regions, so to maximize sensitivity.
echo -e "clrNR_HetVar\tclr_HetVar\tdvNR_HetVar\tdv_HetVar\tclrNR_HetVar_noHmPlm\tclr_HetVar_noHmPlm\tdvNR_HetVar_noHmPlm\tdv_HetVar_noHmPlm" > temp.het.no.tsv
while read -r sample; do
if [ -e clair3/$sample.vcf.gz ]; then
	clrNR_HetVar=$(bcftools view -f 'PASS' -i 'GT="het"' --regions chrX:38268201-38285795,chrX:38286825-38328200 clair3/$sample.vcf.gz | grep -v "^#" | wc -l)
	clr_HetVar=$(bcftools view -f 'PASS' -i 'GT="het"' --regions chrX:38268201-38328200 clair3/$sample.vcf.gz | grep -v "^#" | wc -l)
	dvNR_HetVar=$(bcftools view -f 'PASS' -i 'GT="het"' --regions chrX:38268201-38285795,chrX:38286825-38328200 deepvariant/$sample.RPGR.dv.vcf.gz | grep -v "^#" | wc -l)
	dv_HetVar=$(bcftools view -f 'PASS' -i 'GT="het"' --regions chrX:38268201-38328200 deepvariant/$sample.RPGR.dv.vcf.gz | grep -v "^#" | wc -l)
	
	clrNR_HetVar_noHmPlm=$(bcftools view -f 'PASS' -i 'GT="het"' --regions chrX:38268201-38285413,chrX:38285424-38285795,chrX:38286825-38298268,chrX:38298283-38299738,chrX:38299749-38328200 clair3/$sample.vcf.gz | grep -v "^#" | wc -l)
	clr_HetVar_noHmPlm=$(bcftools view -f 'PASS' -i 'GT="het"' --regions chrX:38268201-38285413,chrX:38285424-38298268,chrX:38298283-38299738,chrX:38299749-38328200 clair3/$sample.vcf.gz | grep -v "^#" | wc -l)
	dvNR_HetVar_noHmPlm=$(bcftools view -f 'PASS' -i 'GT="het"' --regions chrX:38268201-38285413,chrX:38285424-38285795,chrX:38286825-38298268,chrX:38298283-38299738,chrX:38299749-38328200 deepvariant/$sample.RPGR.dv.vcf.gz | grep -v "^#" | wc -l)
	dv_HetVar_noHmPlm=$(bcftools view -f 'PASS' -i 'GT="het"' --regions chrX:38268201-38285413,chrX:38285424-38298268,chrX:38298283-38299738,chrX:38299749-38328200 deepvariant/$sample.RPGR.dv.vcf.gz | grep -v "^#" | wc -l)
	echo -e "$clrNR_HetVar\t$clr_HetVar\t$dvNR_HetVar\t$dv_HetVar\t$clrNR_HetVar_noHmPlm\t$clr_HetVar_noHmPlm\t$dvNR_HetVar_noHmPlm\t$dv_HetVar_noHmPlm" >> temp.het.no.tsv
	if (( clrNR_HetVar_noHmPlm > 0 )); then
		bcftools view -f 'PASS' -i 'QUAL > 0' -Oz -o clair3/$sample.nr.vcf.gz --regions chrX:38268201-38285413,chrX:38285424-38285795,chrX:38286825-38298268,chrX:38298283-38299738,chrX:38299749-38328200 clair3/$sample.vcf.gz
		whatshap phase --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa  clair3/$sample.nr.vcf.gz bam/$sample.bam  | bgzip -f > clair3/$sample.nr.phased.vcf.gz
		tabix -f -p vcf clair3/$sample.nr.phased.vcf.gz
		whatshap haplotag --output-threads 3 -o bam/$sample.RPGR.clr.haplotag.bam --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa   --ignore-read-groups --skip-missing-contigs --output-haplotag-list clair3/$sample.split.tsv clair3/$sample.nr.phased.vcf.gz bam/$sample.bam
		whatshap split --output-h1 bam/$sample.h1.bam --output-h2 bam/$sample.h2.bam bam/$sample.RPGR.clr.haplotag.bam clair3/$sample.split.tsv
		samtools index bam/$sample.RPGR.clr.haplotag.bam
		samtools index bam/$sample.h1.bam
		samtools index bam/$sample.h2.bam
		clair3 --bam_fn=bam/$sample.h1.bam --ref_fn=/data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --threads=8 --platform=ont --bed_fn=/data/OGL/resources/bed/RPGR.bed --model_path=/data/OGL/resources/clair3/ont_r1041/r1041_e82_400bps_sup_v420 --sample_name=$sample --ref_pct_full=1.0 --var_pct_full=0.5 --chunk_size=-1 --output=clair3/$sample/h1
		clair3 --bam_fn=bam/$sample.h2.bam --ref_fn=/data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --threads=8 --platform=ont --bed_fn=/data/OGL/resources/bed/RPGR.bed --model_path=/data/OGL/resources/clair3/ont_r1041/r1041_e82_400bps_sup_v420 --sample_name=$sample --ref_pct_full=1.0 --var_pct_full=0.5 --chunk_size=-1 --output=clair3/$sample/h2
		rm bam/$sample.h1.bam* bam/$sample.h2.bam*
		bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u --no-version clair3/$sample/h1/merge_output.vcf.gz \
			| bcftools annotate --set-id 'clr-h1_%CHROM\:%POS%REF\>%ALT' --output-type z -o clair3/$sample.h1.vcf.gz
		tabix -f -p vcf clair3/$sample.h1.vcf.gz
		bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u --no-version clair3/$sample/h2/merge_output.vcf.gz \
			| bcftools annotate --set-id 'clr-h2_%CHROM\:%POS%REF\>%ALT' --output-type z -o clair3/$sample.h2.vcf.gz
		tabix -f -p vcf clair3/$sample.h2.vcf.gz
		bcftools concat --threads 1 -a --no-version \
			clair3/$sample.vcf.gz clair3/$sample.h1.vcf.gz clair3/$sample.h2.vcf.gz deepvariant/$sample.RPGR.dv.phased.vcf.gz -Oz \
			-o vcf/$sample.clair.dv.vcf.gz
		tabix -f -p vcf vcf/$sample.clair.dv.vcf.gz
		rm -r clair3/$sample
		elif (( clr_HetVar_noHmPlm > 0 )); then
		bcftools view -f 'PASS' -i 'QUAL > 0' -Oz -o clair3/$sample.pass.q.vcf.gz --regions chrX:38268201-38285413,chrX:38285424-38298268,chrX:38298283-38299738,chrX:38299749-38328200 clair3/$sample.vcf.gz
		whatshap phase --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa  clair3/$sample.pass.q.vcf.gz bam/$sample.bam  | bgzip -f > clair3/$sample.pass.q.phased.vcf.gz
		tabix -f -p vcf clair3/$sample.pass.q.phased.vcf.gz
		whatshap haplotag --output-threads 3 -o bam/$sample.RPGR.clr.haplotag.bam --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa   --ignore-read-groups --skip-missing-contigs --output-haplotag-list clair3/$sample.split.tsv clair3/$sample.pass.q.phased.vcf.gz bam/$sample.bam
		whatshap split --output-h1 bam/$sample.h1.bam --output-h2 bam/$sample.h2.bam bam/$sample.RPGR.clr.haplotag.bam clair3/$sample.split.tsv
		samtools index bam/$sample.RPGR.clr.haplotag.bam
		samtools index bam/$sample.h1.bam
		samtools index bam/$sample.h2.bam
		clair3 --bam_fn=bam/$sample.h1.bam --ref_fn=/data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --threads=8 --platform=ont --bed_fn=/data/OGL/resources/bed/RPGR.bed --model_path=/data/OGL/resources/clair3/ont_r1041/r1041_e82_400bps_sup_v420 --sample_name=$sample --ref_pct_full=1.0 --var_pct_full=0.5 --chunk_size=-1 --output=clair3/$sample/h1
		clair3 --bam_fn=bam/$sample.h2.bam --ref_fn=/data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --threads=8 --platform=ont --bed_fn=/data/OGL/resources/bed/RPGR.bed --model_path=/data/OGL/resources/clair3/ont_r1041/r1041_e82_400bps_sup_v420 --sample_name=$sample --ref_pct_full=1.0 --var_pct_full=0.5 --chunk_size=-1 --output=clair3/$sample/h2
		rm bam/$sample.h1.bam* bam/$sample.h2.bam*
		bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u --no-version clair3/$sample/h1/merge_output.vcf.gz \
			| bcftools annotate --set-id 'clr-h1.1_%CHROM\:%POS%REF\>%ALT' --output-type z -o clair3/$sample.h1.vcf.gz
		tabix -f -p vcf clair3/$sample.h1.vcf.gz
		bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u --no-version clair3/$sample/h2/merge_output.vcf.gz \
			| bcftools annotate --set-id 'clr-h2.1_%CHROM\:%POS%REF\>%ALT' --output-type z -o clair3/$sample.h2.vcf.gz
		tabix -f -p vcf clair3/$sample.h2.vcf.gz
		bcftools concat --threads 1 -a --no-version \
			clair3/$sample.vcf.gz clair3/$sample.h1.vcf.gz clair3/$sample.h2.vcf.gz deepvariant/$sample.RPGR.dv.phased.vcf.gz -Oz \
			-o vcf/$sample.clair.dv.vcf.gz
		tabix -f -p vcf vcf/$sample.clair.dv.vcf.gz
		rm -r clair3/$sample
	else
		bcftools concat --threads 1 -a --no-version \
			clair3/$sample.vcf.gz deepvariant/$sample.RPGR.dv.phased.vcf.gz -Oz \
			-o vcf/$sample.clair.dv.vcf.gz
		tabix -p vcf vcf/$sample.clair.dv.vcf.gz
	fi
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
	Rscript ~/git/NGS_genotype_calling/ont/annovar2excel.R annotation/$sample.avinput."$ver"_multianno.txt annotation/$sample.RPGR.clair.dv.anno.xlsx
	rm annotation/$sample.avinput*
else
echo -e "low RPGR coverage\tlow RPGR coverage\tlow RPGR coverage\tlow RPGR coverage\tlow RPGR coverage\tlow RPGR coverage\tlow RPGR coverage\tlow RPGR coverage" >> temp.het.no.tsv
fi
done < manifest.csv
paste $batchname.coverage.summary.tsv temp.het.no.tsv > temp.$batchname.coverage.summary.tsv && mv temp.$batchname.coverage.summary.tsv $batchname.coverage.summary.tsv

echo -e "****** Step 6: Seperate Hap calls and RPGR annotation completed."

#clean up directory

while read -r sample; do
rm -f bam/$sample.RPGR-LW.haplotag.bam* bam/$sample.MW.haplotag.bam*
rm -f bam/$sample.LW.bam* bam/$sample.RPGR-LW.bam* bam/$sample.MW.bam* bam/$sample.RPGR.bam*
rm -rf bam/$sample
rm -f deepvariant/*.html
rm -rf pod5_pass porechop
rm -f temp.het.no.tsv
done < manifest.csv

#nextstep: change to VEP in the next version for better HGVS compliance.
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
# tabix -f -p vcf freebayes/$sample.small.vcf.gz

#perhaps it's better to use the freebayes raw output to represent the variants, as bcftools norm may shorten the variants.

#rm annotation/$sample.avinput."$ver"_multianno.txt


#find LW/fastq -name "*fastq*" | parallel -j 10 -I% --max-args 1 --tag "minimap2 -a -Y -t 3 /data/OGL/resources/genomes/GRCh38/OPN1LW/GRCh38Decoy_OPN1LW.mmi % | sambamba sort -u --compression-level 6 --tmpdir=/lscratch/$SLURM_JOB_ID -t 3 -o %.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin)"

#find MW/fastq -name "*fastq*" | parallel -j 10 -I% --max-args 1 --tag "minimap2 -a -Y -t 3 /data/OGL/resources/genomes/GRCh38/OPN1MW/GRCh38Decoy_OPN1MW.mmi % | sambamba sort -u --compression-level 6 --tmpdir=/lscratch/$SLURM_JOB_ID -t 3 -o MW/bam/%.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin)"

#The minimap2 index were generated by: minimap2 -d genome.mmi  /fdb/igenomes/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa 
#The index is now moved to /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.mmi


