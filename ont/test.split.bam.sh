batchname="20250629_1324"
ml samtools R/4.3.0 whatshap clair3/1.0.10 annovar/2020-06-08
echo -e "clrNR_HetVar\tclr_HetVar\tdvNR_HetVar\tdv_HetVar" > temp.het.no.tsv
while read -r sample; do
if [ -e clair3/$sample.vcf.gz ]; then
	clrNR_HetVar=$(bcftools view -f 'PASS' -i 'GT="het"' --regions chrX:38268201-38285795,chrX:38286825-38328200 clair3/$sample.vcf.gz | grep -v "^#" | wc -l)
	clr_HetVar=$(bcftools view -f 'PASS' -i 'GT="het"' --regions chrX:38268201-38328200 clair3/$sample.vcf.gz | grep -v "^#" | wc -l)
	dvNR_HetVar=$(bcftools view -f 'PASS' -i 'GT="het"' -i 'QUAL > 0' --regions chrX:38268201-38285795,chrX:38286825-38328200 deepvariant/$sample.RPGR.dv.vcf.gz | grep -v "^#" | wc -l)
	dv_HetVar=$(bcftools view -f 'PASS' -i 'GT="het"' -i 'QUAL > 0' --regions chrX:38268201-38328200 deepvariant/$sample.RPGR.dv.vcf.gz | grep -v "^#" | wc -l)
	echo -e "$clrNR_HetVar\t$clr_HetVar\t$dvNR_HetVar\t$dv_HetVar" >> temp.het.no.tsv
	if (( clrNR_HetVar > 0 )); then
		bcftools view -f 'PASS' -i 'QUAL > 0' -Oz -o clair3/$sample.nr.vcf.gz --regions chrX:38268201-38285795,chrX:38286825-38328200 clair3/$sample.vcf.gz
		whatshap phase --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa  clair3/$sample.nr.vcf.gz bam/$sample.bam  | bgzip -f > clair3/$sample.nr.phased.vcf.gz
		tabix -f -p vcf clair3/$sample.nr.phased.vcf.gz
		whatshap haplotag --output-threads 3 -o bam/$sample.RPGR.haplotag.bam --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa   --ignore-read-groups --skip-missing-contigs --output-haplotag-list clair3/$sample.split.tsv clair3/$sample.nr.phased.vcf.gz bam/$sample.bam
		whatshap split --output-h1 bam/$sample.h1.bam --output-h2 bam/$sample.h2.bam bam/$sample.RPGR.haplotag.bam clair3/$sample.split.tsv
		samtools index bam/$sample.RPGR.haplotag.bam
		samtools index bam/$sample.h1.bam
		samtools index bam/$sample.h2.bam
		clair3 --bam_fn=bam/$sample.h1.bam --ref_fn=/data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --threads=8 --platform=ont --bed_fn=/data/OGL/resources/bed/RPGR.bed --model_path=/data/OGL/resources/clair3/ont_r1041/r1041_e82_400bps_sup_v420 --sample_name=$sample --ref_pct_full=1.0 --var_pct_full=0.5 --chunk_size=-1 --output=$sample/h1
		clair3 --bam_fn=bam/$sample.h2.bam --ref_fn=/data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --threads=8 --platform=ont --bed_fn=/data/OGL/resources/bed/RPGR.bed --model_path=/data/OGL/resources/clair3/ont_r1041/r1041_e82_400bps_sup_v420 --sample_name=$sample --ref_pct_full=1.0 --var_pct_full=0.5 --chunk_size=-1 --output=$sample/h2
		mv $sample/h1/merge_output.vcf.gz clair3/$sample.h1.vcf.gz
		mv $sample/h1/merge_output.vcf.gz.tbi clair3/$sample.h1.vcf.gz.tbi
		mv $sample/h2/merge_output.vcf.gz clair3/$sample.h2.vcf.gz
		mv $sample/h2/merge_output.vcf.gz.tbi clair3/$sample.h2.vcf.gz.tbi
		#rm bam/$sample.h1.bam* bam/$sample.h2.bam*
		#rm -r $sample
		bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u --no-version clair3/$sample.h1.vcf.gz \
			| bcftools annotate --set-id 'clr-h1_%CHROM\:%POS%REF\>%ALT' --output-type z -o clair3/$sample.h1.norm.vcf.gz
		tabix -f -p vcf clair3/$sample.h1.norm.vcf.gz
		bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u --no-version clair3/$sample.h2.vcf.gz \
			| bcftools annotate --set-id 'clr-h2_%CHROM\:%POS%REF\>%ALT' --output-type z -o clair3/$sample.h2.norm.vcf.gz
		tabix -f -p vcf clair3/$sample.h2.norm.vcf.gz
		bcftools concat --threads 1 -a --no-version \
			clair3/$sample.norm.vcf.gz clair3/$sample.h1.norm.vcf.gz clair3/$sample.h2.norm.vcf.gz deepvariant/$sample.RPGR.dv.phased.vcf.gz -Oz \
			-o vcf/$sample.clair.dv.vcf.gz
		tabix -f -p vcf vcf/$sample.clair.dv.vcf.gz
		elif (( clr_HetVar > 0 )); then
		bcftools view -f 'PASS' -i 'QUAL > 0' -Oz -o clair3/$sample.pass.q.vcf.gz --regions chrX:38268201-38328200 clair3/$sample.vcf.gz
		whatshap phase --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa  clair3/$sample.pass.q.vcf.gz bam/$sample.bam  | bgzip -f > clair3/$sample.pass.q.phased.vcf.gz
		tabix -f -p vcf clair3/$sample.pass.q.phased.vcf.gz
		whatshap haplotag --output-threads 3 -o bam/$sample.RPGR.haplotag.bam --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa   --ignore-read-groups --skip-missing-contigs --output-haplotag-list clair3/$sample.split.tsv clair3/$sample.pass.q.phased.vcf.gz bam/$sample.bam
		whatshap split --output-h1 bam/$sample.h1.bam --output-h2 bam/$sample.h2.bam bam/$sample.RPGR.haplotag.bam clair3/$sample.split.tsv
		samtools index bam/$sample.RPGR.haplotag.bam
		samtools index bam/$sample.h1.bam
		samtools index bam/$sample.h2.bam
		clair3 --bam_fn=bam/$sample.h1.bam --ref_fn=/data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --threads=8 --platform=ont --bed_fn=/data/OGL/resources/bed/RPGR.bed --model_path=/data/OGL/resources/clair3/ont_r1041/r1041_e82_400bps_sup_v420 --sample_name=$sample --ref_pct_full=1.0 --var_pct_full=0.5 --chunk_size=-1 --output=$sample/h1
		clair3 --bam_fn=bam/$sample.h2.bam --ref_fn=/data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --threads=8 --platform=ont --bed_fn=/data/OGL/resources/bed/RPGR.bed --model_path=/data/OGL/resources/clair3/ont_r1041/r1041_e82_400bps_sup_v420 --sample_name=$sample --ref_pct_full=1.0 --var_pct_full=0.5 --chunk_size=-1 --output=$sample/h2
		mv $sample/h1/merge_output.vcf.gz clair3/$sample.h1.vcf.gz
		mv $sample/h1/merge_output.vcf.gz.tbi clair3/$sample.h1.vcf.gz.tbi
		mv $sample/h2/merge_output.vcf.gz clair3/$sample.h2.vcf.gz
		mv $sample/h2/merge_output.vcf.gz.tbi clair3/$sample.h2.vcf.gz.tbi
		#rm bam/$sample.h1.bam* bam/$sample.h2.bam*
		#rm -r $sample
		bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u --no-version clair3/$sample.h1.vcf.gz \
			| bcftools annotate --set-id 'clr-h1.1_%CHROM\:%POS%REF\>%ALT' --output-type z -o clair3/$sample.h1.norm.vcf.gz
		tabix -f -p vcf clair3/$sample.h1.norm.vcf.gz
		bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u --no-version clair3/$sample.h2.vcf.gz \
			| bcftools annotate --set-id 'clr-h2.1_%CHROM\:%POS%REF\>%ALT' --output-type z -o clair3/$sample.h2.norm.vcf.gz
		tabix -f -p vcf clair3/$sample.h2.norm.vcf.gz
		bcftools concat --threads 1 -a --no-version \
			clair3/$sample.norm.vcf.gz clair3/$sample.h1.norm.vcf.gz clair3/$sample.h2.norm.vcf.gz deepvariant/$sample.RPGR.dv.phased.vcf.gz -Oz \
			-o vcf/$sample.clair.dv.vcf.gz
		tabix -f -p vcf vcf/$sample.clair.dv.vcf.gz
	else
		bcftools concat --threads 1 -a --no-version \
			clair3/$sample.norm.vcf.gz deepvariant/$sample.RPGR.dv.phased.vcf.gz -Oz \
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
	Rscript ~/git/NGS_genotype_calling/dev/annovar2excel.R annotation/$sample.avinput."$ver"_multianno.txt annotation/$sample.RPGR.clair.dv.anno.xlsx
	rm annotation/$sample.avinput*
else
echo "low RPGR coverage" >> temp.het.no.tsv
fi
done < manifest.csv
paste $batchname.coverage.summary.tsv temp.het.no.tsv > temp.$batchname.coverage.summary.tsv && mv temp.$batchname.coverage.summary.tsv $batchname.coverage.summary.tsv

run_deepvariant --model_type ONT_R104 --num_shards 1 \
--ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
--regions /data/OGL/resources/bed/RPGR.bed \
--reads bam/$sample.bam \
--output_vcf deepvariant/$sample.RPGR.dv.vcf.gz \
--sample_name $sample \
--intermediate_results_dir $WORK_DIR/dv/$sample
