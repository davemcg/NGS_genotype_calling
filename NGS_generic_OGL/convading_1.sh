module load $2
module load R/3.5

case $1 in
	1)
		perl ~/git/CoNVaDING/CoNVaDING.pl -mode StartWithBam \
			-inputDir sample_bam/$3 \
			-outputDir /lscratch/$SLURM_JOB_ID \
			-bed $4 \
			-useSampleAsControl \
			-controlsDir $5_male \
			-rmDup
		mkdir -p CoNVaDING/normalized_coverage_male
		cp /lscratch/$SLURM_JOB_ID/$3.b37.aligned.only.normalized.coverage.txt \
			CoNVaDING/normalized_coverage_male/$3.b37.aligned.only.normalized.coverage.txt
		Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/chrRD.R \
			CoNVaDING/normalized_coverage_male/$3.b37.aligned.only.normalized.coverage.txt \
			CoNVaDING/normalized_coverage_male/$3.chrRD.pdf \
			$6 \
			$7 \
			CoNVaDING/normalized_coverage_male/$3.abnormalChr.tsv \
			1
		touch $8
		;;
	2)
		perl ~/git/CoNVaDING/CoNVaDING.pl -mode StartWithBam \
			-inputDir sample_bam/$3 \
			-outputDir /lscratch/$SLURM_JOB_ID \
			-bed $4 \
			-useSampleAsControl \
			-controlsDir $5_female \
			-rmDup
		mkdir -p CoNVaDING/normalized_coverage_female
		cp /lscratch/$SLURM_JOB_ID/$3.b37.aligned.only.normalized.coverage.txt \
			CoNVaDING/normalized_coverage_female/$3.b37.aligned.only.normalized.coverage.txt
		Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/chrRD.R \
			CoNVaDING/normalized_coverage_female/$3.b37.aligned.only.normalized.coverage.txt \
			CoNVaDING/normalized_coverage_female/$3.chrRD.pdf \
			$6 \
			$7 \
			CoNVaDING/normalized_coverage_female/$3.abnormalChr.tsv \
			2
		touch $8
		;;
	*)
		perl ~/git/CoNVaDING/CoNVaDING.pl -mode StartWithBam \
			-inputDir sample_bam/$3 \
			-outputDir /lscratch/$SLURM_JOB_ID \
			-bed $4 \
			-useSampleAsControl \
			-controlsDir $5 \
			-rmDup
		mkdir -p CoNVaDING/normalized_coverage
		cp /lscratch/$SLURM_JOB_ID/$3.b37.aligned.only.normalized.coverage.txt \
			CoNVaDING/normalized_coverage/$3.b37.aligned.only.normalized.coverage.txt
		Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/chrRD.R \
			CoNVaDING/normalized_coverage/$3.b37.aligned.only.normalized.coverage.txt \
			CoNVaDING/normalized_coverage/$3.chrRD.pdf \
			$6 \
			$7 \
			CoNVaDING/normalized_coverage/$3.abnormalChr.tsv \
			0
		touch $8
		;;
esac
