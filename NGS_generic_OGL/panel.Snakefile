from os.path import join
import sys
import datetime
#import os.path # for checking whether a file exist. 7/21/19
#from os import path # for checking whether a file exist. 7/21/19

# the 'metadata_file' if provided, is a csv with three columns
# the first is the sample name (e.g. Patient001)
# the second is the name of the fastq or bam associated with the sample
# the third is the read group you want bwa to use
# 	example: '@RG\\tID:Lineagen_41001412010527\\tSM:Lineagen_41001412010527\\tPL:ILLUMINA'
# a header isn't required, but if it is included it MUST start with #:
# 	#Sample,File
# you can have multiple lines per sample
# most samples are paired end, so there will be at least two files per sample
# often you have a sample sequenced across multiple lanes/machines, so you can have
# upwards of a dozen files for a single sample
SAMPLE_LANEFILE = dict()
LANEFILE_READGROUP = dict()
SAMPLE_SEX = dict()
metadata = open(config['metadata_file'])
for line in metadata:
	read_group = line.split(',')[2][:-1]
	lane_file = line.split(',')[1]
	sample = line.split(',')[0]
	# skip header
	if line.startswith("#"):
		continue
	if sample not in SAMPLE_LANEFILE:
		SAMPLE_LANEFILE[sample] = [lane_file]
		SAMPLE_SEX[sample] = 0
	else:
		old_lane_file = SAMPLE_LANEFILE[sample]
		old_lane_file.append(lane_file)
		SAMPLE_LANEFILE[sample] = old_lane_file
	LANEFILE_READGROUP[lane_file] = [read_group]
metadata.close()
for i in SAMPLE_LANEFILE:
	print (i, SAMPLE_LANEFILE[i], len(SAMPLE_LANEFILE[i]))
#default sample sex is 0 as setting above
#if config['ped'] != '':
if config['ped']:
	with open(config['ped']) as PED_file:
		for line in PED_file:
			if line.startswith("#"):
				continue
			else:
				sample = line.split('\t')[1]
				SAMPLE_SEX[sample] = line.split("\t")[4]
#				SAMPLE_SEX[sample] = [line.split("\t")[4]]
# try if confi['ped']: empty string returns false.
#if line.strip(): empty line returns false.
# for i in SAMPLE_SEX:
#  	print (i, SAMPLE_SEX[i])

if config['analysis_batch_name'] == 'YYYYMMDD':
	currentDT = datetime.datetime.now()
	config['analysis_batch_name'] = currentDT.strftime("%Y%m%d")

def rg(wildcards):
	# returns the read group given in the config['metadata_file']
	lane_file = str(wildcards)
	rg_out = str(LANEFILE_READGROUP[lane_file + config['lane_pair_delim'][0] + '.fastq.gz'][0])
	return(rg_out)

# import CREST hg19 regions
REGIONS_file = config['regions']
if '/home/$USER' in REGIONS_file:
	REGIONS_file = os.environ['HOME'] + REGIONS_file.split('$USER')[-1]
REGIONS = open(REGIONS_file).readlines()
REGIONS = [r.strip() for r in REGIONS]

wildcard_constraints:
	sample='|'.join(list(SAMPLE_LANEFILE.keys())),
	lane = '|'.join(list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in [y for sub in list(SAMPLE_LANEFILE.values()) for y in sub]]))),
	region = '|'.join(REGIONS)

rule all:
	input:
		expand('gvcfs/{sample}.g.vcf.gz', sample=list(SAMPLE_LANEFILE.keys())) if config['GATKgvcf'] == 'TRUE' else 'dummy.txt',
		#'GATK_metrics/multiqc_report' if config['multiqc'] == 'TRUE' else 'dummy.txt',
		'fastqc/multiqc_report' if config['multiqc'] == 'TRUE' else 'dummy.txt',
		expand('picardQC/{sample}.insert_size_metrics.txt', sample=list(SAMPLE_LANEFILE.keys())) if config['picardQC'] == 'TRUE' else 'dummy.txt',
		'CoNVaDING/progress2.done' if config['CoNVaDING'] == 'TRUE' else 'dummy.txt',
		'deepvariant/deepvariantVcf.merge.done.txt' if config['deepvariant'] == 'TRUE' else 'dummy.txt',
		'prioritization/freebayes.merge.done.txt' if config['freebayes_phasing'] == 'TRUE' else 'dummy.txt',
		expand('coverage/{sample}.coverage.xlsx', sample=list(SAMPLE_LANEFILE.keys())) if config['coverage'] == 'TRUE' else 'dummy.txt',
		expand('cram/{sample}.cram', sample=list(SAMPLE_LANEFILE.keys())) if config['cram'] == 'TRUE' else expand('bam/{sample}.bam', sample=list(SAMPLE_LANEFILE.keys())),
		expand('scramble_anno/{sample}.scramble.xlsx', sample=list(SAMPLE_LANEFILE.keys())) if config['SCRAMble'] == 'TRUE' else 'dummy.txt',
		expand('manta/manta.{sample}.annotated.tsv', sample=list(SAMPLE_LANEFILE.keys()))


localrules: dummy
rule dummy:
	input:
		config['metadata_file']
	output:
		temp('dummy.txt')
	shell:
		"""
		touch {output}
		"""

#		expand('INDELseek/{sample}.INDELseek.vcf', sample=list(SAMPLE_LANEFILE.keys())),
#expand('sample_bam/{sample}.bam', sample=list(SAMPLE_LANEFILE.keys())),

# conditinal input:
#rule a:
#    input:
#        name="some/file.txt" if config["condition"] else "other/file.txt"
#    ...

#rule bam_to_fastq:
#	input:
#		PATH + '{lane}.bam'
#	output:
#		forward = temp(PATH  + '{lane}_' + config['fastq_ending_for']),
#		reverse = temp(PATH  + '{lane}_' + config['fastq_ending_rev'])
#	threads: 2
#	shell:
#		"""
#		module load {config[samtools_version]}
#		mkdir -p /scratch/mcgaugheyd/$SLURM_JOB_ID/
#		export REF_CACHE=/scratch/mcgaugheyd/$SLURM_JOB_ID/hts-refcache
#		samtools collate -uOn 128 {wildcards.sample}.bam /scratch/mcgaugheyd/$SLURM_JOB_ID/TMP_{wildcards.sample} | \
#			samtools fastq - -1 {output.forward} -2 {output.reverse}
#		"""

#decided to use sbatch directly
# align with bwa mem

if config['cutadapt'] == 'TRUE':
	rule trim_adatpor:
		input:
			expand('fastq/{{lane}}{pair}.fastq.gz', pair = config['lane_pair_delim'])
		output:
			R1 = temp('trimmed/{lane}_R1_001.fastq.gz'),
			R2 = temp('trimmed/{lane}_R2_001.fastq.gz')
		shell:
			"""
			module load {config[cutadapt_version]}
			cutadapt -a {config[R1_adaptor]} -A {config[R2_adaptor]} --minimum-length 2:2 -o {output.R1} -p {output.R2} {input}
			"""
	rule align:
		input:
			R1 = 'trimmed/{lane}_R1_001.fastq.gz',
			R2 = 'trimmed/{lane}_R2_001.fastq.gz'
		output:
			temp('lane_bam/{lane}.realigned.bam')
		params:
			read_group = rg
		threads: 8
		shell:
			"""
			echo {params.read_group}
			module load {config[bwa-mem2_version]}
			module load {config[samtools_version]}
			bwa-mem2 mem -t {threads} -K 100000000 -M -B 4 -O 6 -E 1 -R {params.read_group} \
				{config[ref_genome]} \
				{input} | \
				samtools view -1 - > \
				{output}
			"""
else:
	rule align:
		input:
			# config['lane_pair_delim'] is the string differentiating
			# the forward from reverse
			# e.g. ['_R1_001', '_R2_001'] if the file names are
			# sample17_R1_001.fastq.gz and sample17_R2_001.fastq.gz
			# for a set of paired end fastq
			# if you don't have a paired fastq set, give as ['']
			expand('fastq/{{lane}}{pair}.fastq.gz', pair = config['lane_pair_delim'])
		output:
			temp('lane_bam/{lane}.realigned.bam')
		params:
			read_group = rg
		threads: 8
		shell:
			"""
			echo {params.read_group}
			module load {config[bwa-mem2_version]}
			module load {config[samtools_version]}
			bwa-mem2 mem -t {threads} -K 100000000 -M -B 4 -O 6 -E 1 -R {params.read_group} \
				{config[ref_genome]} {input} \
				| samtools view -1 - > {output}
			"""
#bwa-mem2 requires 32g mem with 8 threads.
#for WGS use bwa-mem2
			# module load {config[bwa-mem2_version]}
			# module load {config[samtools_version]};
			# bwa-mem2 mem -t {threads} -K 100000000 -M -B 4 -O 6 -E 1 -R {params.read_group} \

rule merge_lane_bam:
	input:
		lambda wildcards: expand('lane_bam/{lane}.realigned.bam', lane = list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in SAMPLE_LANEFILE[wildcards.sample]])))
	output:
		merged_bam = temp('sample_bam/{sample}/{sample}.merged.bam'),
		merged_bai = temp('sample_bam/{sample}/{sample}.merged.bai'),
		bam = temp('sample_bam/{sample}/{sample}.markDup.bam'),
		bai = temp('sample_bam/{sample}/{sample}.markDup.bai'),
		metrics = 'sample_bam/{sample}/{sample}.duplication_metrics.txt'
	threads: 8
	shell:
		"""
		module load {config[picard_version]}
		picard_i=""
		for bam in {input}; do
			picard_i+=" -I $bam"
		done
		java -Xmx16g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			MergeSamFiles \
			-TMP_DIR /lscratch/$SLURM_JOB_ID \
			$picard_i \
			-O {output.merged_bam} \
			--SORT_ORDER coordinate \
			--CREATE_INDEX true
		java -Xmx16g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			MarkDuplicates \
			-TMP_DIR /lscratch/$SLURM_JOB_ID \
			--INPUT {output.merged_bam} \
			--OUTPUT {output.bam} \
			--METRICS_FILE {output.metrics} \
			--CREATE_INDEX true
		"""

#no pipe for picard MergeSamFiles; markduplicates requires indexed bam files.

rule fastqc:
	input:
		bam = 'sample_bam/{sample}/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}/{sample}.markDup.bai'
	output:
		directory('fastqc/{sample}')
	threads: 8
	shell:
		"""
		module load fastqc
		mkdir -p fastqc
		mkdir -p fastqc/{wildcards.sample}
		fastqc -t {threads} -o {output} {input.bam}
		"""
rule multiqc_fastqc:
	input:
		expand('fastqc/{sample}', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		directory('fastqc/multiqc_report')
	shell:
		"""
		module load multiqc
		multiqc -f -o {output} fastqc/
		"""

# rule picard_mark_dups:
# # Mark duplicate reads
# 	input:
# 		bam = 'sample_bam/{sample}/{sample}.b37.bam',
# 		bai = 'sample_bam/{sample}/{sample}.b37.bai'
# 	output:
# 		bam = temp('sample_bam/{sample}.markDup.bam'),
# 		bai = temp('sample_bam/{sample}.markDup.bai'),
# 		metrics = 'GATK_metrics/{sample}.markDup.metrics'
# 	threads: 2
# 	shell:
# 		"""
# 		module load {config[picard_version]}
# 		java -Xmx60g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
# 			MarkDuplicates \
# 			INPUT={input.bam} \
# 			OUTPUT={output.bam} \
# 			METRICS_FILE={output.metrics} \
# 			CREATE_INDEX=true
# 		"""

rule picard_alignmentQC:
#insert size and alignment metrics
	input:
		bam = 'sample_bam/{sample}/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}/{sample}.markDup.bai'
	output:
		insert_size_metrics = 'picardQC/{sample}.insert_size_metrics.txt',
		insert_size_histogram = 'picardQC/{sample}.insert_size_histogram.pdf',
		alignment_metrics = 'picardQC/{sample}.alignment_metrics.txt'
	threads: 2
	shell:
		"""
		module load {config[picard_version]}
		java -Xmx8g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			CollectInsertSizeMetrics \
			INPUT={input.bam} \
			O={output.insert_size_metrics} \
		    H={output.insert_size_histogram} \
		    M=0.5
		java -Xmx8g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			CollectAlignmentSummaryMetrics \
			INPUT={input.bam} \
			R={config[ref_genome]} \
			METRIC_ACCUMULATION_LEVEL=SAMPLE \
			METRIC_ACCUMULATION_LEVEL=READ_GROUP \
			O={output.alignment_metrics}
		"""

#sbatch 8 threads and 32g
localrules: coverage
rule coverage:
	input:
		bam = 'sample_bam/{sample}/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}/{sample}.markDup.bai'
	output:
		thresholds = 'coverage/mosdepth/{sample}.thresholds.bed.gz',
		xlsx = 'coverage/{sample}.coverage.xlsx'
	threads: 8
	shell:
		"""
		module load {config[mosdepth_version]}
		module load {config[R_version]}
		mosdepth -t {threads} --by {config[bed]} --use-median --mapq 1 --thresholds 10,20,30 \
			{wildcards.sample} {input.bam}
		mv {wildcards.sample}.thresholds.bed.gz* coverage/mosdepth/.
		mv {wildcards.sample}.mosdepth* coverage/mosdepth/.
		mv {wildcards.sample}.per-base.bed.gz* coverage/mosdepth/.
		mv {wildcards.sample}.regions.bed.gz* coverage/mosdepth/.
		zcat {output.thresholds} \
			 | sed '1 s/^.*$/chr\tstart\tend\tgene\tcoverageTen\tcoverageTwenty\tcoverageThirty/' \
			 > {output.thresholds}.tsv
		Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/mosdepth_bed_coverage.R \
			{output.thresholds}.tsv {config[OGL_Dx_research_genes]} {output.xlsx}
		"""

# 30% smaller!
rule bam_to_cram:
	input:
		bam = 'sample_bam/{sample}/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}/{sample}.markDup.bai'
	output:
		cram = 'cram/{sample}.cram',
		crai = 'cram/{sample}.crai'
	threads:
		8
	shell:
		"""
		module load {config[samtools_version]}
		samtools sort -O bam -l 0 --threads {threads} -T /lscratch/$SLURM_JOB_ID {input.bam} | \
		samtools view -T {config[ref_genome]} --threads {threads} -C -o {output.cram} -
		samtools index {output.cram} {output.crai}
		"""

localrules: keep_bam
rule keep_bam:
	input:
		bam = 'sample_bam/{sample}/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}/{sample}.markDup.bai'
	output:
		bam = 'bam/{sample}.bam',
		bai = 'bam/{sample}.bai'
	shell:
		"""
		cp -p -l {input.bam} {output.bam}
		cp -p -l {input.bai} {output.bai}
		"""

rule CoNVaDING_1:
	input:
		bam = 'sample_bam/{sample}/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}/{sample}.markDup.bai'
	output:
		temp('CoNVaDING/progress1.{sample}')
	params:
		sex = lambda wildcards: SAMPLE_SEX[wildcards.sample]
	shell:
		"""
		module load {config[samtools_version]}
		module load {config[R_version]}
		case "{params.sex}" in
			"1")
				perl ~/git/CoNVaDING/CoNVaDING.pl -mode StartWithBam \
					-inputDir sample_bam/{wildcards.sample} \
					-outputDir /lscratch/$SLURM_JOB_ID \
					-bed {config[bed]} \
					-useSampleAsControl \
					-controlsDir {config[CoNVaDING_ctr_dir]}_male \
					-rmDup
				mkdir -p CoNVaDING/normalized_coverage_male
				cp /lscratch/$SLURM_JOB_ID/{wildcards.sample}.b37.aligned.only.normalized.coverage.txt \
					CoNVaDING/normalized_coverage_male/{wildcards.sample}.b37.aligned.only.normalized.coverage.txt
				Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/chrRD.R \
					CoNVaDING/normalized_coverage_male/{wildcards.sample}.b37.aligned.only.normalized.coverage.txt \
					CoNVaDING/normalized_coverage_male/{wildcards.sample}.chrRD.pdf \
					{config[chrRD_highcutoff]} \
					{config[chrRD_lowcutoff]} \
					CoNVaDING/normalized_coverage_male/{wildcards.sample}.abnormalChr.tsv \
					1
				touch {output}
				;;
			"2")
				perl ~/git/CoNVaDING/CoNVaDING.pl -mode StartWithBam \
					-inputDir sample_bam/{wildcards.sample} \
					-outputDir /lscratch/$SLURM_JOB_ID \
					-bed {config[bed]} \
					-useSampleAsControl \
					-controlsDir {config[CoNVaDING_ctr_dir]}_female \
					-rmDup
				mkdir -p CoNVaDING/normalized_coverage_female
				cp /lscratch/$SLURM_JOB_ID/{wildcards.sample}.b37.aligned.only.normalized.coverage.txt \
					CoNVaDING/normalized_coverage_female/{wildcards.sample}.b37.aligned.only.normalized.coverage.txt
				Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/chrRD.R \
					CoNVaDING/normalized_coverage_female/{wildcards.sample}.b37.aligned.only.normalized.coverage.txt \
					CoNVaDING/normalized_coverage_female/{wildcards.sample}.chrRD.pdf \
					{config[chrRD_highcutoff]} \
					{config[chrRD_lowcutoff]} \
					CoNVaDING/normalized_coverage_female/{wildcards.sample}.abnormalChr.tsv \
					2
				touch {output}
				;;
			*)
				perl ~/git/CoNVaDING/CoNVaDING.pl -mode StartWithBam \
					-inputDir sample_bam/{wildcards.sample} \
					-outputDir /lscratch/$SLURM_JOB_ID \
					-bed {config[bed]} \
					-useSampleAsControl \
					-controlsDir {config[CoNVaDING_ctr_dir]} \
					-rmDup
				mkdir -p CoNVaDING/normalized_coverage
				cp /lscratch/$SLURM_JOB_ID/{wildcards.sample}.b37.aligned.only.normalized.coverage.txt \
					CoNVaDING/normalized_coverage/{wildcards.sample}.b37.aligned.only.normalized.coverage.txt
				Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/chrRD.R \
					CoNVaDING/normalized_coverage/{wildcards.sample}.b37.aligned.only.normalized.coverage.txt \
					CoNVaDING/normalized_coverage/{wildcards.sample}.chrRD.pdf \
					{config[chrRD_highcutoff]} \
					{config[chrRD_lowcutoff]} \
					CoNVaDING/normalized_coverage/{wildcards.sample}.abnormalChr.tsv \
					0
				touch {output}
				;;
		esac
		"""

				# cp /lscratch/$SLURM_JOB_ID/{wildcards.sample}.b37.aligned.only.normalized.coverage.txt \
				# 	{config[CoNVaDING_ctr_dir]}_male/.

				# cp /lscratch/$SLURM_JOB_ID/{wildcards.sample}.b37.aligned.only.normalized.coverage.txt \
				# 	{config[CoNVaDING_ctr_dir]}_female/.

				# cp /lscratch/$SLURM_JOB_ID/{wildcards.sample}.b37.aligned.only.normalized.coverage.txt \
				# 	{config[CoNVaDING_ctr_dir]}/.

#localrules: CoNVaDING_2
rule CoNVaDING_2:
	input:
		expand('CoNVaDING/progress1.{sample}', sample=list(SAMPLE_LANEFILE.keys())),
	output:
		temp('CoNVaDING/progress2.done')
	shell:
		"""
 		filetest0=$((ls CoNVaDING/normalized_coverage/*.b37.aligned.only.normalized.coverage.txt >> /dev/null 2>&1 && echo TRUE) || echo FALSE)
		if [ $filetest0 == "TRUE" ];
		then
			perl ~/git/CoNVaDING/CoNVaDING.pl -mode StartWithMatchScore \
				-inputDir CoNVaDING/normalized_coverage \
				-outputDir  CoNVaDING/MatchScore \
				-controlsDir {config[CoNVaDING_ctr_dir]}
			perl ~/git/CoNVaDING/CoNVaDING.pl \
  				-mode StartWithBestScore \
  				-inputDir CoNVaDING/MatchScore \
  				-outputDir CoNVaDING/CNV_hiSens \
  				-controlsDir {config[CoNVaDING_ctr_dir]} \
  				-ratioCutOffLow 0.71 \
  				-ratioCutOffHigh 1.35
			# perl ~/git/CoNVaDING/CoNVaDING.pl \
  			# 	-mode GenerateTargetQcList \
  			# 	-inputDir {config[CoNVaDING_ctr_dir]} \
  			# 	-outputDir CoNVaDING/TargetQcList \
  			# 	-controlsDir {config[CoNVaDING_ctr_dir]} \
  			# 	-ratioCutOffLow 0.71 \
  			# 	-ratioCutOffHigh 1.35
			# perl ~/git/CoNVaDING/CoNVaDING.pl \
			# 	-mode CreateFinalList \
			# 	-inputDir CoNVaDING/CNV_hiSens \
  			# 	-targetQcList CoNVaDING/TargetQcList \
  			# 	--outputDir CoNVaDING/finalList
			for i in CoNVaDING/CNV_hiSens/*.shortlist.txt; do awk -F "\t" '{{print FILENAME"\t"$0}}' $i >> CoNVaDING/shortlist.temp; done
			awk -F"\t" 'BEGIN{{OFS="\t"}} {{sub(/CoNVaDING\/CNV_hiSens\//,""); sub(/.b37.aligned.only.best.score.shortlist.txt/,""); print }}' CoNVaDING/shortlist.temp \
				| grep -v -P 'CHR\tSTART' - > CoNVaDING/SHORTlist.txt && \
				echo -e "SAMPLE\tCHR\tSTART\tSTOP\tGENE\tNUMBER_OF_TARGETS\tNUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST\tABBERATION" \
				| cat - CoNVaDING/SHORTlist.txt > CoNVaDING/tmpout && mv CoNVaDING/tmpout CoNVaDING/SHORTlist.txt
		fi
		filetest1=$((ls CoNVaDING/normalized_coverage_male/*.b37.aligned.only.normalized.coverage.txt >> /dev/null 2>&1 && echo TRUE) || echo FALSE)
		if [ $filetest1 == "TRUE" ];
		then
			perl ~/git/CoNVaDING/CoNVaDING.pl -mode StartWithMatchScore \
				-inputDir CoNVaDING/normalized_coverage_male \
				-outputDir  CoNVaDING/MatchScore_male \
				-controlsDir {config[CoNVaDING_ctr_dir]}_male \
				-sexChr
			perl ~/git/CoNVaDING/CoNVaDING.pl \
  				-mode StartWithBestScore \
  				-inputDir CoNVaDING/MatchScore_male \
  				-outputDir CoNVaDING/CNV_hiSens_male \
  				-controlsDir {config[CoNVaDING_ctr_dir]}_male \
  				-ratioCutOffLow 0.71 \
  				-ratioCutOffHigh 1.35 \
				-sexChr
			for i in CoNVaDING/CNV_hiSens_male/*.shortlist.txt; do awk -F "\t" '{{print FILENAME"\t"$0}}' $i >> CoNVaDING/shortlist_male.temp; done
			awk -F"\t" 'BEGIN{{OFS="\t"}} {{sub(/CoNVaDING\/CNV_hiSens_male\//,""); sub(/.b37.aligned.only.best.score.shortlist.txt/,""); print }}' CoNVaDING/shortlist_male.temp \
				| grep -v -P 'CHR\tSTART' - > CoNVaDING/SHORTlist_male.txt && \
				echo -e "SAMPLE\tCHR\tSTART\tSTOP\tGENE\tNUMBER_OF_TARGETS\tNUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST\tABBERATION" \
				| cat - CoNVaDING/SHORTlist_male.txt > CoNVaDING/tmpout_male && mv CoNVaDING/tmpout_male CoNVaDING/SHORTlist_male.txt
		fi
		filetest2=$((ls CoNVaDING/normalized_coverage_female/*.b37.aligned.only.normalized.coverage.txt >> /dev/null 2>&1 && echo TRUE) || echo FALSE)
		if [ $filetest2 == "TRUE" ];
		then
			perl ~/git/CoNVaDING/CoNVaDING.pl -mode StartWithMatchScore \
				-inputDir CoNVaDING/normalized_coverage_female \
				-outputDir  CoNVaDING/MatchScore_female \
				-controlsDir {config[CoNVaDING_ctr_dir]}_female \
				-sexChr
			perl ~/git/CoNVaDING/CoNVaDING.pl \
  				-mode StartWithBestScore \
  				-inputDir CoNVaDING/MatchScore_female \
  				-outputDir CoNVaDING/CNV_hiSens_female \
  				-controlsDir {config[CoNVaDING_ctr_dir]}_female \
  				-ratioCutOffLow 0.71 \
  				-ratioCutOffHigh 1.35 \
				-sexChr
			for i in CoNVaDING/CNV_hiSens_female/*.shortlist.txt; do awk -F "\t" '{{print FILENAME"\t"$0}}' $i >> CoNVaDING/shortlist_female.temp; done
			awk -F"\t" 'BEGIN{{OFS="\t"}} {{sub(/CoNVaDING\/CNV_hiSens_female\//,""); sub(/.b37.aligned.only.best.score.shortlist.txt/,""); print }}' CoNVaDING/shortlist_female.temp \
				| grep -v -P 'CHR\tSTART' - > CoNVaDING/SHORTlist_female.txt && \
				echo -e "SAMPLE\tCHR\tSTART\tSTOP\tGENE\tNUMBER_OF_TARGETS\tNUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST\tABBERATION" \
				| cat - CoNVaDING/SHORTlist_female.txt > CoNVaDING/tmpout_female && mv CoNVaDING/tmpout_female CoNVaDING/SHORTlist_female.txt
		fi
		touch {output}
		"""

### Consider performing CoNVaDING QC and creating final list.
### 				-controlSamples 20 for "StartWithMatchScore" for the male samples because low sample no.
#need 30 samples for step 2 above.

rule deepvariant:
	input:
		bam = 'sample_bam/{sample}/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}/{sample}.markDup.bai'
	output:
		vcf = 'deepvariant/vcf/{sample}.dv.vcf.gz',
		gvcf = 'deepvariant/vcf/{sample}.dv.g.vcf.gz',
		filteredvcf = temp('deepvariant/vcf/{sample}.dv.filtered.vcf.gz'),
		filteretbi = temp('deepvariant/vcf/{sample}.dv.filtered.vcf.gz.tbi'),
		phasedvcf = 'deepvariant/vcf/{sample}.dv.phased.vcf.gz'
	shell:
		"""
		module load {config[deepvariant_version]}
		run_deepvariant --model_type WES --num_shards 4 \
			--ref {config[ref_genome]} \
			--regions {config[padded_bed]} \
			--reads {input.bam} \
			--output_vcf {output.vcf} \
			--output_gvcf {output.gvcf}
		module unload {config[deepvariant_version]}
		module load {config[samtools_version]}
		bcftools norm --multiallelics -any --output-type v {output.vcf} \
			| bcftools norm -d none --output-type v - \
			| bcftools filter --include 'QUAL>0' --output-type z --output {output.filteredvcf}
		sleep 2
		tabix -f -p vcf {output.filteredvcf}
		module load {config[whatshap_version]}
		whatshap phase --reference {config[ref_genome]} --indels {output.filteredvcf} {input.bam} | bgzip -f > {output.phasedvcf}
		tabix -f -p vcf {output.phasedvcf}
		"""
#FILTER="PASS"


localrules: merge_deepvariant_vcf
rule merge_deepvariant_vcf:
	input:
		vcf = expand('deepvariant/vcf/{sample}.dv.phased.vcf.gz', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		'deepvariant/deepvariantVcf.merge.done.txt'
	threads: 8
	shell:
		"""
		module load {config[samtools_version]}
		case "{input.vcf}" in
			*\ *)
				bcftools merge --merge none --missing-to-ref --output-type z --threads {threads} {input.vcf} \
				> deepvariant/{config[analysis_batch_name]}.dv.phased.vcf.gz
				;;
			*)
				cp {input.vcf} deepvariant/{config[analysis_batch_name]}.dv.phased.vcf.gz
				;;
		esac
		sleep 2
		tabix -f -p vcf deepvariant/{config[analysis_batch_name]}.dv.phased.vcf.gz
		touch {output}
		"""

# localrules: merge_deepvariant_gvcf
# rule merge_deepvariant_gvcf:
# 	input:
# 		vcf = expand('deepvariant/vcf/{sample}.dv.g.vcf.gz', sample=list(SAMPLE_LANEFILE.keys()))
# 	#	tbi = expand('freebayes/{sample}.phased.vcf.gz.tbi', sample=list(SAMPLE_LANEFILE.keys()))
# 	output:
# 		'deepvariant/deepvariant.txt'
# 	threads: 8
# 	shell:
# 		"""
# 		module load {config[glnexus_version]}
# 		module load {config[samtools_version]}
# 		glnexus --config DeepVariant --bed {config[padded_bed]} \
# 			--threads {threads} \
# 			{input.vcf} \
# 			| bcftools view - | bgzip -c > deepvariant/{config[analysis_batch_name]}.deepvariant.vcf.gz
# 		tabix -f -p vcf deepvariant/{config[analysis_batch_name]}.deepvariant.vcf.gz
# 		touch {output}
# 		"""
#try phasing gvcf files

#freebayes avoids indel_realign, base-quality_recalibration, better to run on all of the bam files
#freebayes can identify MNP and complex indels.
#need to have unique RG ID fields for each sample. Use the RG made by the wrapper.
#Using CleanSam and FixMateInformation, freebayes was able to call one additional ins variant using the panel NA12878 data.

rule freebayes_phasing:
	input:
		bam = 'sample_bam/{sample}/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}/{sample}.markDup.bai'
	output:
		vcf = 'freebayes/{sample}.vcf.gz',
		filteredvcf = temp('freebayes/{sample}.filtered.vcf.gz'),
		tbi = temp('freebayes/{sample}.filtered.vcf.gz.tbi'),
		phasedvcf = 'freebayes/{sample}.phased.vcf.gz',
		phasedvcf_tbi = 'freebayes/{sample}.phased.vcf.gz.tbi'
	threads: 16
	shell:
		"""
		module load {config[freebayes_version]}
		module load {config[vcflib_version]}
		module load {config[samtools_version]}
		module load {config[vt_version]}
		freebayes-parallel {config[freebayes_region]} {threads} -f {config[ref_genome]} \
			--limit-coverage 1000 {input.bam} --min-alternate-fraction 0.05 \
			--min-mapping-quality 1 --genotype-qualities --strict-vcf --use-mapping-quality \
			| bgzip -f > {output.vcf}
		sleep 2
		tabix -f -p vcf {output.vcf}
		bcftools norm --multiallelics -any --output-type v {output.vcf} \
			| vt decompose_blocksub -p -m -d 2 - \
			| bcftools norm --check-ref s --fasta-ref {config[ref_genome]} --output-type v - \
			| bcftools norm -d none --output-type v - \
			| vcffilter -f "( QUAL > 15 & QA / AO > 15 & SAF > 0 & SAR > 0 & RPR > 0 & RPL > 0 & AO > 2 & DP > 3 ) | ( QUAL > 30 & QA / AO > 25 & ( SAF = 0 | SAR = 0 | RPR = 0 | RPL = 0 ) & AO > 2 & DP > 3 )" \
			| bgzip -f > {output.filteredvcf}
		sleep 2
		tabix -f -p vcf {output.filteredvcf}
		module unload {config[freebayes_version]}
		module unload {config[vcflib_version]}
		module unload {config[vt_version]}
		module load {config[whatshap_version]}
		whatshap phase --reference {config[ref_genome]} --indels {output.filteredvcf} {input.bam} | bgzip -f > {output.phasedvcf}
		tabix -f -p vcf {output.phasedvcf}
		"""
# bcftools norm -d none means remove duplicates if they are identical (keep the first instance)
#vt decompose_blocksub -a separated inframe insertion to fs. thus do not use.
# --gvcf: after gvcf,I tried to pipe it to vcffilter, which removed the reference regions | vcffilter -f "QUAL > 1"
#freebayes -f /data/OGVFB/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta --gvcf --limit-coverage 1000 --min-coverage 4 sample_bam/14_NA12878.markDup.bam > 14_NA12878.test.gvcf
#Tag with "PASS" worked as shown below. It's possible to use --gvcf then tag with PASS?
#vcffilter -t "PASS" -f "( QUAL > 20 & QA / AO > 10 & SAF > 0 & SAR > 0 & RPR > 0 & RPL > 0 )  & AO > 2 & DP > 3 | ( QUAL > 30 & QA / AO > 25 & ( SAF = 0 | SAR = 0 | RPR = 0 | RPL = 0 ) & AO > 2 & DP > 3 )" 14_NA12878.freebayes.vcf.gz | bgzip > 14_NA12878.freebayes.filter7.vcf.gz
#vcffilter -F "filter2" -f ... did not work
##WGS vcffilter
#vcffilter -f "QUAL > 20 & QA / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" {output.vcf} | bgzip > {output.filteredvcf}
# QUAL > 1 removes horrible sites, testing 11/2/19
# QUAL / AO > 10 additional contribution of each obs should be 10 log units (~ Q10 per read), Used QUAL/AO > 5 before 10/30/19;
# QA / AO > 10 (QA: Sum of quality of the alternate observations), starting 11/1/19 - 99% specificity and 98% sensitivity
# SAF > 0 & SAR > 0 reads on both strands, Number of alternate observations on the forward/reverse strand, remove this one for panel/exome data, use AO > 2 instead.
# RPR > 1 & RPL > 1 at least two reads “balanced” to each side of the site


rule merge_freebayes:
	input:
		vcf = expand('freebayes/{sample}.phased.vcf.gz', sample=list(SAMPLE_LANEFILE.keys())),
		tbi = expand('freebayes/{sample}.phased.vcf.gz.tbi', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		'prioritization/freebayes.merge.done.txt'
	threads: 8
	shell:
		"""
		module load {config[samtools_version]}
		case "{input.vcf}" in
			*\ *)
				bcftools merge --merge none --missing-to-ref --output-type z --threads {threads} {input.vcf} \
				> prioritization/{config[analysis_batch_name]}.freebayes.vcf.gz
				;;
			*)
				cp {input.vcf} prioritization/{config[analysis_batch_name]}.freebayes.vcf.gz
				;;
		esac
		sleep 2
		tabix -f -p vcf prioritization/{config[analysis_batch_name]}.freebayes.vcf.gz
		touch {output}
		"""

# rule freebayes:
# 	input:
# 		bam = expand('sample_bam/{sample}.markDup.bam', sample=list(SAMPLE_LANEFILE.keys())),
# 		bai = expand('sample_bam/{sample}.markDup.bai', sample=list(SAMPLE_LANEFILE.keys()))
# 	output:
# 		temp('freebayes.vcf')
# 	threads: 36
# 	shell:
# 		"""
# 		module load {config[freebayes_version]}
# 		module load {config[vcflib_version]}
# 		module load {config[samtools_version]}
# 		freebayes-parallel {config[freebayes_region]} {threads} -f {config[ref_genome]} \
# 			--limit-coverage 1000 {input.bam} | vcffilter -f "QUAL > 1" | bgzip > {config[analysis_batch_name]}.freebayes.vcf.gz
# 		sleep 2
# 		tabix -f -p vcf {config[analysis_batch_name]}.freebayes.vcf.gz
# 		vcffilter -f "QUAL > 20 & QUAL / AO > 5 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1"  {config[analysis_batch_name]}.freebayes.vcf.gz | bgzip > {config[analysis_batch_name]}.freebayes.filtered.vcf.gz
# 		sleep 2
# 		tabix -f -p vcf {config[analysis_batch_name]}.freebayes.filtered.vcf.gz
# 		touch freebayes.vcf
# 		"""
#For filtering: tried QUAL > 20 suggested in GitHub freebayes
#consider change the filter option to according to Eric Garrison's Univ Iowa hardfilter suggestion: vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1"
#For a single file: freebayes -f {config[ref_genome]} {input.bam} | vcffilter -f "QUAL > 20" > {output}
#24 OGLv1 panel on MiSeq 7/17/19, finished in 70 min, when running on 32 threads and 128 gb mem: 500 regions: freebayes-parallel /data/OGVFB/OGL_NGS/bed/freebayes.OGLv1.500.region 32 -f /data/OGVFB/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta 24bam | vcffilter -f "QUAL > 20" > freebayes.vcf
#48 OGLv1 panel on MiSeq 7/17/19, failed after one batch of writing when running on 36 threads and 720gb mem, when on 500 regions, when not setting "--use-best-n-alleles 4"
#48 samples above worked when using --use-best-n-alleles 4 7/21/19 took ~ 3 hours.
#removed --use-best-n-alleles 4 on 7/23/2019 When working with 12 samples, and use "QUAL > 20 & QUAL / AO > 5 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" retained both orf15 variants.

localrules: scramble
rule scramble:
	input:
		bam = 'sample_bam/{sample}/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}/{sample}.markDup.bai'
	output:
		cluster = temp('scramble/{sample}.cluster.txt'),
		mei = 'scramble/{sample}_MEIs.txt',
		deletion = 'scramble/{sample}_PredictedDeletions.txt'
	shell:
		"""
		module load {config[scramble_version]}
		scramble cluster_identifier {input.bam} > {output.cluster}
		scramble Rscript --vanilla /app/cluster_analysis/bin/SCRAMble.R \
			--out-name ${{PWD}}/scramble/{wildcards.sample} \
			--cluster-file ${{PWD}}/{output.cluster} \
			--install-dir /app/cluster_analysis/bin \
			--mei-refs /app/cluster_analysis/resources/MEI_consensus_seqs.fa \
			--ref {config[ref_genome]} \
			--eval-meis \
			--eval-dels \
			--no-vcf
		"""

# change to smaller?       --min-del-len=MIN-DEL-LEN
#                minimum deletion length to call [default 50]
# module load ncbi-toolkit
# makeblastdb -in human_g1k_v37_decoy.fasta -input_type fasta  -dbtype nucl
#--bind /gpfs,/spin1,/data,/lscratch,/scratch,/fdb
localrules: scramble_annotation
rule scramble_annotation:
	input:
		mei = 'scramble/{sample}_MEIs.txt',
		deletion = 'scramble/{sample}_PredictedDeletions.txt'
	output:
		avinput = temp('scramble_anno/{sample}.avinput'),
		annovar = temp('scramble_anno/{sample}.hg19_multianno.txt'),
		annovarR = temp('scramble_anno/{sample}.forR.txt'),
		anno = 'scramble_anno/{sample}.scramble.xlsx'
	shell:
		"""
		module load {config[R_version]}
		module load {config[annovar_version]}
		if [[ $(wc -l {input.mei} | cut -d " " -f 1) == 1 ]]
		then
			touch {output.avinput}
			touch {output.annovar}
			touch {output.annovarR}
			touch {output.anno}
		else
			cut -f 1 {input.mei} | awk -F ":" 'BEGIN{{OFS="\t"}} NR>1 {{print $1,$2,$2,"0","-"}}' > {output.avinput}
			table_annovar.pl {output.avinput} \
				$ANNOVAR_DATA/hg19 \
				-buildver hg19 \
				-remove \
				-out scramble_anno/{wildcards.sample} \
				--protocol refGene \
				-operation  g \
				--argument '-splicing 100 -hgvs' \
				--polish -nastring . \
				--thread 1
			awk -F"\t" 'BEGIN{{OFS="\t"}} NR==1 {{print "Func_refGene","Gene","Intronic","AA"}} NR>1 {{print $6,$7,$8,$10}}' {output.annovar} | paste {input.mei} - > {output.annovarR}
			Rscript /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/scramble_anno.R {output.annovarR} {config[SCRAMBLEdb]} {config[OGL_Dx_research_genes]} {config[HGMDtranscript]} {output.anno} {wildcards.sample}
		fi
		"""
#--intronhgvs 100
localrules: configManta
rule configManta:
	input:
 		bam = 'sample_bam/{sample}/{sample}.markDup.bam',
 		bai = 'sample_bam/{sample}/{sample}.markDup.bai'
	output:
		'manta/{sample}/runWorkflow.py'
	shell:
		"""
		module load {config[manta_version]}
		configManta.py --referenceFasta {config[ref_genome]} --exome --runDir manta/{wildcards.sample} --bam {input.bam}
		"""

localrules: manta
rule manta:
	input:
		'manta/{sample}/runWorkflow.py'
	output:
		'manta/manta.{sample}.annotated.tsv'
	threads: 8
	shell:
		"""
		module load {config[manta_version]}
		manta/{wildcards.sample}/runWorkflow.py -m local -j {threads} -g $((SLURM_MEM_PER_NODE / 1024))
		module load {config[annotsv_version]}
		AnnotSV -SVinputFile manta/{wildcards.sample}/results/variants/diploidSV.vcf.gz \
			-SVinputInfo 1 \
			-outputFile manta/manta.{wildcards.sample}.annotated.tsv
		"""

# localrules: configManta_cohort
# rule configManta_cohort:
# 	input:
#  		bam = expand('sample_bam/{sample}.markDup.bam', sample=list(SAMPLE_LANEFILE.keys())),
#  		bai = expand('sample_bam/{sample}.markDup.bai', sample=list(SAMPLE_LANEFILE.keys()))
# 	output:
# 		'manta_cohort/runWorkflow.py'
# 	shell:
# 		"""
# 		module load {config[manta_version]}
# 		bam_i=""
# 		for bam in {input.bam}; do
# 			bam_i+=" --bam $bam"
# 		done
# 		configManta.py --referenceFasta {config[ref_genome]} --exome --runDir manta_cohort $bam_i
# 		"""
#
# localrules: manta_cohort
# rule manta_cohort:
# 	input:
# 		'manta_cohort/runWorkflow.py'
# 	output:
# 		temp('manta_cohort/manta.cohort.done.txt')
# 	threads: 8
# 	shell:
# 		"""
# 		module load {config[manta_version]}
# 		manta/runWorkflow.py -m local -j {threads} -g $((SLURM_MEM_PER_NODE / 1024))
# 		module load {config[annotsv_version]}
# 		AnnotSV -SVinputFile manta_cohort/results/variants/diploidSV.vcf.gz \
# 			-SVinputInfo 1 \
# 			-outputFile manta_cohort/manta.cohort.{config[analysis_batch_name]}.annotated.tsv
# 		touch {output}
# 		"""
#threads 8; mem=32g as in the list
rule gatk_realigner_target:
# identify regions which need realignment
	input:
		bam = 'sample_bam/{sample}/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}/{sample}.markDup.bai'
	output:
		temp('sample_bam/{sample}.forIndexRealigner.intervals')
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 8g RealignerTargetCreator  \
			-R {config[ref_genome]}  \
			-I {input.bam} \
			--known {config[1000g_indels]} \
			--known {config[mills_gold_indels]} \
			-o {output} \
			--interval_padding 100 \
    		-L {config[bed]}
		"""

rule gatk_indel_realigner:
# realigns indels to improve quality
	input:
		bam = 'sample_bam/{sample}/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}/{sample}.markDup.bai',
		targets = 'sample_bam/{sample}.forIndexRealigner.intervals'
	output:
		bam = temp('sample_bam/{sample}.gatk_realigner.bam'),
		bai = temp('sample_bam/{sample}.gatk_realigner.bai')
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 8g IndelRealigner \
			-R {config[ref_genome]} \
			-I {input.bam} \
			--knownAlleles {config[1000g_indels]} \
			--knownAlleles {config[mills_gold_indels]} \
			-targetIntervals {input.targets} \
			-o {output.bam} \
			--interval_padding 100 \
    		-L {config[bed]}
		"""

rule gatk_base_recalibrator:
# recalculate base quality scores
	input:
		bam = 'sample_bam/{sample}.gatk_realigner.bam',
		bai = 'sample_bam/{sample}.gatk_realigner.bai'
	output:
		'GATK_metrics/{sample}.recal_data.table1'
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 15g BaseRecalibrator  \
			-R {config[ref_genome]} \
			-I {input.bam} \
			--knownSites {config[1000g_indels]} \
			--knownSites {config[mills_gold_indels]} \
			--knownSites {config[dbsnp_var]} \
			-o {output} \
			--interval_padding 100 \
    		-L {config[bed]}
		"""

rule gatk_print_reads:
# print out new bam with recalibrated scoring
	input:
		bam = 'sample_bam/{sample}.gatk_realigner.bam',
		bai = 'sample_bam/{sample}.gatk_realigner.bai',
		bqsr = 'GATK_metrics/{sample}.recal_data.table1'
	output:
		bam = temp('sample_bam/{sample}.recal.bam'),
		bai = temp('sample_bam/{sample}.recal.bai')
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 15g PrintReads \
			-R {config[ref_genome]} \
			-I {input.bam} \
			-BQSR {input.bqsr} \
			-o {output.bam}
		"""

rule gatk_base_recalibrator2:
# recalibrate again
	input:
		bam = 'sample_bam/{sample}.gatk_realigner.bam',
		bai = 'sample_bam/{sample}.gatk_realigner.bai',
		bqsr = 'GATK_metrics/{sample}.recal_data.table1'
	output:
		'GATK_metrics/{sample}.recal_data.table2'
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 15g BaseRecalibrator  \
			-R {config[ref_genome]} \
			-I {input.bam} \
			--knownSites {config[1000g_indels]} \
			--knownSites {config[mills_gold_indels]} \
			--knownSites {config[dbsnp_var]} \
			-BQSR {input.bqsr} \
			-o {output}
			"""

rule gatk_analyze_covariates:
	input:
		one = 'GATK_metrics/{sample}.recal_data.table1',
		two = 'GATK_metrics/{sample}.recal_data.table2'
	output:
		'GATK_metrics/{sample}.BQSRplots.pdf'
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 8g AnalyzeCovariates \
			-R {config[ref_genome]} \
			-before {input.one} \
			-after {input.two} \
			-plots {output}
		"""

rule gatk_haplotype_caller:
# call gvcf
	input:
		bam = 'sample_bam/{sample}.recal.bam',
		bai	= 'sample_bam/{sample}.recal.bai',
		bqsr = 'GATK_metrics/{sample}.recal_data.table1'
	output:
		'gvcfs/{sample}.g.vcf.gz'
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 8g HaplotypeCaller \
			-R {config[ref_genome]} \
			-I {input.bam} \
			--emitRefConfidence GVCF \
			-BQSR {input.bqsr} \
			-o {output} \
			--interval_padding 100 \
    		-L {config[bed]}
		"""


rule multiqc_gatk:
# run multiqc on recalibrator metrics
	input:
		expand('GATK_metrics/{sample}.recal_data.table1', sample=list(SAMPLE_LANEFILE.keys())),
		expand('GATK_metrics/{sample}.recal_data.table2', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		directory('GATK_metrics/multiqc_report')
	shell:
		"""
		module load multiqc
		multiqc -f -o {output} GATK_metrics
		"""
