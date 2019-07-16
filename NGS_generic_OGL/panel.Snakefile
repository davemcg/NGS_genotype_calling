from os.path import join
import sys

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
	else:
		old_lane_file = SAMPLE_LANEFILE[sample]
		old_lane_file.append(lane_file)
		SAMPLE_LANEFILE[sample] = old_lane_file
	LANEFILE_READGROUP[lane_file] = [read_group]

def rg(wildcards):
	# returns the read group given in the config['metadata_file']
	lane_file = str(wildcards)
	rg_out = str(LANEFILE_READGROUP[lane_file + config['lane_pair_delim'][0] + '.fastq.gz'][0])
	return(rg_out)


wildcard_constraints:
	sample='|'.join(list(SAMPLE_LANEFILE.keys())),
	lane = '|'.join(list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in [y for sub in list(SAMPLE_LANEFILE.values()) for y in sub]])))

rule all:
	input:
		expand('CREST/hg19.{sample}.markDup.bam.predSV.txt', sample=list(SAMPLE_LANEFILE.keys())),
		expand('gvcfs/{sample}.g.vcf.gz', sample=list(SAMPLE_LANEFILE.keys())),
		'GATK_metrics/multiqc_report',
		'fastqc/multiqc_report',
		'CoNVaDING/SHORTlist.txt'
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
	rule align_hg19:
		input:
			R1 = 'trimmed/{lane}_R1_001.fastq.gz',
			R2 = 'trimmed/{lane}_R2_001.fastq.gz'
		output:
			temp('lane_bam/hg19bam/hg19.{lane}.realigned.bam')
		params:
			read_group = rg
		threads: 8
		shell:
			"""
			echo {params.read_group}
			module load {config[bwa_version]};
			module load {config[samtools_version]};
			bwa mem -t {threads} -K 100000000 -Y -B 4 -O 6 -E 1 -R {params.read_group} \
				/data/OGVFB/OGL_NGS/genomes/hg19/hg19.fa \
				{input} | \
				samtools view -1 - > \
				{output}
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
			module load {config[bwa_version]};
			module load {config[samtools_version]};
			bwa mem -t {threads} -K 100000000 -M -B 4 -O 6 -E 1 -R {params.read_group} \
				{config[bwa_genome]} \
				{input} | \
				samtools view -1 - > \
				{output}
			"""
else:
	rule align_hg19:
		input:
			expand('fastq/{{lane}}{pair}.fastq.gz', pair = config['lane_pair_delim'])
		output:
			temp('lane_bam/hg19bam/hg19.{lane}.realigned.bam')
		params:
			read_group = rg
		threads: 8
		shell:
			"""
			echo {params.read_group}
			module load {config[bwa_version]};
			module load {config[samtools_version]};
			bwa mem -t {threads} -K 100000000 -Y -B 4 -O 6 -E 1 -R {params.read_group} \
				/data/OGVFB/OGL_NGS/genomes/hg19/hg19.fa \
				{input} | \
				samtools view -1 - > \
				{output}
			"""
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
			module load {config[bwa_version]};
			module load {config[samtools_version]};
			bwa mem -t {threads} -K 100000000 -M -B 4 -O 6 -E 1 -R {params.read_group} \
				{config[bwa_genome]} \
				{input} | \
				samtools view -1 - > \
				{output}
			"""


rule merge_lane_bam_hg19:
	input:
		lambda wildcards: expand('lane_bam/hg19bam/hg19.{lane}.realigned.bam', lane = list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in SAMPLE_LANEFILE[wildcards.sample]])))
	output:
		bam = temp('sample_bam/hg19bam/hg19.{sample}.bam'),
		bai = temp('sample_bam/hg19bam/hg19.{sample}.bai')
	shell:
		"""
		module load {config[picard_version]}
		picard_i=""
		for bam in {input}; do
			picard_i+=" I=$bam"
		done
		java -Xmx8g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			MergeSamFiles \
			TMP_DIR=/scratch/$SLURM_JOB_ID \
			$picard_i \
			O={output.bam} \
			SORT_ORDER=coordinate \
			CREATE_INDEX=true
		"""

rule picard_mark_dups_hg19:
# Mark duplicate reads
	input:
		'sample_bam/hg19bam/hg19.{sample}.bam'
	output:
		bam = temp('sample_bam/hg19bam/{sample}/hg19.{sample}.markDup.bam'),
		bai1 = temp('sample_bam/hg19bam/{sample}/hg19.{sample}.markDup.bai'),
		bai2 = temp('sample_bam/hg19bam/{sample}/hg19.{sample}.markDup.bam.bai'),
		metrics = temp('sample_bam/hg19bam/{sample}/GATK_metrics/{sample}.markDup.metrics')
	threads: 2
	shell:
		"""
		module load {config[picard_version]}
		java -Xmx16g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			MarkDuplicates \
			INPUT={input} \
			OUTPUT={output.bam} \
			METRICS_FILE={output.metrics} \
			CREATE_INDEX=true
		cp {output.bai1} {output.bai2}
		"""


rule CREST:
	input:
		bam = 'sample_bam/hg19bam/{sample}/hg19.{sample}.markDup.bam',
		bai = 'sample_bam/hg19bam/{sample}/hg19.{sample}.markDup.bam.bai'
	output:
		cover = temp('hg19.{sample}.markDup.bam.cover'),
		predSV = 'hg19.{sample}.markDup.bam.predSV.txt',
		sclip = temp('hg19.{sample}.markDup.bam.sclip.txt')
	shell:
		"""
		BLAT_PORT=50000
		function safe_start_blat {{
    		local genome_2bit=$1
    		# sleep random amount up to 30s to avoid risk of race condition
    		sleep $((RANDOM % 30))
			# count how many gfServers there are already running
			local other_gfservers=$(ps -eo comm | grep -c gfServer)
			echo "there are $other_gfservers other gfServers running"
			BLAT_PORT=$((BLAT_PORT + other_gfservers))
			# start gfServer
			echo "startig gfServer on port $BLAT_PORT"
			gfServer start localhost $BLAT_PORT $genome_2bit -canStop \
					-log=blatServer_${{BLAT_PORT}}.log &> /dev/null &
			# wait until gfServer is running
			while [[ $(gfServer files localhost $BLAT_PORT 2>&1) =~ "Error in TCP" ]]; do
			 	echo "Waiting for BLAT server to start..."
				sleep 10
			done
			echo "BLAT server is running!"
		}}
		module load CREST
		tmp=$(mktemp -d ./XXXX)
		export TMPDIR=${{tmp}}
		# Start up BLAT server:
		safe_start_blat /data/OGVFB/OGL_NGS/genomes/hg19/hg19.2bit
		trap "gfServer stop localhost $BLAT_PORT; rm -rf ${{tmp}}" EXIT
		extractSClip.pl -i {input.bam} --ref_genome /data/OGVFB/OGL_NGS/genomes/hg19/hg19.fa
		CREST.pl -f {output.cover} \
			-d {input.bam} \
			--ref_genome /data/OGVFB/OGL_NGS/genomes/hg19/hg19.fa \
			-t /data/OGVFB/OGL_NGS/genomes/hg19/hg19.2bit --blatserver localhost \
			--blatport $BLAT_PORT
		"""

localrules: mvCREST
rule mvCREST:
	input:
	 	expand('hg19.{sample}.markDup.bam.predSV.txt', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		expand('CREST/hg19.{sample}.markDup.bam.predSV.txt', sample=list(SAMPLE_LANEFILE.keys()))
	shell:
		"""
		mv hg19.*.markDup.bam.predSV.txt CREST/.
		rm blatServer*.log
		"""
#### Use R to remove the common calls found in too many samples, also annotate the region???

#
# rule align:
# 	input:
# 		# config['lane_pair_delim'] is the string differentiating
# 		# the forward from reverse
# 		# e.g. ['_R1_001', '_R2_001'] if the file names are
# 		# sample17_R1_001.fastq.gz and sample17_R2_001.fastq.gz
# 		# for a set of paired end fastq
# 		# if you don't have a paired fastq set, give as ['']
# 		expand('fastq/{{lane}}{pair}.fastq.gz', pair = config['lane_pair_delim'])
# 	output:
# 		temp('lane_bam/{lane}.realigned.bam')
# 	params:
# 		read_group = rg
# 	threads: 8
# 	shell:
# 		"""
# 		echo {params.read_group}
# 		module load {config[bwa_version]};
# 		module load {config[samtools_version]};
# 		bwa mem -t {threads} -K 100000000 -M -B 4 -O 6 -E 1 -R {params.read_group} \
# 			{config[bwa_genome]} \
# 			{input} | \
# 			samtools view -1 - > \
# 			{output}
# 		"""
# bwa mem -K 100000000 : process input bases in each batch reardless of nThreads (for reproducibility));
# -Y : use soft clipping for supplementary alignments. This is necessary for CREST.
# -M : mark shorter split hits as secondary. David used -M
# flag -M is compatible with lumpy: https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-6-r84
#-B 4 -O 6 -E 1 : these are bwa mem default.

# if multiple sets of fastq/bam provided for a sample, now merge together

rule merge_lane_bam:
	input:
		lambda wildcards: expand('lane_bam/{lane}.realigned.bam', lane = list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in SAMPLE_LANEFILE[wildcards.sample]])))
	output:
		bam = temp('sample_bam/{sample}/{sample}.b37.bam'),
		bai = temp('sample_bam/{sample}/{sample}.b37.bai')
	shell:
		"""
		module load {config[picard_version]}
		picard_i=""
		for bam in {input}; do
			picard_i+=" I=$bam"
		done
		java -Xmx8g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			MergeSamFiles \
			TMP_DIR=/scratch/$SLURM_JOB_ID \
			$picard_i \
			O={output.bam} \
			SORT_ORDER=coordinate \
			CREATE_INDEX=true
		"""

rule fastqc:
	input:
		bam = 'sample_bam/{sample}/{sample}.b37.bam',
		bai = 'sample_bam/{sample}/{sample}.b37.bai'
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

rule CoNVaDING_1:
	input:
		bam = 'sample_bam/{sample}/{sample}.b37.bam',
		bai = 'sample_bam/{sample}/{sample}.b37.bai'
	output:
		'CoNVaDING/normalized_coverage/{sample}.b37.aligned.only.normalized.coverage.txt'
	shell:
		"""
		module load {config[samtools_version]}
		perl ~/git/CoNVaDING/CoNVaDING.pl -mode StartWithBam \
			-inputDir sample_bam/{wildcards.sample} \
			-outputDir /lscratch/$SLURM_JOB_ID \
			-bed {config[bed]} \
			-useSampleAsControl \
			-controlsDir {config[CoNVaDING_ctr_dir]} \
			-rmDup
		cp /lscratch/$SLURM_JOB_ID/*.b37.aligned.only.normalized.coverage.txt CoNVaDING/normalized_coverage/.
		module load {config[R_version]}
		Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/chrRD.R \
			{output} \
			CoNVaDING/normalized_coverage/{wildcards.sample}.chrRD.pdf \
			{config[chrRD_highcutoff]} \
			{config[chrRD_lowcutoff]} \
			CoNVaDING/normalized_coverage/{wildcards.sample}.abnormalChr.tsv
		"""

# rule CoNVaDING_1:
# 	input:
# 		bam = 'sample_bam/{sample}/{sample}.b37.bam',
# 		bai = 'sample_bam/{sample}/{sample}.b37.bai'
# 	output:
# 		'CoNVaDING/normalized_coverage/{sample}.b37.aligned.only.normalized.coverage.txt'
# 	shell:
# 		"""
# 		module load {config[samtools_version]}
# 		perl ~/git/CoNVaDING/CoNVaDING.pl -mode StartWithBam \
# 			-inputDir sample_bam/{sample} \
# 			-outputDir /lscratch/$SLURM_JOB_ID \
# 			-bed {config[bed]} \
# 			-useSampleAsControl \
# 			-controlsDir /lscratch/$SLURM_JOB_ID/controls \
# 			-rmDup
# 		cp -a /lscratch/$SLURM_JOB_ID/{sample}.b37.aligned.only.normalized.coverage.txt CoNVaDING/normalized_coverage/.
# 		cp -a /lscratch/$SLURM_JOB_ID/{sample}.b37.aligned.only.normalized.coverage.txt {config[CoNVaDING_ctr_dir]}/.
# 		module unload {config[samtools_version]}
# 		"""
#
#			-controlsDir {config[CoNVaDING_ctr_dir]} \
### <30 min per sample.
### When sufficient number of male and females present, will have 2 control folders for M and F,
### Will add -sexChr option at that time.

localrules: CoNVaDING_2
rule CoNVaDING_2:
	input:
		expand('CoNVaDING/normalized_coverage/{sample}.b37.aligned.only.normalized.coverage.txt', sample=list(SAMPLE_LANEFILE.keys())),
	output:
		MatchScore = directory('CoNVaDING/MatchScore'),
		hiSens = directory('CoNVaDING/CNV_hiSens'),
		shortlist = 'CoNVaDING/SHORTlist.txt'
	shell:
		"""
		perl ~/git/CoNVaDING/CoNVaDING.pl -mode StartWithMatchScore \
			-inputDir CoNVaDING/normalized_coverage \
			-outputDir  {output.MatchScore} \
			-controlsDir {config[CoNVaDING_ctr_dir]}
		perl ~/git/CoNVaDING/CoNVaDING.pl \
  			-mode StartWithBestScore \
  			-inputDir {output.MatchScore} \
  			-outputDir {output.hiSens} \
  			-controlsDir {config[CoNVaDING_ctr_dir]} \
  			-ratioCutOffLow 0.71 \
  			-ratioCutOffHigh 1.35
		for i in {output.hiSens}/*.shortlist.txt; do awk -F "\t" '{{print FILENAME"\t"$0}}' $i >> CoNVaDING/shortlist.temp; done
		awk -F"\t" 'BEGIN{{OFS="\t"}} {{ sub(/.b37.aligned.only.best.score.shortlist.txt/,""); print }}' CoNVaDING/shortlist.temp \
			| grep -v -P 'CHR\tSTART' - > CoNVaDING/SHORTlist.txt && \
			echo -e "SAMPLE\tCHR\tSTART\tSTOP\tGENE\tNUMBER_OF_TARGETS\tNUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST\tABBERATION" \
			| cat - CoNVaDING/SHORTlist.txt > CoNVaDING/tmpout && mv CoNVaDING/tmpout {output.shortlist}
		"""
#CoNVaDING detect 100% match samples and remove the 100% match samples from control.

#need 30 samples for step 2 above.

#############
##IndelSeek
############

rule picard_clean_sam:
# "Soft-clipping beyond-end-of-reference alignments and setting MAPQ to 0 for unmapped reads"
	input:
		bam = 'sample_bam/{sample}/{sample}.b37.bam',
		bai = 'sample_bam/{sample}/{sample}.b37.bai'
	output:
		temp('sample_bam/{sample}.CleanSam.bam')
	threads: 2
	shell:
		"""
		module load {config[picard_version]}
		java -Xmx60g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			CleanSam \
			TMP_DIR=/lscratch/$SLURM_JOB_ID \
			INPUT={input.bam} \
			OUTPUT={output}
		"""

rule picard_fix_mate_information:
# "Verify mate-pair information between mates and fix if needed."
# also coord sorts
	input:
		'sample_bam/{sample}.CleanSam.bam'
	output:
		temp('sample_bam/{sample}.sorted.bam')
	threads: 2
	shell:
		"""
		module load {config[picard_version]}
		java -Xmx60g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
		FixMateInformation \
			SORT_ORDER=coordinate \
			INPUT={input} \
			OUTPUT={output}
		"""

rule picard_mark_dups:
# Mark duplicate reads
	input:
		'sample_bam/{sample}.sorted.bam'
	output:
		bam = temp('sample_bam/{sample}.markDup.bam'),
		bai = temp('sample_bam/{sample}.markDup.bai'),
		metrics = 'GATK_metrics/{sample}.markDup.metrics'
	threads: 2
	shell:
		"""
		module load {config[picard_version]}
		java -Xmx16g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			MarkDuplicates \
			INPUT={input} \
			OUTPUT={output.bam} \
			METRICS_FILE={output.metrics} \
			CREATE_INDEX=true
		"""
#	bai2 = temp('sample_bam/{sample}.markDup.bam.bai'),
#	cp {output.bai1} {output.bai2}

# rule picard_bam_index:
# # Build bam index
# 	input:
# 		'sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.bam'
# 	output:
# 		temp('sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.bam.bai')
# 	threads: 2
# 	shell:
# 		"""
# 		module load {config[picard_version]}
# 		java -Xmx60g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
# 		BuildBamIndex \
# 			INPUT={input} \
# 			OUTPUT={output}
# 		"""

rule gatk_realigner_target:
# identify regions which need realignment
	input:
		bam = 'sample_bam/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}.markDup.bai'
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
		bam = 'sample_bam/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}.markDup.bai',
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
		bam = 'sample_bam/{sample}.bam',
		bai1 = 'sample_bam/{sample}.bai',
		bai2 = 'sample_bam/{sample}.bam.bai'
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 15g PrintReads \
			-R {config[ref_genome]} \
			-I {input.bam} \
			-BQSR {input.bqsr} \
			-o {output.bam}
		cp {output.bai1} {output.bai2}
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
		bam = 'sample_bam/{sample}.bam',
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
		expand('GATK_metrics/{sample}.recal_data.table1',sample=list(SAMPLE_LANEFILE.keys())),
		expand('GATK_metrics/{sample}.recal_data.table2', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		directory('GATK_metrics/multiqc_report')
	shell:
		"""
		module load multiqc
		multiqc -f -o {output} GATK_metrics
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
