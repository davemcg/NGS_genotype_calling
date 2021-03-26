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


wildcard_constraints:
	sample='|'.join(list(SAMPLE_LANEFILE.keys())),
	lane = '|'.join(list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in [y for sub in list(SAMPLE_LANEFILE.values()) for y in sub]])))

rule all:
	input:
		expand('CREST/hg19.{sample}.markDup.bam.predSV.txt', sample=list(SAMPLE_LANEFILE.keys())) if config['CREST'] == 'TRUE' else 'dummy.txt',
		expand('gvcfs/{sample}.g.vcf.gz', sample=list(SAMPLE_LANEFILE.keys())) if config['GATKgvcf'] == 'TRUE' else 'dummy.txt',
		'GATK_metrics/multiqc_report' if config['multiqc'] == 'TRUE' else 'dummy.txt',
		'fastqc/multiqc_report' if config['multiqc'] == 'TRUE' else 'dummy.txt',
		expand('picardQC/{sample}.insert_size_metrics.txt', sample=list(SAMPLE_LANEFILE.keys())) if config['picardQC'] == 'TRUE' else 'dummy.txt',
		'CoNVaDING/progress2.done' if config['CoNVaDING'] == 'TRUE' else 'dummy.txt',
		'freebayes.vcf' if config['freebayes'] == 'TRUE' else 'dummy.txt',
		expand('freebayes/{sample}.freebayes.filtered.vcf.gz', sample=list(SAMPLE_LANEFILE.keys())) if config['freebayes_individual'] == 'TRUE' else 'dummy.txt',
		expand('sample_cram/{sample}.cram', sample=list(SAMPLE_LANEFILE.keys())) if config['cram'] == 'TRUE' else expand('sample_bam/{sample}.bam', sample=list(SAMPLE_LANEFILE.keys())),

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
#Try samtools rmdup instead in next version? CREST may not read the markDup reads.
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
			REMOVE_DUPLICATES=true \
			CREATE_INDEX=true
		cp {output.bai1} {output.bai2}
		"""

#may be better run with lscratch.
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
#### Create blat hyperlink (can be done in the previous step), use R to remove the common calls found in too many samples, also annotate the region???

localrules: mvCREST
rule mvCREST:
	input:
	 	expand('hg19.{sample}.markDup.bam.predSV.txt', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		expand('CREST/hg19.{sample}.markDup.bam.predSV.txt', sample=list(SAMPLE_LANEFILE.keys()))
	shell:
		"""
		cp -p -l hg19.*.markDup.bam.predSV.txt CREST/.
		rm hg19.*.markDup.bam.predSV.txt
		rm blatServer*.log
		"""
#cp -l: use hard-links instead of copying data of the regular files
#cp -p: preserve attributes;
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
# 30% smaller!
rule bam_to_cram:
	input:
		bam = 'sample_bam/{sample}/{sample}.b37.bam',
		bai = 'sample_bam/{sample}/{sample}.b37.bai'
	output:
		cram = 'sample_cram/{sample}.cram',
		crai = 'sample_cram/{sample}.crai'
	threads:
		8
	shell:
		"""
		module load {config[samtools_version]}
		samtools sort -O bam -l 0 --threads {threads} -T /lscratch/$SLURM_JOB_ID {input.bam} | \
		samtools view -T {config[bwa_genome]} --threads {threads} -C -o {output.cram} -
		samtools index {output.cram} {output.crai}
		"""

localrules: keep_bam
rule keep_bam:
	input:
		bam = 'sample_bam/{sample}/{sample}.b37.bam',
		bai = 'sample_bam/{sample}/{sample}.b37.bai'
	output:
		bam = 'sample_bam/{sample}.bam',
		bai = 'sample_bam/{sample}.bai'
	shell:
		"""
		cp -p -l {input.bam} {output.bam}
		cp -p -l {input.bai} {output.bai}
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

# rule CoNVaDING_1:
# 	input:
# 		bam = 'sample_bam/{sample}/{sample}.b37.bam',
# 		bai = 'sample_bam/{sample}/{sample}.b37.bai'
# 	output:
# 		temp('CoNVaDING/progress1.{sample}')
# 	params:
# 		sex = lambda wildcards: SAMPLE_SEX[wildcards.sample]
# 	run:
# 		shell("bash ~/git/NGS_genotype_calling/NGS_generic_OGL/convading_1.sh {params.sex} {config[samtools_version]} {wildcards.sample} {config[bed]} {config[CoNVaDING_ctr_dir]} {config[chrRD_highcutoff]} {config[chrRD_lowcutoff]} {output}")


rule CoNVaDING_1:
	input:
		bam = 'sample_bam/{sample}/{sample}.b37.bam',
		bai = 'sample_bam/{sample}/{sample}.b37.bai'
	output:
		temp('CoNVaDING/progress1.{sample}')
	params:
		sex = lambda wildcards: SAMPLE_SEX[wildcards.sample]
	shell:
		"""
		module load {config[samtools_version]}
		module load {config[R_version]}
		case {params.sex} in
			1)
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
				cp /lscratch/$SLURM_JOB_ID/{wildcards.sample}.b37.aligned.only.normalized.coverage.txt \
					{config[CoNVaDING_ctr_dir]}
				Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/chrRD.R \
					CoNVaDING/normalized_coverage_male/{wildcards.sample}.b37.aligned.only.normalized.coverage.txt \
					CoNVaDING/normalized_coverage_male/{wildcards.sample}.chrRD.pdf \
					{config[chrRD_highcutoff]} \
					{config[chrRD_lowcutoff]} \
					CoNVaDING/normalized_coverage_male/{wildcards.sample}.abnormalChr.tsv \
					1
				touch {output}
				;;
			2)
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
				cp /lscratch/$SLURM_JOB_ID/{wildcards.sample}.b37.aligned.only.normalized.coverage.txt \
					{config[CoNVaDING_ctr_dir]}
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

localrules: CoNVaDING_2
rule CoNVaDING_2:
	input:
		expand('CoNVaDING/progress1.{sample}', sample=list(SAMPLE_LANEFILE.keys())),
	output:
		temp('CoNVaDING/progress2.done')
	shell:
		"""
 		filetest=$((ls CoNVaDING/normalized_coverage/*.b37.aligned.only.normalized.coverage.txt >> /dev/null 2>&1 && echo TRUE) || echo FALSE)
		if [ $filetest == "TRUE" ];
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
			for i in CoNVaDING/CNV_hiSens/*.shortlist.txt; do awk -F "\t" '{{print FILENAME"\t"$0}}' $i >> CoNVaDING/shortlist.temp; done
			awk -F"\t" 'BEGIN{{OFS="\t"}} {{sub(/CoNVaDING\/CNV_hiSens\//,""); sub(/.b37.aligned.only.best.score.shortlist.txt/,""); print }}' CoNVaDING/shortlist.temp \
				| grep -v -P 'CHR\tSTART' - > CoNVaDING/SHORTlist.txt && \
				echo -e "SAMPLE\tCHR\tSTART\tSTOP\tGENE\tNUMBER_OF_TARGETS\tNUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST\tABBERATION" \
				| cat - CoNVaDING/SHORTlist.txt > CoNVaDING/tmpout && mv CoNVaDING/tmpout CoNVaDING/SHORTlist.txt
		fi
		filetest=$((ls CoNVaDING/normalized_coverage_male/*.b37.aligned.only.normalized.coverage.txt >> /dev/null 2>&1 && echo TRUE) || echo FALSE)
		if [ $filetest == "TRUE" ];
		then
			perl ~/git/CoNVaDING/CoNVaDING.pl -mode StartWithMatchScore \
				-inputDir CoNVaDING/normalized_coverage_male \
				-outputDir  CoNVaDING/MatchScore_male \
				-controlsDir {config[CoNVaDING_ctr_dir]}_male \
				-sexChr \
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
		filetest=$((ls CoNVaDING/normalized_coverage_female/*.b37.aligned.only.normalized.coverage.txt >> /dev/null 2>&1 && echo TRUE) || echo FALSE)
		if [ $filetest == "TRUE" ];
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
rule picard_alignmentQC:
#insert size and alignment metrics
	input:
		bam = 'sample_bam/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}.markDup.bai'
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
			R={config[bwa_genome]} \
			METRIC_ACCUMULATION_LEVEL=SAMPLE \
			METRIC_ACCUMULATION_LEVEL=READ_GROUP \
			O={output.alignment_metrics}
		"""

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

#freebayes avoids indel_realign, base-quality_recalibration, better to run on all of the bam files
#freebayes can identify MNP and complex indels.
#need to have unique RG ID fields for each sample. Use the RG made by the wrapper.

rule freebayes_individual:
	input:
		bam = 'sample_bam/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}.markDup.bai'
	output:
		vcf = 'freebayes/{sample}.freebayes.vcf.gz',
		filteredvcf = 'freebayes/{sample}.freebayes.filtered.vcf.gz'
#	threads: 2
	shell:
		"""
		module load {config[freebayes_version]}
		module load {config[vcflib_version]}
		module load {config[samtools_version]}
		freebayes -f {config[bwa_genome]} --gvcf \
			--limit-coverage 2000 {input.bam} | vcffilter -f "QUAL > 1" | bgzip > {output.vcf}
		sleep 2
		tabix -f -p vcf {output.vcf}
		vcffilter -f "QUAL > 20 & QUAL / AO > 5 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" {output.vcf} | bgzip > {output.filteredvcf}
		sleep 2
		tabix -f -p vcf {output.filteredvcf}
		"""

rule freebayes:
	input:
		bam = expand('sample_bam/{sample}.markDup.bam', sample=list(SAMPLE_LANEFILE.keys())),
		bai = expand('sample_bam/{sample}.markDup.bai', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		temp('freebayes.vcf')
	threads: 36
	shell:
		"""
		module load {config[freebayes_version]}
		module load {config[vcflib_version]}
		module load {config[samtools_version]}
		freebayes-parallel {config[freebayes_region]} {threads} -f {config[bwa_genome]} \
			--limit-coverage 1000 {input.bam} | vcffilter -f "QUAL > 1" | bgzip > {config[analysis_batch_name]}.freebayes.vcf.gz
		sleep 2
		tabix -f -p vcf {config[analysis_batch_name]}.freebayes.vcf.gz
		vcffilter -f "QUAL > 20 & QUAL / AO > 5 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1"  {config[analysis_batch_name]}.freebayes.vcf.gz | bgzip > {config[analysis_batch_name]}.freebayes.filtered.vcf.gz
		sleep 2
		tabix -f -p vcf {config[analysis_batch_name]}.freebayes.filtered.vcf.gz
		touch freebayes.vcf
		"""
#For filtering: tried QUAL > 20 suggested in GitHub freebayes
#consider change the filter option to according to Eric Garrison's Univ Iowa hardfilter suggestion: vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1"
#For a single file: freebayes -f {config[ref_genome]} {input.bam} | vcffilter -f "QUAL > 20" > {output}
#24 OGLv1 panel on MiSeq 7/17/19, finished in 70 min, when running on 32 threads and 128 gb mem: 500 regions: freebayes-parallel /data/OGVFB/OGL_NGS/bed/freebayes.OGLv1.500.region 32 -f /data/OGVFB/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta 24bam | vcffilter -f "QUAL > 20" > freebayes.vcf
#48 OGLv1 panel on MiSeq 7/17/19, failed after one batch of writing when running on 36 threads and 720gb mem, when on 500 regions, when not setting "--use-best-n-alleles 4"
#48 samples above worked when using --use-best-n-alleles 4 7/21/19 took ~ 3 hours.
#removed --use-best-n-alleles 4 on 7/23/2019 When working with 12 samples, and use "QUAL > 20 & QUAL / AO > 5 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" retained both orf15 variants.
