from os.path import join
import sys
import datetime

def chr_GVCF_to_single_GVCF(wildcards):
	# creates the filenames for the chr level GVCFs to use to concatenate to a single file
	# ensures that input GVCF chrs are provided in order (same as CHRS) below
	sample = str(wildcards)
	sample_by_chr = []
	for chrom in CHRS:
		sample_by_chr.append('gvcfs/chr_split/' + sample + '/' + sample + '__' + str(chrom) + '.g.vcf.gz')
	return(sample_by_chr)

def chr_bam_to_single_bam(wildcards):
	# creates the filenames for the chr level bams to use to concatenate to a single file
	sample = str(wildcards)
	sample_by_chr = []
	for chrom in CHRS:
		sample_by_chr.append('sample_bam/chr_split/' + sample + '/' + sample + '__' + str(chrom) + '.CleanSam.sorted.markDup.gatk_realigner.recalibrated.bam')
	return(sample_by_chr)

# the 'metadata_file' is a csv with three columns
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

if config['analysis_batch_name'] == 'YYYYMMDD':
	currentDT = datetime.datetime.now()
	config['analysis_batch_name'] = currentDT.strftime("%Y%m%d")

def rg(wildcards):
	# returns the read group given in the config['metadata_file']
	lane_file = str(wildcards)
	rg_out = str(LANEFILE_READGROUP[lane_file + config['lane_pair_delim'][0] + '.fastq.gz'][0])
	return(rg_out)

CHRS=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT_contigs"]
MT_CONTIGS="MT GL000207.1 GL000226.1 GL000229.1 GL000231.1 GL000210.1 GL000239.1 GL000235.1 GL000201.1 GL000247.1 GL000245.1 GL000197.1 GL000203.1 GL000246.1 GL000249.1 GL000196.1 GL000248.1 GL000244.1 GL000238.1 GL000202.1 GL000234.1 GL000232.1 GL000206.1 GL000240.1 GL000236.1 GL000241.1 GL000243.1 GL000242.1 GL000230.1 GL000237.1 GL000233.1 GL000204.1 GL000198.1 GL000208.1 GL000191.1 GL000227.1 GL000228.1 GL000214.1 GL000221.1 GL000209.1 GL000218.1 GL000220.1 GL000213.1 GL000211.1 GL000199.1 GL000217.1 GL000216.1 GL000215.1 GL000205.1 GL000219.1 GL000224.1 GL000223.1 GL000195.1 GL000212.1 GL000222.1 GL000200.1 GL000193.1 GL000194.1 GL000225.1 GL000192.1 NC_007605"

wildcard_constraints:
	sample='|'.join(list(SAMPLE_LANEFILE.keys())),
	chr = '|'.join(CHRS),
	lane = '|'.join(list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in [y for sub in list(SAMPLE_LANEFILE.values()) for y in sub]])))


rule all:
	input:
		expand('gvcfs/{sample}.g.vcf.gz', sample=list(SAMPLE_LANEFILE.keys())) if config['GATKgvcf'] == 'TRUE' else 'dummy.txt',
		expand('recal_bam/{sample}.recal.bam', sample=list(SAMPLE_LANEFILE.keys())) if config['recal_bam'] == 'TRUE' else 'dummy.txt',
		expand('sample_cram/{sample}.cram', sample=list(SAMPLE_LANEFILE.keys())) if config['cram'] == 'TRUE' else expand('sample_bam/{sample}.bam', sample=list(SAMPLE_LANEFILE.keys())),
		'GATK_metrics/multiqc_report' if config['multiqc'] == 'TRUE' else 'dummy.txt',
		'fastqc/multiqc_report' if config['multiqc'] == 'TRUE' else 'dummy.txt',
		expand('picardQC/{sample}.insert_size_metrics.txt', sample=list(SAMPLE_LANEFILE.keys())) if config['picardQC'] == 'TRUE' else 'dummy.txt',
		'freebayes_prioritization/freebayes.merged.vcf' if config['freebayes_individual'] == 'TRUE' else 'dummy.txt',
		'freebayes.vcf' if config['freebayes'] == 'TRUE' else 'dummy.txt'
#freebayes_individual: consider gvcf files.

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

# conditinal input:
#rule a:
#    input:
#        name="some/file.txt" if config["condition"] else "other/file.txt"
#    ...

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
	threads: 32
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
		bam = temp('sample_bam/{sample}.b37.bam'),
		bai = temp('sample_bam/{sample}.b37.bai')
	shell:
		"""
		module load {config[samtools_version]}
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
#cp {sample}.bai {sample}.bam.bai

# rule sort:
# 	input:
# 		bam = 'sample_bam/{sample}.bam',
# 		bai = 'sample_bam/{sample}.bam.bai'
# 	output:
# 		temp('sample_bam/{sample}.sorted.bam')
# 	threads: 8
# 	shell:
# 		"""
# 		export REF_CACHE=/scratch/$SLURM_JOB_ID/hts-refcache
# 		module load {config[samtools_version]}
# 		samtools sort {input.bam} -@ {threads} -T /scratch/$SLURM_JOB_ID -o {output}
# 		"""
#
# rule build_index:
# 	input:
# 		'sample_bam/{sample}.sorted.bam'
# 	output:
# 		temp('sample_bam/{sample}.sorted.bam.bai')
# 	threads: 2
# 	shell:
# 		"""
# 		export REF_CACHE=/lscratch/$SLURM_JOB_ID/hts-refcache
# 		module load {config[samtools_version]}
# 		samtools index {input} {output}
# 		"""

rule fastqc:
	input:
		'sample_bam/{sample}.b37.bam'
	output:
		directory('fastqc/{sample}')
	threads: 8
	shell:
		"""
		module load fastqc
		mkdir -p fastqc
		mkdir -p fastqc/{wildcards.sample}
		fastqc -t {threads} -o {output} {input}
		"""

rule picard_mark_dups_allchr:
# Mark duplicate reads
	input:
		bam = 'sample_bam/{sample}.b37.bam',
		bai = 'sample_bam/{sample}.b37.bai'
	output:
		bam = temp('sample_bam/{sample}.markDup.bam'),
		bai = temp('sample_bam/{sample}.markDup.bai'),
		metrics = 'GATK_metrics/{sample}.markDup.metrics'
	threads: 2
	shell:
		"""
		module load {config[picard_version]}
		java -Xmx60g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			MarkDuplicates \
			INPUT={input.bam} \
			OUTPUT={output.bam} \
			METRICS_FILE={output.metrics} \
			CREATE_INDEX=true
		"""

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
		freebayes -f {config[bwa_genome]} \
			--limit-coverage 1000 {input.bam} | vcffilter -f "QUAL > 1" | bgzip > {output.vcf}
		sleep 2
		tabix -f -p vcf {output.vcf}
		vcffilter -f "QUAL > 20 & QUAL / AO > 5 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" {output.vcf} | bgzip > {output.filteredvcf}
		sleep 2
		tabix -f -p vcf {output.filteredvcf}
		"""

rule merge_freebayes:
	input:
		expand('freebayes/{sample}.freebayes.filtered.vcf.gz', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		temp('freebayes_prioritization/freebayes.merged.vcf')
	shell:
		"""
		module load {config[vcftools_version]}
		vcf-merge freebayes/*.freebayes.filtered.vcf.gz | bgzip -c > {config[analysis_batch_name]}.freebayes.filtered.merged.vcf.gz
		sleep 2
		tabix -f -p vcf {config[analysis_batch_name]}.freebayes.filtered.merged.vcf.gz
		touch {output}
		"""
#	module load {config[samtools_version]} - samtools automaticcally loaded when loading vcftools.	

rule create_freebayes_region:
	input:
		bam = expand('sample_bam/{sample}.markDup.bam', sample=list(SAMPLE_LANEFILE.keys())[0]),
		bai = expand('sample_bam/{sample}.markDup.bai', sample=list(SAMPLE_LANEFILE.keys())[0])
	output:
		'freebayes/freebayes.exome.10000.region'
	shell:
		"""
		module load {config[freebayes_version]}
		module load {config[samtools_version]}
		samtools depth {input.bam} | coverage_to_regions.py {config[bwa_genome]}.fai 10000 > {output}
		"""

rule freebayes:
	input:
		bam = expand('sample_bam/{sample}.markDup.bam', sample=list(SAMPLE_LANEFILE.keys())),
		bai = expand('sample_bam/{sample}.markDup.bai', sample=list(SAMPLE_LANEFILE.keys())),
		region = 'freebayes/freebayes.exome.10000.region'
	output:
		temp('freebayes.vcf')
	threads: 36
	shell:
		"""
		module load {config[freebayes_version]}
		module load {config[vcflib_version]}
		module load {config[samtools_version]}
		freebayes-parallel {input.region} {threads} -f {config[bwa_genome]} \
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

# 30% smaller!
rule bam_to_cram:
	input:
		bam = 'sample_bam/{sample}.b37.bam',
		bai = 'sample_bam/{sample}.b37.bai'
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
		bam = 'sample_bam/{sample}.b37.bam',
		bai = 'sample_bam/{sample}.b37.bai'
	output:
		bam = 'sample_bam/{sample}.bam',
		bai = 'sample_bam/{sample}.bai'
	shell:
		"""
		cp -p -l {input.bam} {output.bam}
		cp -p -l {input.bai} {output.bai}
		"""

rule split_bam_by_chr:
	input:
		bam = 'sample_bam/{sample}.b37.bam',
		bai = 'sample_bam/{sample}.b37.bai'
	output:
		temp('sample_bam/chr_split/{sample}/{sample}__{chr}.bam')
	threads: 2
	shell:
		"""
		module load {config[samtools_version]}
		if [[ {wildcards.chr} != "MT_contigs" ]]; then
			samtools view -bh {input.bam} {wildcards.chr} > {output}
		else
			samtools view -bh {input.bam} {MT_CONTIGS}  > {output}
		fi
		"""

rule picard_clean_sam:
# "Soft-clipping beyond-end-of-reference alignments and setting MAPQ to 0 for unmapped reads"
	input:
		'sample_bam/chr_split/{sample}/{sample}__{chr}.bam'
	output:
		temp('sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.bam')
	threads: 2
	shell:
		"""
		module load {config[picard_version]}
		java -Xmx60g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			CleanSam \
			TMP_DIR=/lscratch/$SLURM_JOB_ID \
			INPUT={input} \
			OUTPUT={output}
		"""

rule picard_fix_mate_information:
# "Verify mate-pair information between mates and fix if needed."
# also coord sorts
	input:
		'sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.bam'
	output:
		temp('sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.bam')
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
		'sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.bam'
	output:
		bam = temp('sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.bam'),
		bai = temp('sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.bai'),
		metrics = 'GATK_metrics/{sample}__{chr}.markDup.metrics'
	threads: 2
	shell:
		"""
		module load {config[picard_version]}
		java -Xmx60g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			MarkDuplicates \
			INPUT={input} \
			OUTPUT={output.bam} \
			METRICS_FILE={output.metrics} \
			CREATE_INDEX=true
		"""

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
		bam = 'sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.bam',
		bai = 'sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.bai'
	output:
		temp('sample_bam/chr_split/{sample}/{sample}__{chr}.forIndexRealigner.intervals')
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 8g RealignerTargetCreator  \
			-R {config[ref_genome]}  \
			-I {input.bam} \
			--known {config[1000g_indels]} \
			--known {config[mills_gold_indels]} \
			-o {output}
		"""

rule gatk_indel_realigner:
# realigns indels to improve quality
	input:
		bam = 'sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.bam',
		bai = 'sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.bai',
		targets = 'sample_bam/chr_split/{sample}/{sample}__{chr}.forIndexRealigner.intervals'
	output:
		temp('sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.gatk_realigner.bam')
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
			-o {output}
		"""

rule gatk_base_recalibrator:
# recalculate base quality scores
	input:
		'sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.gatk_realigner.bam'
	output:
		'GATK_metrics/{sample}__{chr}.recal_data.table1'
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 15g BaseRecalibrator  \
			-R {config[ref_genome]} \
			-I {input} \
			--knownSites {config[1000g_indels]} \
			--knownSites {config[mills_gold_indels]} \
			--knownSites {config[dbsnp_var]} \
			-o {output}
		"""

rule gatk_print_reads:
# print out new bam with recalibrated scoring
	input:
		bam = 'sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.gatk_realigner.bam',
		bqsr = 'GATK_metrics/{sample}__{chr}.recal_data.table1'
	output:
		temp('sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.gatk_realigner.recalibrated.bam')
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 15g PrintReads \
			-R {config[ref_genome]} \
			-I {input.bam} \
			-BQSR {input.bqsr} \
			-o {output}
		"""

rule gatk_base_recalibrator2:
# recalibrate again
	input:
		bam = 'sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.gatk_realigner.bam',
		bqsr = 'GATK_metrics/{sample}__{chr}.recal_data.table1'
	output:
		'GATK_metrics/{sample}__{chr}.recal_data.table2'
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
		one = 'GATK_metrics/{sample}__{chr}.recal_data.table1',
		two = 'GATK_metrics/{sample}__{chr}.recal_data.table2'
	output:
		'GATK_metrics/{sample}__{chr}.BQSRplots.pdf'
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
		bam = 'sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.gatk_realigner.recalibrated.bam',
		bqsr = 'GATK_metrics/{sample}__{chr}.recal_data.table1'
	output:
		temp('gvcfs/chr_split/{sample}/{sample}__{chr}.g.vcf.gz')
	threads: 2
	shell:
		"""
		module load {config[gatk_version]}
		GATK -p {threads} -m 8g HaplotypeCaller \
			-R {config[ref_genome]} \
			-I {input.bam} \
			--emitRefConfidence GVCF \
			-BQSR {input.bqsr} \
			-o {output}
		"""

rule picard_merge_bams:
# merge chr split bams into one bam per sample
	input:
		chr_bam_to_single_bam
	output:
		bam = 'recal_bam/{sample}.recal.bam'
	threads: 2
	shell:
		"""
		module load {config[picard_version]}
		cat_inputs_i=""
		for bam in {input}; do
			cat_inputs_i+="I=$bam "; done
		java -Xmx15g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			MergeSamFiles \
			SORT_ORDER=coordinate \
			CREATE_INDEX=true \
			$cat_inputs_i \
			O={output.bam}
		"""
# java -Xmx15g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
# 	BuildBamIndex \
# 	INPUT={output.bam} \
# 	OUTPUT={output.bai}


rule picard_merge_gvcfs:
# merge chr split gvcf back into one gvcf per sample
	input:
		chr_GVCF_to_single_GVCF
	output:
		'gvcfs/{sample}.g.vcf.gz'
	threads: 2
	shell:
		"""
		module load {config[picard_version]}
		cat_inputs_i=""
		for gvcf in {input}; do
			cat_inputs_i+="I=$gvcf "; done
		java -Xmx15g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			MergeVcfs \
			$cat_inputs_i \
			O={output}
		"""

rule multiqc_gatk:
# run multiqc on recalibrator metrics
	input:
		expand('GATK_metrics/{sample}__{chr}.recal_data.table1',sample=list(SAMPLE_LANEFILE.keys()), chr=CHRS),
		expand('GATK_metrics/{sample}__{chr}.recal_data.table2', sample=list(SAMPLE_LANEFILE.keys()), chr=CHRS)
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
