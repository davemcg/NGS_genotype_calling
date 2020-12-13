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

def chr_fbVCF_to_single_fbVCF_filtered(wildcards):
	# creates the filenames for the chr level GVCFs to use to concatenate to a single file
	# ensures that input GVCF chrs are provided in order (same as CHRS) below
	sample = str(wildcards)
	sample_by_chr = []
	for chrom in CHRS:
		sample_by_chr.append('freebayes/chr_split/' + sample + '/' + sample + '__' + str(chrom) + '.filtered.vcf.gz')
	return(sample_by_chr)

def chr_fbVCF_to_single_fbVCF_phased(wildcards):
	# creates the filenames for the chr level GVCFs to use to concatenate to a single file
	# ensures that input GVCF chrs are provided in order (same as CHRS) below
	sample = str(wildcards)
	sample_by_chr = []
	for chrom in CHRS:
		sample_by_chr.append('freebayes/chr_split/' + sample + '/' + sample + '__' + str(chrom) + '.phased.vcf.gz')
	return(sample_by_chr)

def chr_fbVCF_to_single_fbVCF(wildcards):
	# creates the filenames for the chr level GVCFs to use to concatenate to a single file
	# ensures that input GVCF chrs are provided in order (same as CHRS) below
	sample = str(wildcards)
	sample_by_chr = []
	for chrom in CHRS:
		sample_by_chr.append('freebayes/chr_split/' + sample + '/' + sample + '__' + str(chrom) + '.vcf.gz')
	return(sample_by_chr)

def chr_bam_to_single_bam(wildcards):
	# creates the filenames for the chr level bams to use to concatenate to a single file
	sample = str(wildcards)
	sample_by_chr = []
	for chrom in CHRS:
		sample_by_chr.append('sample_bam/chr_split/' + sample + '/' + sample + '__' + str(chrom) + '.CleanSam.sorted.markDup.bam')
	return(sample_by_chr)

def chr_scramble_to_single_scramble(wildcards):
	# creates the filenames for the chr level GVCFs to use to concatenate to a single file
	# ensures that input GVCF chrs are provided in order (same as CHRS) below
	sample = str(wildcards)
	sample_by_chr = []
	for chrom in CHRS:
		sample_by_chr.append('scramble/chr_split/' + sample + '/' + sample + '__' + str(chrom) + '.scramble.tsv')
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

if config['inputFileType'] == 'fastq':
	def rg(wildcards):
		# returns the read group given in the config['metadata_file']
		lane_file = str(wildcards)
		rg_out = str(LANEFILE_READGROUP[lane_file + config['lane_pair_delim'][0] + '.fastq.gz'][0])
		return(rg_out)
else:
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
		expand('cram/{sample}.cram', sample=list(SAMPLE_LANEFILE.keys())) if config['cram'] == 'TRUE' else expand('bam/{sample}.bam', sample=list(SAMPLE_LANEFILE.keys())),
		# 'GATK_metrics/multiqc_report' if config['multiqc'] == 'TRUE' else 'dummy.txt',
		'fastqc/multiqc_report' if config['multiqc'] == 'TRUE' else 'dummy.txt',
		# expand('picardQC/{sample}.insert_size_metrics.txt', sample=list(SAMPLE_LANEFILE.keys())) if config['picardQC'] == 'TRUE' else 'dummy.txt',
		'prioritization/freebayes.merge.done.txt' if config['freebayes_phasing'] == 'TRUE' else 'dummy.txt',
		expand('deepvariant/{sample}/{sample}.vcf.gz', sample=list(SAMPLE_LANEFILE.keys())) if config['deepvariant'] == 'TRUE' else 'dummy.txt',
		expand('manta/manta.{sample}.annotated.tsv', sample=list(SAMPLE_LANEFILE.keys())),
		expand('scramble_anno/{sample}.scramble.xlsx', sample=list(SAMPLE_LANEFILE.keys())) if config['SCRAMble'] == 'TRUE' else 'dummy.txt'
#		expand('freebayes/{sample}.freebayes.filtered.vcf.gz', sample=list(SAMPLE_LANEFILE.keys())) if config['freebayes_individual'] == 'TRUE' else 'dummy.txt',
#		'freebayes.vcf' if config['freebayes'] == 'TRUE' else 'dummy.txt'
#freebayes_individual can consider gvcf files.

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

##may reduce to 36 threads: https://blog.dnanexus.com/2020-03-10-bwa-mem2-review/
##tested bwa-mem2, sambamba or with biobambam2. Sambamba shown here is 1.5-2x faster.
if config['inputFileType'] == 'single_lane_fastq':
	rule align_markdup:
		input:
			expand('fastq/{{lane}}{pair}.fastq.gz', pair = config['lane_pair_delim'])
		output:
			bam = temp('lane_bam/{lane}.bam'),
			bai = temp('lane_bam/{lane}.bai'),
			metrics = temp('lane_bam/{lane}.duplication_metrics.txt')
		params:
			read_group = rg
		threads: 56
		shell:
			"""
			export TMPDIR=/lscratch/$SLURM_JOB_ID
			echo {params.read_group}
			module load {config[bwa-mem2_version]} {config[samblaster_version]} {config[sambamba_version]}
			bwa-mem2 mem -t {threads} -K 100000000 -M -Y -B 4 -O 6 -E 1 -R {params.read_group} \
				{config[bwa-mem2_ref]} {input} \
			 	| samblaster -M --addMateTags --quiet \
				| sambamba sort -u --tmpdir=/lscratch/$SLURM_JOB_ID -t {threads} -o {output.bam} \
					<(sambamba view -S -f bam -l 0 -t $SLURM_CPUS_PER_TASK /dev/stdin)
			"""
	localrules: cp_lane_bam
	rule cp_lane_bam:
		input:
			bam = lambda wildcards: expand('lane_bam/{lane}.bam', lane = list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in SAMPLE_LANEFILE[wildcards.sample]]))),
			bai = lambda wildcards: expand('lane_bam/{lane}.bai', lane = list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in SAMPLE_LANEFILE[wildcards.sample]]))),
			metrics = lambda wildcards: expand('lane_bam/{lane}.duplication_metrics.txt', lane = list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in SAMPLE_LANEFILE[wildcards.sample]])))
		output:
			bam = temp('sample_bam/{sample}.markDup.bam'),
			bai = temp('sample_bam/{sample}.markDup.bai'),
			metrics = 'sample_bam/{sample}.duplication_metrics.txt'
		shell:
			"""
			cp -p -l {input.bam} {output.bam}
			cp -p -l {input.bai} {output.bai}
			cp -p -l {input.metrics} {output.metrics}
			"""
elif config['inputFileType'] == 'bam':
	rule realign:
		input:
			lambda wildcards: join('old_bam/', str(SAMPLE_LANEFILE[wildcards.sample][0]))
		output:
			bam = temp('sample_bam/{sample}.markDup.bam'),
			bai = temp('sample_bam/{sample}.markDup.bai'),
			metrics = 'sample_bam/{sample}.duplication_metrics'
		threads: 36
		shell:
			"""
			export TMPDIR=/lscratch/$SLURM_JOB_ID
			module load {config[bazam_version]} {config[bwa-mem2_version]} {config[biobambam2_version]}
 			case "{input}" in
				*bam)
					java -Xmx12g -jar $BAZAMPATH/bazam.jar -bam {input} \
					| bwa-mem2 mem -t {threads} -K 100000000 -M -B 4 -O 6 -E 1 -p {config[bwa-mem2_ref]} - \
					| bamsormadup SO=coordinate threads={threads} level=6 inputformat=sam \
					tmpfile=/lscratch/$SLURM_JOB_ID/bamsormadup \
					indexfilename={output.bai} M={output.metrics} > {output.bam}
				*cram)
					java -Xmx12g -Dsamjdk.reference_fasta={config[old_cram_ref]} -jar $BAZAMPATH/bazam.jar -bam {input} \
					| bwa-mem2 mem -t {threads} -K 100000000 -M -B 4 -O 6 -E 1 -p {config[bwa-mem2_ref]} - \
					| bamsormadup SO=coordinate threads={threads} level=6 inputformat=sam \
					tmpfile=/lscratch/$SLURM_JOB_ID/bamsormadup \
					indexfilename={output.bai} M={output.metrics} > {output.bam}
			esac
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
			temp('lane_bam/{lane}.bam')
		params:
			read_group = rg
		threads: 36
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
	rule merge_lane_bam:
		input:
			lambda wildcards: expand('lane_bam/{lane}.bam', lane = list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in SAMPLE_LANEFILE[wildcards.sample]])))
		output:
			merged_bam = temp('sample_bam/{sample}.merged.bam'),
			merged_bai = temp('sample_bam/{sample}.merged.bai'),
			bam = temp('sample_bam/{sample}.markDup.bam'),
			bai = temp('sample_bam/{sample}.markDup.bai'),
			metrics = 'sample_bam/{sample}.duplication_metrics.txt'
		threads: 16
		shell:
			"""
			module load {config[picard_version]}
			picard_i=""
			for bam in {input}; do
				picard_i+=" -I $bam"
			done
			java -Xmx32g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
				MergeSamFiles \
				-TMP_DIR /lscratch/$SLURM_JOB_ID \
				$picard_i \
				-O {output.merged_bam} \
				--SORT_ORDER coordinate \
				--CREATE_INDEX true
			java -Xmx32g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
				MarkDuplicates \
				-TMP_DIR /lscratch/$SLURM_JOB_ID \
				--INPUT {output.merged_bam} \
				--OUTPUT {output.bam} \
				--METRICS_FILE {output.metrics} \
				--CREATE_INDEX true
			"""
#samblaster --addMateTags for SV calling later? samblaster also MarkDuplicates
#The following bwa-mem2 and biobambam2 worked, but one WGS sample did not get processed well (28 threads < 25g mememory used, <9 hours);
# if using 56 or more threads, it likely will be faster.
			# module load {config[bwa-mem2_version]}
			# module load {config[biobambam2_version]}
			# bwa-mem2 mem -t {threads} -K 100000000 -M -B 4 -O 6 -E 1 -R {params.read_group} \
			# 	{config[ref_genome]} {input} \
			# 	| bamsormadup SO=coordinate threads={threads} level=6 inputformat=sam \
			# 	tmpfile=/lscratch/$SLURM_JOB_ID/bamsormadup \
			# 	indexfilename={output.bai} M={output.metrics} > {output.bam}

# export TMPDIR=/lscratch/$SLURM_JOB_ID necessary?

# bwa mem -K 100000000 : process input bases in each batch reardless of nThreads (for reproducibility));
# -Y : use soft clipping for supplementary alignments. This is necessary for CREST.
# -M : mark shorter split hits as secondary. David used -M
# flag -M is compatible with lumpy: https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-6-r84
#-B 4 -O 6 -E 1 : these are bwa mem default.

rule fastqc:
	input:
		bam = 'sample_bam/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}.markDup.bai'
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
			R={config[ref_genome]} \
			METRIC_ACCUMULATION_LEVEL=SAMPLE \
			METRIC_ACCUMULATION_LEVEL=READ_GROUP \
			O={output.alignment_metrics}
		"""

# 30% smaller!
rule bam_to_cram:
	input:
		bam = 'sample_bam/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}.markDup.bai'
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
		bam = 'sample_bam/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}.markDup.bai'
	output:
		bam = 'bam/{sample}.bam',
		bai = 'bam/{sample}.bai'
	shell:
		"""
		cp -p -l {input.bam} {output.bam}
		cp -p -l {input.bai} {output.bai}
		"""

localrules: deepvariantS1
rule deepvariantS1:
	input:
		bam = 'sample_bam/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}.markDup.bai'
	output:
		examples = directory('deepvariant/{sample}/examples'),
		gvcf_tfrecords = directory('deepvariant/{sample}/gvcf_tfrecords')
	threads: 56
	shell:
		"""
		module load {config[deepvariant_version]} parallel
		mkdir -p deepvariant/{wildcards.sample}/examples deepvariant/{wildcards.sample}/gvcf_tfrecords
		BASE=$PWD
		PROJECT="${{BASE}}/deepvariant/{wildcards.sample}"
		N_SHARDS=$SLURM_CPUS_PER_TASK
		WORK_DIR="/lscratch/${{SLURM_JOB_ID}}"
		DV_OUTPUT_DIR="${{WORK_DIR}}/output"
		DV_INPUT_DIR="${{WORK_DIR}}/input"
		TMPDIR="${{WORK_DIR}}/temp"
		LOG="${{WORK_DIR}}/logs"
		EXAMPLES="${{DV_OUTPUT_DIR}}/{wildcards.sample}.examples.tfrecord@${{N_SHARDS}}.gz"
		GVCF_TFRECORDS="${{DV_OUTPUT_DIR}}/{wildcards.sample}.gvcf.tfrecord@${{N_SHARDS}}.gz"
		rm -rf $WORK_DIR/input $WORK_DIR/output $WORK_DIR/logs $WORK_DIR/temp
		mkdir -p $WORK_DIR/input $WORK_DIR/output $WORK_DIR/logs $WORK_DIR/temp
		cp {input} {config[ref_genome]} {config[ref_genome]}.fai $DV_INPUT_DIR || exit 1
		( time seq 0 $((N_SHARDS-1)) | parallel -j {threads} --eta --halt 2 --line-buffer \
			make_examples --mode calling \
			--ref  $DV_INPUT_DIR/$(basename {config[ref_genome]}) \
			--reads $DV_INPUT_DIR/$(basename {input.bam}) \
			--examples "${{EXAMPLES}}" \
			--gvcf "${{GVCF_TFRECORDS}}" \
			--task {{}} \
		) 2>&1 | tee "${{LOG}}/{wildcards.sample}.make_examples.log"
		echo "Done." || exit 1
		cp ${{DV_OUTPUT_DIR}}/{wildcards.sample}.examples.* ${{PROJECT}}/examples
		cp ${{DV_OUTPUT_DIR}}/{wildcards.sample}.gvcf.* ${{PROJECT}}/gvcf_tfrecords || exit 1
		"""
	# output:
	# 	tfrecord = temp('deepvariant/{sample}.tfrecord.gz'),
	# 	gvcf_tfrecord = temp('deepvariant/{sample}.gvcf.tfrecord.gz'),
	# 	call_variants = temp('deepvariant/{sample}.call_variants.tfrecord.gz'),
	# 	vcf = 'deepvariant/vcf/{sample}.dv.vcf.gz',
	# 	gvcf = 'deepvariant/gvcf/{sample}.dv.g.vcf.gz'
	# threads: 56
	# shell:
	# 	"""
	# 	module load {config[deepvariant_version]} parallel
	# 	BASE=$PWD
	# 	N_SHARDS="${SLURM_CPUS_PER_TASK}"
	# 	time seq 0 $((N_SHARDS-1)) | parallel -j {threads} --eta --halt 2 --line-buffer \
	# 		make_examples --mode calling \
	# 		--ref {config[ref_genome]} \
	# 		--reads {input.bam} \
	# 		--examples {output.tfrecord} \
	# 		--gvcf {output.gvcf_tfrecord} \
	# 		--task {}
	# 	call_variants \
	# 		--outfile {output.call_variants} \
	# 		--examples {output.tfrecord} \
	# 		--checkpoint /opt/models/wes/model.ckpt
	# 	postprocess_variants \
	# 		--ref {config[ref_genome]} \
	# 		--infile {output.call_variants} \
	# 		--outfile {output.vcf} \
	# 		--nonvariant_site_tfrecord_path {output.gvcf_tfrecord} \
	# 		--gvcf_outfile {output.gvcf}
	# 	"""
	#
localrules: deepvariantS2
rule deepvariantS2:
	input:
		'deepvariant/{sample}/examples'
	output:
		'deepvariant/{sample}/call_variants/{sample}.call_variants_output.tfrecord.gz'
	threads: 24
	shell:
		"""
		module load {config[deepvariant_version]}
		BASE=$PWD
		PROJECT="${{BASE}}/deepvariant/{wildcards.sample}"
		N_SHARDS="56"
		WORK_DIR="/lscratch/${{SLURM_JOB_ID}}"
		DV_OUTPUT_DIR="${{WORK_DIR}}/output"
		DV_INPUT_DIR="${{WORK_DIR}}/input"
		TMPDIR="${{WORK_DIR}}/temp"
		LOG="${{WORK_DIR}}/logs"
		EXAMPLES="${{DV_OUTPUT_DIR}}/{wildcards.sample}.examples.tfrecord@${{N_SHARDS}}.gz"
		GVCF_TFRECORDS="${{DV_OUTPUT_DIR}}/{wildcards.sample}.gvcf.tfrecord@${{N_SHARDS}}.gz"
		CALL_VARIANTS_OUTPUT="${{DV_OUTPUT_DIR}}/{wildcards.sample}.call_variants_output.tfrecord.gz"
		rm -rf $WORK_DIR/input $WORK_DIR/output $WORK_DIR/logs $WORK_DIR/temp
		mkdir -p $WORK_DIR/input $WORK_DIR/output $WORK_DIR/logs $WORK_DIR/temp
		cp ${{PROJECT}}/examples/* ${{DV_OUTPUT_DIR}} || exit 1
		MODEL="/opt/models/wgs/model.ckpt"
		( time \
    		call_variants \
    		--outfile "${{CALL_VARIANTS_OUTPUT}}" \
    		--examples "${{EXAMPLES}}" \
    		--checkpoint "${{MODEL}}" \
		) 2>&1 | tee "${{LOG}}/{wildcards.sample}.call_variants.log"
		echo "Done." || exit 1
		cp ${{CALL_VARIANTS_OUTPUT}} ${{PROJECT}}/call_variants || exit 1
		"""
localrules: deepvariantS3
rule deepvariantS3:
	input:
		'deepvariant/{sample}/call_variants/{sample}.call_variants_output.tfrecord.gz'
	output:
		vcf = 'deepvariant/{sample}/{sample}.vcf.gz',
		gvcf = 'deepvariant/{sample}/{sample}.g.vcf.gz'
	shell:
		"""
		module load {config[deepvariant_version]} parallel
		BASE=$PWD
		PROJECT="${{BASE}}/deepvariant/{wildcards.sample}"
		N_SHARDS="56"
		WORK_DIR="/lscratch/${{SLURM_JOB_ID}}"
		DV_OUTPUT_DIR="${{WORK_DIR}}/output"
		DV_INPUT_DIR="${{WORK_DIR}}/input"
		TMPDIR="${{WORK_DIR}}/temp"
		LOG="${{WORK_DIR}}/logs"
		EXAMPLES="${{DV_OUTPUT_DIR}}/{wildcards.sample}.examples.tfrecord@${{N_SHARDS}}.gz"
		GVCF_TFRECORDS="${{DV_OUTPUT_DIR}}/{wildcards.sample}.gvcf.tfrecord@${{N_SHARDS}}.gz"
		CALL_VARIANTS_OUTPUT="${{DV_OUTPUT_DIR}}/{wildcards.sample}.call_variants_output.tfrecord.gz"
		CV_FILE="{wildcards.sample}.call_variants_output.tfrecord.gz"
		OUTPUT_VCF="${{DV_OUTPUT_DIR}}/{wildcards.sample}.vcf.gz"
		OUTPUT_GVCF="${{DV_OUTPUT_DIR}}/{wildcards.sample}.g.vcf.gz"
		rm -rf $WORK_DIR/input $WORK_DIR/output $WORK_DIR/logs $WORK_DIR/temp
		mkdir -p $WORK_DIR/input $WORK_DIR/output $WORK_DIR/logs $WORK_DIR/temp
		cp ${{PROJECT}}/gvcf_tfrecords/* {config[ref_genome]} {config[ref_genome]}.fai ${{PROJECT}}/call_variants/${{CV_FILE}} ${{DV_OUTPUT_DIR}} || exit 1
		( time \
    		postprocess_variants \
      		--ref $DV_OUTPUT_DIR/$(basename {config[ref_genome]}) \
      		--infile "${{CALL_VARIANTS_OUTPUT}}" \
      		--outfile "${{OUTPUT_VCF}}" \
      		--nonvariant_site_tfrecord_path "${{GVCF_TFRECORDS}}" \
      		--gvcf_outfile "${{OUTPUT_GVCF}}" \
		) 2>&1 | tee "${{LOG}}/{wildcards.sample}.postprocess_variants.log" || exit 1
		cp ${{OUTPUT_VCF}}* ${{OUTPUT_GVCF}}* ${{DV_OUTPUT_DIR}}/*.html ${{PROJECT}} || exit 1
		"""
#deepvariant separates MNPs, thus no need to use decompose_blocksub.
#if too slow, then seperate by chr as above.
#localrules: scramble
rule scramble:
	input:
		bam = 'sample_bam/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}.markDup.bai'
	output:
		cluster = temp('scramble/{sample}.cluster.txt'),
		mei = 'scramble/{sample}.mei.txt',
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
#--bind /gpfs,/spin1,/data,/lscratch,/scratch,/fdb
#localrules: scramble_annotation
rule scramble_annotation:
	input:
		mei = 'scramble/{sample}.mei.txt'
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

#Try this strategy,if not working well, consider split by chr as in the GATK.
localrules: configManta
rule configManta:
	input:
 		bam = 'sample_bam/{sample}.markDup.bam',
 		bai = 'sample_bam/{sample}.markDup.bai'
	output:
		'manta/{sample}/runWorkflow.py'
	shell:
		"""
		module load {config[manta_version]}
		configManta.py --referenceFasta {config[ref_genome]} --exome --runDir manta/{wildcards.sample} --bam {input.bam}
		"""

rule manta:
	input:
		'manta/{sample}/runWorkflow.py'
	output:
		'manta/manta.{sample}.annotated.tsv'
	threads: 32
	shell:
		"""
		module load {config[manta_version]}
		manta/{wildcards.sample}/runWorkflow.py -m local -j {threads} -g $((SLURM_MEM_PER_NODE / 1024))
		module load {config[annotsv_version]}
		AnnotSV -SVinputFile manta/{wildcards.sample}/results/variants/diploidSV.vcf.gz \
			-SVinputInfo 1 \
			-outputFile manta/manta.{wildcards.sample}.annotated.tsv
		"""

#SV, use 16 threads and 64 g
rule sve:
	input:
		bam = 'sample_bam/{sample}.bam',
		bai = 'sample_bam/{sample}.bai',
		bai2 = 'sample_bam/{sample}.bam.bai'
	output:
		vcf = 'sve/{sample}/{sample}_S4.vcf'
	shell:
		"""
		module load sve/0.1.0
		for method in breakdancer breakseq cnvnator delly lumpy cnmops; do
		 	sve call -r {config[ref_genome]} -o sve/{wildcards.sample} \
			-t 16 -M 4 -g hg19 -a $method {input.bam}; done
		"""


rule split_bam_by_chr:
	input:
		bam = 'sample_bam/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}.markDup.bai'
	output:
		bam = temp('sample_bam/chr_split/{sample}/{sample}__{chr}.bam'),
		bai = temp('sample_bam/chr_split/{sample}/{sample}__{chr}.bai'),
		metrics = temp('sample_bam/chr_split/{sample}/{sample}__{chr}.markdup.metrics.txt')
	threads: 2
	shell:
		"""
		module load {config[samtools_version]} {config[biobambam2_version]}
		if [[ {wildcards.chr} != "MT_contigs" ]]; then
			samtools view -h --threads {threads} --output-fmt SAM {input.bam} {wildcards.chr} \
			| bamsormadup SO=coordinate threads={threads} level=6 inputformat=sam \
			tmpfile=/lscratch/$SLURM_JOB_ID/bamsormadup \
			indexfilename={output.bai} M={output.metrics} > {output.bam}
		else
			samtools view -h --threads {threads} --output-fmt SAM {input.bam} {MT_CONTIGS} \
			| bamsormadup SO=coordinate threads={threads} level=6 inputformat=sam \
			tmpfile=/lscratch/$SLURM_JOB_ID/bamsormadup \
			indexfilename={output.bai} M={output.metrics} > {output.bam}
		fi
		"""
rule freebayes_phasing:
	input:
		bam = 'sample_bam/chr_split/{sample}/{sample}__{chr}.bam',
		bai = 'sample_bam/chr_split/{sample}/{sample}__{chr}.bai'
	output:
		vcf = temp('freebayes/chr_split/{sample}/{sample}__{chr}.vcf.gz'),
		tbi = temp('freebayes/chr_split/{sample}/{sample}__{chr}.vcf.gz.tbi'),
		filteredvcf = temp('freebayes/chr_split/{sample}/{sample}__{chr}.filtered.vcf.gz'),
		filteredtbi = temp('freebayes/chr_split/{sample}/{sample}__{chr}.filtered.vcf.gz.tbi'),
		phasedvcf = temp('freebayes/chr_split/{sample}/{sample}__{chr}.phased.vcf.gz'),
		phasedvcf_tbi = temp('freebayes/chr_split/{sample}/{sample}__{chr}.phased.vcf.gz.tbi')
	threads: 32
	shell:
		"""
		module load {config[freebayes_version]}
		module load {config[vcflib_version]}
		module load {config[samtools_version]}
		module load {config[vt_version]}
		freebayes-parallel {config[freebayes_wgs_region]}.{wildcards.chr} {threads} -f {config[ref_genome]} \
			--limit-coverage 500 {input.bam} --min-alternate-fraction 0.05 \
			--min-mapping-quality 1 --genotype-qualities --strict-vcf --use-mapping-quality \
			| bcftools norm --multiallelics -any --output-type v - \
			| vt decompose_blocksub -p -m -d 2 - \
			| bcftools norm --check-ref s --fasta-ref {config[ref_genome]} --output-type v - \
			| bcftools norm -d none --output-type v - \
			| bgzip -f > {output.vcf}
		sleep 2
		tabix -f -p vcf {output.vcf}
		vcffilter -f "QUAL > 15 & QA / AO > 15 & SAF > 0 & SAR > 0 & RPR > 0 & RPL > 0 & AO > 2 & DP > 3" {output.vcf} | bgzip -f > {output.filteredvcf}
		sleep 2
		tabix -f -p vcf {output.filteredvcf}
		module unload {config[freebayes_version]}
		module unload {config[vcflib_version]}
		module unload {config[vt_version]}
		module load {config[whatshap_version]}
		whatshap phase --reference {config[ref_genome]} --indels {output.filteredvcf} {input.bam} | bgzip -f > {output.phasedvcf}
		tabix -f -p vcf {output.phasedvcf}
		"""

# try rtg eval wihtout vt and see how it goes.
	# | vt decompose -s - \
	# | vt normalize -m -r {config[ref_genome]} - \
	#
#2/11/20: --min-coverage 3, UUAL > 0
#freebayes high sensitivity and benchmark: https://europepmc.org/article/PMC/6500473
# -g "GQ > 1"
#			| sed -e "s|1/.|0/1|" -e "s|./1|0/1|" \

rule cat_fbvcfs:
	input:
		vcf = chr_fbVCF_to_single_fbVCF,
		filteredvcf = chr_fbVCF_to_single_fbVCF_filtered,
		phasedvcf = chr_fbVCF_to_single_fbVCF_phased
	output:
		vcf = 'freebayes/{sample}.vcf.gz',
		filteredvcf = 'freebayes/{sample}.filtered.vcf.gz',
		phasedvcf = 'freebayes/{sample}.phased.vcf.gz'
	threads: 8
	shell:
		"""
		module load {config[samtools_version]}
		bcftools concat --threads {threads} --output-type z {input.vcf} > {output.vcf}
		tabix -f -p vcf {output.vcf}
		bcftools concat --threads {threads} --output-type z {input.filteredvcf} > {output.filteredvcf}
		tabix -f -p vcf {output.filteredvcf}
		bcftools concat --threads {threads} --output-type z {input.phasedvcf} > {output.phasedvcf}
		tabix -f -p vcf {output.phasedvcf}
		"""

rule merge_freebayes:
	input:
		vcf = expand('freebayes/{sample}.phased.vcf.gz', sample=list(SAMPLE_LANEFILE.keys()))
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
		# bcftools merge --merge none --output-type z --threads {threads} {input.vcf} \
		#  	> prioritization/{config[analysis_batch_name]}.freebayes.vcf.gz

#sort necessary before mark_dup? it was sorted during the merge lane_bam step.
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

# if config['cram']  == 'TRUE':
# 	rule picard_merge_bams:
# # merge chr split bams into one bam per sample
# 		input:
# 			chr_bam_to_single_bam
# 		output:
# 			temp('sample_bam/{sample}.recalibrated.bam')
# 		threads: 2
# 		shell:
# 			"""
# 			module load {config[picard_version]}
# 			cat_inputs_i=""
# 			for bam in {input}; do
# 				cat_inputs_i+="I=$bam "; done
# 			java -Xmx15g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
# 				MergeSamFiles \
# 				$cat_inputs_i \
# 				O={output}
# 			"""
# else:
# 	rule picard_merge_bams_index:
# # merge chr split bams into one bam per sample
# 		input:
# 			chr_bam_to_single_bam
# 		output:
# 			bam = 'sample_bam/{sample}.recalibrated.bam',
# 			bai1 = 'sample_bam/{sample}.recalibrated.bai',
# 			bai2 = 'sample_bam/{sample}.recalibrated.bam.bai'
# 		threads: 2
# 		shell:
# 			"""
# 			module load {config[picard_version]}
# 			cat_inputs_i=""
# 			for bam in {input}; do
# 				cat_inputs_i+="I=$bam "; done
# 			java -Xmx15g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
# 				MergeSamFiles \
# 				SORT_ORDER=coordinate \
# 				CREATE_INDEX=true \
# 				$cat_inputs_i \
# 				O={output.bam}
# 			cp {output.bai1} {output.bai2}
# 			"""
#java -Xmx15g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
#	BuildBamIndex \
#	INPUT={output.bam} \
#	OUTPUT={output.bai}


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


##freebayes
#freebayes-parallel <(fasta_generate_regions.py ref.fa.fai 100000) 36 -f ref.fa aln.bam >out.vcf"


# rule scramble:
# 	input:
# 		bam = 'sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.bam',
# 		bai = 'sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.bai'
# 	output:
# 		cluster = temp('scramble/chr_split/{sample}/{sample}__{chr}.cluster.txt'),
# 		mei = 'scramble/chr_split/{sample}/{sample}__{chr}.mei.txt'
# 	shell:
# 		"""
# 		module load scramble
# 		scramble cluster_identifier {input.bam} > {output.cluster}
# 		scramble Rscript --vanilla /app/cluster_analysis/bin/SCRAMble-MEIs.R \
# 			--out-name ${{PWD}}/{output.mei} \
# 			--cluster-file ${{PWD}}/{output.cluster} \
# 			--install-dir /app/cluster_analysis/bin \
# 			--mei-refs /app/cluster_analysis/resources/MEI_consensus_seqs.fa
# 		"""
#
# #localrules: scramble_annotation
# rule scramble_annotation:
# 	input:
# 		mei = 'scramble/chr_split/{sample}/{sample}__{chr}.mei.txt'
# 	output:
# 		avinput = temp('scramble/chr_split/{sample}/{sample}__{chr}.avinput'),
# 		annovar = temp('scramble/chr_split/{sample}/{sample}__{chr}.hg19_multianno.txt'),
# 		annovarR = temp('scramble/chr_split/{sample}/{sample}__{chr}.forR.txt'),
# 		anno = 'scramble/chr_split/{sample}/{sample}__{chr}.scramble.tsv'
# 	shell:
# 		"""
# 		module load {config[R_version]}
# 		module load {config[annovar_version]}
# 		if [[ $(wc -l {input.mei} | cut -d " " -f 1) == 1 ]]
# 		then
# 			touch {output.avinput}
# 			touch {output.annovar}
# 			touch {output.annovarR}
# 			touch {output.anno}
# 		else
# 			cut -f 1 {input.mei} | awk -F ":" 'BEGIN{{OFS="\t"}} NR>1 {{print $1,$2,$2,"0","-"}}' > {output.avinput}
# 			table_annovar.pl {output.avinput} \
# 				$ANNOVAR_DATA/hg19 \
# 				-buildver hg19 \
# 				-remove \
# 				-out scramble_anno/{wildcards.sample} \
# 				--protocol refGene \
# 				-operation  g \
# 				--argument '-splicing 100 -hgvs' \
# 				--polish -nastring . \
# 				--thread 1
# 			awk -F"\t" 'BEGIN{{OFS="\t"}} NR==1 {{print "Func_refGene","Gene","Intronic","AA"}} NR>1 {{print $6,$7,$8,$10}}' {output.annovar} | paste {input.mei} - > {output.annovarR}
# 			Rscript /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/scramble_CHRanno.R {output.annovarR} {config[SCRAMBLEdb]} {config[OGL_Dx_research_genes]} {config[HGMDtranscript]} {output.anno} {wildcards.sample}
# 		fi
# 		"""
#
# localrules: cat_scramble
# rule cat_scramble:
# 	input:
# 		tsv = chr_scramble_to_single_scramble
# 	output:
# 		tsv = temp('scramble_anno/{sample}.scramble.tsv'),
# 		xlsx = 'scramble_anno/{sample}.scramble.xlsx'
# 	shell:
# 		"""
# 		module load {config[R_version]}
# 		cat {input} | grep -v "^eyeGene" > {output.tsv}
# 		sed -i "1 i\eyeGene\tInsertion\tMEI_Family\tInsertion_Direction\tClipped_Reads_In_Cluster\tAlignment_Score\tAlignment_Percent_Length\tAlignment_Percent_Identity\tClipped_Sequence\tClipped_Side\tStart_In_MEI\tStop_In_MEI\tpolyA_Position\tpolyA_Seq\tpolyA_SupportingReads\tTSD\tTSD_length\tpanel_class\tGene\tIntronic\tAA\tclassification\tpopAF\tsample\tnote" {output.tsv}
# 		Rscript /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/scramble_sortanno.R {output.tsv} {output.xlsx}
# 		"""
