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

if config['inputFileType'].upper() in ['BAM', 'CRAM']:
	def rg(wildcards):
		# returns the read group given in the config['metadata_file']
		lane_file = str(wildcards)
		rg_out = str(LANEFILE_READGROUP[str(SAMPLE_LANEFILE[lane_file][0])]).replace("['", "").replace("']","")
		return(rg_out)
else:
	def rg(wildcards):
		# returns the read group given in the config['metadata_file']
		lane_file = str(wildcards)
		rg_out = str(LANEFILE_READGROUP[lane_file + config['lane_pair_delim'][0] + '.fastq.gz'][0])
		return(rg_out)

# import CREST hg19 regions, chr1 to chrY
REGIONS_file = config['regions']
if '/home/$USER' in REGIONS_file:
	REGIONS_file = os.environ['HOME'] + REGIONS_file.split('$USER')[-1]
REGIONS = open(REGIONS_file).readlines()
REGIONS = [r.strip() for r in REGIONS]

CHRS=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT_contigs"]
MT_CONTIGS="MT GL000207.1 GL000226.1 GL000229.1 GL000231.1 GL000210.1 GL000239.1 GL000235.1 GL000201.1 GL000247.1 GL000245.1 GL000197.1 GL000203.1 GL000246.1 GL000249.1 GL000196.1 GL000248.1 GL000244.1 GL000238.1 GL000202.1 GL000234.1 GL000232.1 GL000206.1 GL000240.1 GL000236.1 GL000241.1 GL000243.1 GL000242.1 GL000230.1 GL000237.1 GL000233.1 GL000204.1 GL000198.1 GL000208.1 GL000191.1 GL000227.1 GL000228.1 GL000214.1 GL000221.1 GL000209.1 GL000218.1 GL000220.1 GL000213.1 GL000211.1 GL000199.1 GL000217.1 GL000216.1 GL000215.1 GL000205.1 GL000219.1 GL000224.1 GL000223.1 GL000195.1 GL000212.1 GL000222.1 GL000200.1 GL000193.1 GL000194.1 GL000225.1 GL000192.1 NC_007605"

wildcard_constraints:
	sample='|'.join(list(SAMPLE_LANEFILE.keys())),
	chr = '|'.join(CHRS),
	lane = '|'.join(list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in [y for sub in list(SAMPLE_LANEFILE.values()) for y in sub]]))),
	region = '|'.join(REGIONS)

rule all:
	input:
		# expand('CRESTanno/{sample}.predSV.xlsx', sample=list(SAMPLE_LANEFILE.keys())) if config['CREST'] == 'TRUE' else 'dummy.txt',
		expand('gvcfs/{sample}.g.vcf.gz', sample=list(SAMPLE_LANEFILE.keys())) if config['GATKgvcf'] == 'TRUE' else 'dummy.txt',
		# expand('recal_bam/{sample}.recal.bam', sample=list(SAMPLE_LANEFILE.keys())) if config['recal_bam'] == 'TRUE' else 'dummy.txt',
		expand('cram/{sample}.cram', sample=list(SAMPLE_LANEFILE.keys())) if config['cram'] == 'TRUE' else expand('bam/{sample}.bam', sample=list(SAMPLE_LANEFILE.keys())),
		# 'GATK_metrics/multiqc_report' if config['multiqc'] == 'TRUE' else 'dummy.txt',
		'fastqc/multiqc_report' if config['multiqc'] == 'TRUE' else 'dummy.txt',
		# expand('picardQC/{sample}.insert_size_metrics.txt', sample=list(SAMPLE_LANEFILE.keys())) if config['picardQC'] == 'TRUE' else 'dummy.txt',
		'deepvariant/deepvariantVcf.merge.done.txt' if config['deepvariant'] == 'TRUE' else 'dummy.txt',
		'prioritization/dv_fb.merge.done.txt' if config['freebayes_phasing'] == 'TRUE' else 'dummy.txt',
		expand('coverage/{sample}.coverage.xlsx', sample=list(SAMPLE_LANEFILE.keys())) if config['coverage'] == 'TRUE' else 'dummy.txt',
		expand('manta/manta.{sample}.annotated.tsv', sample=list(SAMPLE_LANEFILE.keys())),
		expand('scramble_anno/{sample}.scramble.xlsx', sample=list(SAMPLE_LANEFILE.keys())) if config['SCRAMble'] == 'TRUE' else 'dummy.txt'


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

if config['inputFileType'] == 'single_lane_fastq':
	rule align_markdup:
		input:
			expand('fastq/{{lane}}{pair}.fastq.gz', pair = config['lane_pair_delim'])
		output:
			bam = temp('lane_bam/{lane}.bam'),
			bai = temp('lane_bam/{lane}.bam.bai')
		params:
			read_group = rg
		threads: 16
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
			bai = lambda wildcards: expand('lane_bam/{lane}.bam.bai', lane = list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in SAMPLE_LANEFILE[wildcards.sample]])))
		output:
			bam = temp('sample_bam/{sample}.markDup.bam'),
			bai = temp('sample_bam/{sample}.markDup.bai')
		shell:
			"""
			cp -p -l {input.bam} {output.bam}
			cp -p -l {input.bai} {output.bai}
			"""
elif config['inputFileType'].upper() in ['BAM', 'CRAM']:
	rule realign:
		input:
			lambda wildcards: join('old_bam/', str(SAMPLE_LANEFILE[wildcards.sample][0]))
		output:
			bam = temp('sample_bam/{sample}.markDup.bam'),
			bai = temp('sample_bam/{sample}.markDup.bai')
		threads: 16
		params:
			read_group = rg
		shell:
			"""
			export TMPDIR=/lscratch/$SLURM_JOB_ID
			module load {config[bazam_version]}
			module load {config[bwa-mem2_version]} {config[samblaster_version]} {config[sambamba_version]}
 			case "{input}" in
				*bam)
					java -Xmx12g -jar $BAZAMPATH/bazam.jar -bam {input} \
					| bwa-mem2 mem -t {threads} -K 100000000 -M -Y -B 4 -O 6 -E 1 -p -R {params.read_group} {config[bwa-mem2_ref]} - \
			 		| samblaster -M --addMateTags --quiet \
					| sambamba sort -u --tmpdir=/lscratch/$SLURM_JOB_ID -t {threads} -o {output.bam} \
						<(sambamba view -S -f bam -l 0 -t $SLURM_CPUS_PER_TASK /dev/stdin)
					mv {output.bam}.bai {output.bai}
					#sambamba index -t {threads} {output.bam} {output.bai}
					;;
				*cram)
					java -Xmx12g -Dsamjdk.reference_fasta={config[old_cram_ref]} -jar $BAZAMPATH/bazam.jar -bam {input} \
					| bwa-mem2 mem -t {threads} -K 100000000 -M -Y -B 4 -O 6 -E 1 -p -R {params.read_group} {config[bwa-mem2_ref]} - \
			 		| samblaster -M --addMateTags --quiet \
					| sambamba sort -u --tmpdir=/lscratch/$SLURM_JOB_ID -t {threads} -o {output.bam} \
						<(sambamba view -S -f bam -l 0 -t $SLURM_CPUS_PER_TASK /dev/stdin)
					mv {output.bam}.bai {output.bai}
					;;
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
			bam = temp('lane_bam/{lane}.bam'),
			bai = temp('lane_bam/{lane}.bam.bai')
		params:
			read_group = rg
		threads: 16
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
	rule merge_lane_bam:
		input:
			lambda wildcards: expand('lane_bam/{lane}.bam', lane = list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in SAMPLE_LANEFILE[wildcards.sample]]))),
			lambda wildcards: expand('lane_bam/{lane}.bam.bai', lane = list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in SAMPLE_LANEFILE[wildcards.sample]])))
		output:
			merged_bam = temp('sample_bam/{sample}.merged.bam'),
			merged_bai = temp('sample_bam/{sample}.merged.bam.bai'),
			bam = temp('sample_bam/{sample}.markDup.bam'),
			bai = temp('sample_bam/{sample}.markDup.bai'),
			metrics = 'sample_bam/{sample}.duplication_metrics.txt'
		threads: 16
		shell:
			"""
			module load {config[sambamba_version]} {config[picard_version]}
			case "{input.bam}" in
				*\ *)
					sambamba merge -t {threads} {output.merged_bam} {input.bam}
					java -Xmx16g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
						MarkDuplicates \
						--INPUT {output.merged_bam} \
						--OUTPUT {output.bam} \
						--METRICS_FILE {output.metrics} \
						--COMPRESSION_LEVEL 6 \
						--CREATE_INDEX true \
						--ASSUME_SORT_ORDER coordinate
					;;
				*)
					touch {output.merged_bam} {output.merged_bai} {output.metrics}
					cp -p -l {input.bam} {output.bam}
					cp -p -l {input.bai} {output.bai}
					;;
			esac
			"""

# rule merge_lane_bam_hg19:
# 	input:
# 		lambda wildcards: expand('lane_bam/hg19bam/hg19.{lane}.realigned.bam', lane = list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in SAMPLE_LANEFILE[wildcards.sample]])))
# 	output:
# 		bam = temp('sample_bam/hg19bam/hg19.{sample}.bam'),
# 		bai = temp('sample_bam/hg19bam/hg19.{sample}.bai')
# 	shell:
# 		"""
# 		module load {config[picard_version]}
# 		picard_i=""
# 		for bam in {input}; do
# 			picard_i+=" I=$bam"
# 		done
# 		java -Xmx8g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
# 			MergeSamFiles \
# 			TMP_DIR=/scratch/$SLURM_JOB_ID \
# 			$picard_i \
# 			O={output.bam} \
# 			SORT_ORDER=coordinate \
# 			CREATE_INDEX=true
# 		"""

# if multiple sets of fastq/bam provided for a sample, now merge together
# rule merge_lane_bam:
# 	input:
# 		lambda wildcards: expand('lane_bam/{lane}.bam', lane = list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in SAMPLE_LANEFILE[wildcards.sample]])))
# 	output:
# 		bam = temp('sample_bam/{sample}/{sample}.b37.bam'),
# 		bai = temp('sample_bam/{sample}/{sample}.b37.bai')
# 	shell:
# 		"""
# 		module load {config[picard_version]}
# 		picard_i=""
# 		for bam in {input}; do
# 			picard_i+=" I=$bam"
# 		done
# 		java -Xmx8g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
# 			MergeSamFiles \
# 			TMP_DIR=/scratch/$SLURM_JOB_ID \
# 			$picard_i \
# 			O={output.bam} \
# 			SORT_ORDER=coordinate \
# 			CREATE_INDEX=true
# 		"""

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
		module load {config[picard_version]} {config[R_version]}
		java -Xmx8g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			CollectInsertSizeMetrics \
			-TMP_DIR /lscratch/$SLURM_JOB_ID \
			--INPUT {input.bam} \
			-O {output.insert_size_metrics} \
		    -H {output.insert_size_histogram} \
		    -M 0.5
		java -Xmx8g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			CollectAlignmentSummaryMetrics \
			-TMP_DIR /lscratch/$SLURM_JOB_ID \
			--INPUT {input.bam} \
			-R {config[ref_genome]} \
			--METRIC_ACCUMULATION_LEVEL SAMPLE \
			--METRIC_ACCUMULATION_LEVEL READ_GROUP \
			-O {output.alignment_metrics}
		"""

# rule picard_mark_dups_allchr:
# # Mark duplicate reads
# 	input:
# 		bam = 'sample_bam/{sample}/{sample}.b37.bam',
# 		bai = 'sample_bam/{sample}/{sample}.b37.bai'
# 	output:
# 		bam = temp('sample_bam/{sample}.markDup.bam'),
# 		bai = temp('sample_bam/{sample}.markDup.bai'),
# 		metrics = temp('GATK_metrics/{sample}.markDup.metrics')
# 	threads: 16
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
		mosdepth -t {threads} --no-per-base --by {config[bed]} --use-median --mapq 1 --thresholds 10,20,30 \
			{wildcards.sample} {input.bam}
		mv {wildcards.sample}.thresholds.bed.gz* coverage/mosdepth/.
		mv {wildcards.sample}.mosdepth* coverage/mosdepth/.
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

#if too slow, then seperate by chr as above.
#localrules: scramble
rule scramble:
	input:
		bam = 'sample_bam/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}.markDup.bai'
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
#--bind /gpfs,/spin1,/data,/lscratch,/scratch,/fdb
#localrules: scramble_annotation
rule scramble_annotation:
	input:
		mei = 'scramble/{sample}_MEIs.txt',
		deletion = 'scramble/{sample}_PredictedDeletions.txt'
	output:
		avinput = temp('scramble_anno/{sample}.avinput'),
		annovar = temp('scramble_anno/{sample}.hg19_multianno.txt'),
		annovarR = temp('scramble_anno/{sample}.forR.txt'),
		anno = 'scramble_anno/{sample}.scramble.xlsx',
		del_anno = 'scramble_anno/{sample}.scramble.del.tsv'
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
		if [[ $(wc -l {input.deletion} | cut -d " " -f 1) == 1 ]]
		then
			touch {output.del_anno}
		else
			module load {config[annotsv_version]}
			tail -n +2 {input.deletion} | awk -F"\t" 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,"DEL"}}' > {input.deletion}.bed
			AnnotSV -SVinputFile {input.deletion}.bed -SVinputInfo 0 -svtBEDcol 4 -outputFile {output.del_anno}
		fi
		"""
#Try this strategy,if not working well, consider split by chr as in the GATK.
rule freebayes_phasing:
	input:
		bam = 'sample_bam/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}.markDup.bai'
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
		freebayes-parallel {config[freebayes_exome_region]} {threads} -f {config[ref_genome]} \
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

		# """
		# module load {config[vcftools_version]}
		# vcf-merge {input.vcf} | bgzip -c > prioritization/{config[analysis_batch_name]}.freebayes.vcf.gz
		# sleep 60
		# tabix -f -p vcf prioritization/{config[analysis_batch_name]}.freebayes.vcf.gz
		# touch {output}
		# """
#	module load {config[samtools_version]} - samtools automaticcally loaded when loading vcftools.

#
# rule create_freebayes_region:
# 	input:
# 		bam = expand('sample_bam/{sample}.markDup.bam', sample=list(SAMPLE_LANEFILE.keys())[0]),
# 		bai = expand('sample_bam/{sample}.markDup.bai', sample=list(SAMPLE_LANEFILE.keys())[0])
# 	output:
# 		'freebayes/freebayes.exome.10000.region'
# 	shell:
# 		"""
# 		module load {config[freebayes_version]}
# 		module load {config[samtools_version]}
# 		samtools depth {input.bam} | coverage_to_regions.py {config[ref_genome]}.fai 10000 > {output}
# 		"""
rule deepvariant:
	input:
		bam = 'sample_bam/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}.markDup.bai'
	output:
		vcf = 'deepvariant/vcf/{sample}.dv.vcf.gz',
		gvcf = 'deepvariant/gvcf/{sample}.dv.g.vcf.gz',
		filteredvcf = temp('deepvariant/vcf/{sample}.dv.filtered.vcf.gz'),
		filteretbi = temp('deepvariant/vcf/{sample}.dv.filtered.vcf.gz.tbi'),
		phasedvcf = 'deepvariant/vcf/{sample}.dv.phased.vcf.gz',
		phasedtbi = 'deepvariant/vcf/{sample}.dv.phased.vcf.gz.tbi'
	threads: 32
	shell:
		"""
		module load {config[deepvariant_version]}
		PROJECT_WD=$PWD
		N_SHARDS="16"
		mkdir -p /lscratch/$SLURM_JOB_ID/{wildcards.sample}
		WORK_DIR=/lscratch/$SLURM_JOB_ID/{wildcards.sample}
		cd $WORK_DIR
		run_deepvariant --model_type WES --num_shards $N_SHARDS \
			--ref {config[ref_genome]} \
			--regions {config[padded_bed]} \
			--reads $PROJECT_WD/{input.bam} \
			--output_vcf $PROJECT_WD/{output.vcf} \
			--output_gvcf $PROJECT_WD/{output.gvcf} \
			--sample_name {wildcards.sample} \
			--intermediate_results_dir $WORK_DIR \
			--call_variants_extra_args="use_openvino=true"
		cd $PROJECT_WD
		module unload {config[deepvariant_version]}
		module load {config[samtools_version]}
		bcftools norm --multiallelics -any --output-type v {output.vcf} \
			| bcftools norm -d none --output-type v - \
			| bcftools filter --include 'FILTER="PASS"' --output-type z --output {output.filteredvcf}
		sleep 2
		tabix -f -p vcf {output.filteredvcf}
		module load {config[whatshap_version]}
		whatshap phase --reference {config[ref_genome]} --indels {output.filteredvcf} {input.bam} | bgzip -f > {output.phasedvcf}
		tabix -f -p vcf {output.phasedvcf}
		"""

rule merge_deepvariant_vcf:
	input:
		vcf = expand('deepvariant/vcf/{sample}.dv.phased.vcf.gz', sample=list(SAMPLE_LANEFILE.keys())),
		tbi = expand('deepvariant/vcf/{sample}.dv.phased.vcf.gz.tbi', sample=list(SAMPLE_LANEFILE.keys()))
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
				sleep 2
				tabix -f -p vcf deepvariant/{config[analysis_batch_name]}.dv.phased.vcf.gz
				;;
			*)
				cp -p -l {input.vcf} deepvariant/{config[analysis_batch_name]}.dv.phased.vcf.gz
				cp -p -l {input.tbi} deepvariant/{config[analysis_batch_name]}.dv.phased.vcf.gz.tbi
				;;
		esac
		touch {output}
		"""

localrules: merge_dv_fb_vcfs
rule merge_dv_fb_vcfs:
	input:
		'deepvariant/deepvariantVcf.merge.done.txt',
		'freebayes/freebayes.merge.done.txt'
	output:
		'prioritization/dv_fb.merge.done.txt'
	threads: 8
	shell:
		"""
		module load {config[samtools_version]}
		WORK_DIR=/lscratch/$SLURM_JOB_ID
		bcftools isec -p $WORK_DIR --collapse none -Ov \
			deepvariant/{config[analysis_batch_name]}.dv.phased.vcf.gz \
			freebayes/{config[analysis_batch_name]}.freebayes.vcf.gz
		rm $WORK_DIR/0003.vcf &
		bcftools annotate --threads {threads} --set-id 'dv' -x FORMAT/VAF,FORMAT/PL \
			--no-version $WORK_DIR/0000.vcf -Oz -o $WORK_DIR/dv.vcf.gz
		rm $WORK_DIR/0000.vcf &
		bcftools annotate --threads {threads} --set-id 'fb' -x INFO,FORMAT/RO,FORMAT/QR,FORMAT/AO,FORMAT/QA,FORMAT/GL \
			--no-version $WORK_DIR/0001.vcf -Oz -o $WORK_DIR/fb.vcf.gz
		rm $WORK_DIR/0001.vcf &
		bcftools annotate --threads {threads} --set-id 'dvFb' -x FORMAT/VAF,FORMAT/PL \
			--no-version $WORK_DIR/0002.vcf -Oz -o $WORK_DIR/dvFb.vcf.gz
		rm $WORK_DIR/0002.vcf &
		tabix -f -p vcf $WORK_DIR/dv.vcf.gz
		tabix -f -p vcf $WORK_DIR/fb.vcf.gz
		tabix -f -p vcf $WORK_DIR/dvFb.vcf.gz
		bcftools concat --threads {threads} -a --rm-dups none --no-version \
			$WORK_DIR/dvFb.vcf.gz $WORK_DIR/dv.vcf.gz $WORK_DIR/fb.vcf.gz \
			-Oz -o prioritization/{config[analysis_batch_name]}.vcf.gz
		tabix -f -p vcf prioritization/{config[analysis_batch_name]}.vcf.gz
		touch {output}
		"""

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
	threads: 8
	shell:
		"""
		module load {config[manta_version]}
		mkdir -p /lscratch/$SLURM_JOB_ID/manta/{wildcards.sample}
		RUNDIR="/lscratch/$SLURM_JOB_ID/manta/{wildcards.sample}"
		configManta.py --referenceFasta {config[ref_genome]} \
			--exome --runDir $RUNDIR --bam {input.bam}
		$RUNDIR/runWorkflow.py -m local -j {threads} -g $((SLURM_MEM_PER_NODE / 1024))
		cp $RUNDIR/results/variants/diploidSV.vcf.gz manta/{wildcards.sample}.diploidSV.vcf.gz
		cp $RUNDIR/results/variants/diploidSV.vcf.gz.tbi manta/{wildcards.sample}.diploidSV.vcf.gz.tbi
		module load {config[annotsv_version]}
		AnnotSV -SVinputFile $RUNDIR/results/variants/diploidSV.vcf.gz \
			-SVinputInfo 1 -genomeBuild GRCh37\
			-outputDir $RUNDIR \
			-outputFile $RUNDIR/manta.{wildcards.sample}.annotated.tsv
		cp $RUNDIR/manta.{wildcards.sample}.annotated.tsv {output}
		"""


rule split_bam_by_chr:
	input:
		bam = 'sample_bam/{sample}/{sample}.b37.bam',
		bai = 'sample_bam/{sample}/{sample}.b37.bai'
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

rule picard_mark_dups:
# Mark duplicate reads
	input:
		'sample_bam/chr_split/{sample}/{sample}__{chr}.bam'
	output:
		clean_bam = temp('sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.bam'),
		sorted_bam = temp('sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.bam'),
		bam = temp('sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.bam'),
		bai = temp('sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.bai'),
		metrics = 'GATK_metrics/{sample}__{chr}.markDup.metrics'
	threads: 2
	shell:
		"""
		module load {config[picard_version]}
		java -Xmx60g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			CleanSam \
			TMP_DIR=/lscratch/$SLURM_JOB_ID \
			INPUT={input} \
			OUTPUT={output.clean_bam}
		java -Xmx60g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			FixMateInformation \
			SORT_ORDER=coordinate \
			INPUT={output.clean_bam} \
			OUTPUT={output.sorted_bam}
		java -Xmx60g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
			MarkDuplicates \
			INPUT={output.sorted_bam}	 \
			OUTPUT={output.bam} \
			METRICS_FILE={output.metrics} \
			CREATE_INDEX=true
		"""

# rule picard_clean_sam:
# # "Soft-clipping beyond-end-of-reference alignments and setting MAPQ to 0 for unmapped reads"
# 	input:
# 		'sample_bam/chr_split/{sample}/{sample}__{chr}.bam'
# 	output:
# 		temp('sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.bam')
# 	threads: 2
# 	shell:
# 		"""
# 		module load {config[picard_version]}
# 		java -Xmx60g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
# 			CleanSam \
# 			TMP_DIR=/lscratch/$SLURM_JOB_ID \
# 			INPUT={input} \
# 			OUTPUT={output}
# 		"""
#
# rule picard_fix_mate_information:
# # "Verify mate-pair information between mates and fix if needed."
# # also coord sorts
# 	input:
# 		'sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.bam'
# 	output:
# 		temp('sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.bam')
# 	threads: 2
# 	shell:
# 		"""
# 		module load {config[picard_version]}
# 		java -Xmx60g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
# 		FixMateInformation \
# 			SORT_ORDER=coordinate \
# 			INPUT={input} \
# 			OUTPUT={output}
# 		"""
#
# rule picard_mark_dups:
# # Mark duplicate reads
# 	input:
# 		'sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.bam'
# 	output:
# 		bam = temp('sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.bam'),
# 		bai = temp('sample_bam/chr_split/{sample}/{sample}__{chr}.CleanSam.sorted.markDup.bai'),
# 		metrics = 'GATK_metrics/{sample}__{chr}.markDup.metrics'
# 	threads: 2
# 	shell:
# 		"""
# 		module load {config[picard_version]}
# 		java -Xmx60g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
# 			MarkDuplicates \
# 			INPUT={input} \
# 			OUTPUT={output.bam} \
# 			METRICS_FILE={output.metrics} \
# 			CREATE_INDEX=true
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

localrules: picard_merge_gvcfs
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
