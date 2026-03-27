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
# for i in SAMPLE_LANEFILE:
# 	print (i, SAMPLE_LANEFILE[i], len(SAMPLE_LANEFILE[i]))
# for i in LANEFILE_READGROUP:
# 	print (i, LANEFILE_READGROUP[i], len(LANEFILE_READGROUP[i]))

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

# if config['analysis_batch_name'] == 'YYYYMMDD':
# 	currentDT = datetime.datetime.now()
# 	config['analysis_batch_name'] = currentDT.strftime("%Y%m%d")

if config['inputFileType'].upper() in ['BAM', 'CRAM', 'BAM-REALIGN', 'CRAM-REALIGN']:
	def rg(wildcards):
		# returns the read group given in the config['metadata_file']
		lane_file = str(wildcards)
		rg_out = str(LANEFILE_READGROUP[str(SAMPLE_LANEFILE[lane_file][0])]).replace("['", "").replace("']","")
		return(rg_out)
else:
	def rg(wildcards):
		# returns the read group given in the config['metadata_file']
		lane_file = str(wildcards)
		rg_out = str(LANEFILE_READGROUP[lane_file + config['lane_pair_delim'][0] + '.gz'][0])
		return(rg_out)


# import CREST hg19 regions
#REGIONS_file = config['regions']
#if '/home/$USER' in REGIONS_file:
#	REGIONS_file = os.environ['HOME'] + REGIONS_file.split('$USER')[-1]
#REGIONS = open(REGIONS_file).readlines()
#REGIONS = [r.strip() for r in REGIONS]

if config['genomeBuild'].upper() in ['GRCH37', 'HG19']:
	config['ref_genome'] = '/data/OGL/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta'
	config['bwa-mem2_ref'] = '/data/OGL/resources/1000G_phase2_GRCh37/bwa-mem2/human_g1k_v37_decoy.fasta'
	config['SCRAMBLEdb'] = '/data/OGL/resources/SCRAMBLEvariantClassification.xlsx'
	CHRS=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT_contigs"]
	MT_CONTIGS = "MT GL000207.1 GL000226.1 GL000229.1 GL000231.1 GL000210.1 GL000239.1 GL000235.1 GL000201.1 GL000247.1 GL000245.1 GL000197.1 GL000203.1 GL000246.1 GL000249.1 GL000196.1 GL000248.1 GL000244.1 GL000238.1 GL000202.1 GL000234.1 GL000232.1 GL000206.1 GL000240.1 GL000236.1 GL000241.1 GL000243.1 GL000242.1 GL000230.1 GL000237.1 GL000233.1 GL000204.1 GL000198.1 GL000208.1 GL000191.1 GL000227.1 GL000228.1 GL000214.1 GL000221.1 GL000209.1 GL000218.1 GL000220.1 GL000213.1 GL000211.1 GL000199.1 GL000217.1 GL000216.1 GL000215.1 GL000205.1 GL000219.1 GL000224.1 GL000223.1 GL000195.1 GL000212.1 GL000222.1 GL000200.1 GL000193.1 GL000194.1 GL000225.1 GL000192.1 NC_007605"
elif config['genomeBuild'].upper() in ['GRCH38', 'HG38']:
	CHRS=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","MT_contigs"]
	MT_CONTIGS = "chrM chr1_KI270706v1_random chr1_KI270707v1_random chr1_KI270708v1_random chr1_KI270709v1_random chr1_KI270710v1_random chr1_KI270711v1_random chr1_KI270712v1_random chr1_KI270713v1_random chr1_KI270714v1_random chr2_KI270715v1_random chr2_KI270716v1_random chr3_GL000221v1_random chr4_GL000008v2_random chr5_GL000208v1_random chr9_KI270717v1_random chr9_KI270718v1_random chr9_KI270719v1_random chr9_KI270720v1_random chr11_KI270721v1_random chr14_GL000009v2_random chr14_GL000225v1_random chr14_KI270722v1_random chr14_GL000194v1_random chr14_KI270723v1_random chr14_KI270724v1_random chr14_KI270725v1_random chr14_KI270726v1_random chr15_KI270727v1_random chr16_KI270728v1_random chr17_GL000205v2_random chr17_KI270729v1_random chr17_KI270730v1_random chr22_KI270731v1_random chr22_KI270732v1_random chr22_KI270733v1_random chr22_KI270734v1_random chr22_KI270735v1_random chr22_KI270736v1_random chr22_KI270737v1_random chr22_KI270738v1_random chr22_KI270739v1_random chrY_KI270740v1_random chrUn_KI270302v1 chrUn_KI270304v1 chrUn_KI270303v1 chrUn_KI270305v1 chrUn_KI270322v1 chrUn_KI270320v1 chrUn_KI270310v1 chrUn_KI270316v1 chrUn_KI270315v1 chrUn_KI270312v1 chrUn_KI270311v1 chrUn_KI270317v1 chrUn_KI270412v1 chrUn_KI270411v1 chrUn_KI270414v1 chrUn_KI270419v1 chrUn_KI270418v1 chrUn_KI270420v1 chrUn_KI270424v1 chrUn_KI270417v1 chrUn_KI270422v1 chrUn_KI270423v1 chrUn_KI270425v1 chrUn_KI270429v1 chrUn_KI270442v1 chrUn_KI270466v1 chrUn_KI270465v1 chrUn_KI270467v1 chrUn_KI270435v1 chrUn_KI270438v1 chrUn_KI270468v1 chrUn_KI270510v1 chrUn_KI270509v1 chrUn_KI270518v1 chrUn_KI270508v1 chrUn_KI270516v1 chrUn_KI270512v1 chrUn_KI270519v1 chrUn_KI270522v1 chrUn_KI270511v1 chrUn_KI270515v1 chrUn_KI270507v1 chrUn_KI270517v1 chrUn_KI270529v1 chrUn_KI270528v1 chrUn_KI270530v1 chrUn_KI270539v1 chrUn_KI270538v1 chrUn_KI270544v1 chrUn_KI270548v1 chrUn_KI270583v1 chrUn_KI270587v1 chrUn_KI270580v1 chrUn_KI270581v1 chrUn_KI270579v1 chrUn_KI270589v1 chrUn_KI270590v1 chrUn_KI270584v1 chrUn_KI270582v1 chrUn_KI270588v1 chrUn_KI270593v1 chrUn_KI270591v1 chrUn_KI270330v1 chrUn_KI270329v1 chrUn_KI270334v1 chrUn_KI270333v1 chrUn_KI270335v1 chrUn_KI270338v1 chrUn_KI270340v1 chrUn_KI270336v1 chrUn_KI270337v1 chrUn_KI270363v1 chrUn_KI270364v1 chrUn_KI270362v1 chrUn_KI270366v1 chrUn_KI270378v1 chrUn_KI270379v1 chrUn_KI270389v1 chrUn_KI270390v1 chrUn_KI270387v1 chrUn_KI270395v1 chrUn_KI270396v1 chrUn_KI270388v1 chrUn_KI270394v1 chrUn_KI270386v1 chrUn_KI270391v1 chrUn_KI270383v1 chrUn_KI270393v1 chrUn_KI270384v1 chrUn_KI270392v1 chrUn_KI270381v1 chrUn_KI270385v1 chrUn_KI270382v1 chrUn_KI270376v1 chrUn_KI270374v1 chrUn_KI270372v1 chrUn_KI270373v1 chrUn_KI270375v1 chrUn_KI270371v1 chrUn_KI270448v1 chrUn_KI270521v1 chrUn_GL000195v1 chrUn_GL000219v1 chrUn_GL000220v1 chrUn_GL000224v1 chrUn_KI270741v1 chrUn_GL000226v1 chrUn_GL000213v1 chrUn_KI270743v1 chrUn_KI270744v1 chrUn_KI270745v1 chrUn_KI270746v1 chrUn_KI270747v1 chrUn_KI270748v1 chrUn_KI270749v1 chrUn_KI270750v1 chrUn_KI270751v1 chrUn_KI270752v1 chrUn_KI270753v1 chrUn_KI270754v1 chrUn_KI270755v1 chrUn_KI270756v1 chrUn_KI270757v1 chrUn_GL000214v1 chrUn_KI270742v1 chrUn_GL000216v2 chrUn_GL000218v1 chrEBV"
else:
	print("ref_genome is ", config['ref_genome'])
	print("bwa-mem2_ref is", config['bwa-mem2_ref'])
	print("SCRAMBLEdb is", config['SCRAMBLEdb'])
	CHRS=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT_contigs"]
	MT_CONTIGS = "MT GL000207.1 GL000226.1 GL000229.1 GL000231.1 GL000210.1 GL000239.1 GL000235.1 GL000201.1 GL000247.1 GL000245.1 GL000197.1 GL000203.1 GL000246.1 GL000249.1 GL000196.1 GL000248.1 GL000244.1 GL000238.1 GL000202.1 GL000234.1 GL000232.1 GL000206.1 GL000240.1 GL000236.1 GL000241.1 GL000243.1 GL000242.1 GL000230.1 GL000237.1 GL000233.1 GL000204.1 GL000198.1 GL000208.1 GL000191.1 GL000227.1 GL000228.1 GL000214.1 GL000221.1 GL000209.1 GL000218.1 GL000220.1 GL000213.1 GL000211.1 GL000199.1 GL000217.1 GL000216.1 GL000215.1 GL000205.1 GL000219.1 GL000224.1 GL000223.1 GL000195.1 GL000212.1 GL000222.1 GL000200.1 GL000193.1 GL000194.1 GL000225.1 GL000192.1 NC_007605"

wildcard_constraints:
	sample='|'.join(list(SAMPLE_LANEFILE.keys())),
	lane = '|'.join(list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in [y for sub in list(SAMPLE_LANEFILE.values()) for y in sub]])))
#	region = '|'.join(REGIONS)

rule all:
	input:
		expand('gvcfs/{sample}.g.vcf.gz', sample=list(SAMPLE_LANEFILE.keys())) if config['GATKgvcf'] == 'TRUE' else 'dummy.txt',
		#'GATK_metrics/multiqc_report' if config['multiqc'] == 'TRUE' else 'dummy.txt',
		'fastqc/multiqc_report' if config['multiqc'] == 'TRUE' else 'dummy.txt',
		expand('picardQC/{sample}.insert_size_metrics.txt', sample=list(SAMPLE_LANEFILE.keys())) if config['picardQC'] == 'TRUE' else 'dummy.txt',
		'CoNVaDING/progress2.done' if config['CoNVaDING'] == 'TRUE' else 'dummy.txt',
		#'deepvariant/deepvariantVcf.merge.done.txt' if config['deepvariant'] == 'TRUE' else 'dummy.txt',
		'prioritization/dv_fb.merge.done.txt' if config['freebayes_phasing'] == 'TRUE' else 'dummy.txt',
		'coverage/mean.coverage.done.txt' if config['coverage'] == 'TRUE' else 'dummy.txt',
		expand('bam/{sample}.cram', sample=list(SAMPLE_LANEFILE.keys())) if config['cram'] == 'TRUE' else expand('bam/{sample}.bam', sample=list(SAMPLE_LANEFILE.keys())),
		expand('scramble_anno/{sample}.scramble.tsv', sample=list(SAMPLE_LANEFILE.keys())) if config['SCRAMble'] == 'TRUE' else 'dummy.txt',
		expand('manta/manta.{sample}.annotated.tsv', sample=list(SAMPLE_LANEFILE.keys())),
		'bcmlocus/combine.bcmlocus.done.txt' if config['bcm_locus'] == 'TRUE' else 'dummy.txt' #expand('bcmlocus/{sample}.bcmlocus.txt', sample=list(SAMPLE_LANEFILE.keys()))
		#expand('AutoMap/{sample}/{sample}.HomRegions.annot.tsv', sample=list(SAMPLE_LANEFILE.keys()))

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
#sambamba default compression level is 6

if config['cutadapt'] == 'TRUE':
	rule trim_adatpor:
		input:
			expand('fastq/{{lane}}{pair}.gz', pair = config['lane_pair_delim'])
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
			bam = temp('lane_bam/{lane}.realigned.bam'),
			bai = temp('lane_bam/{lane}.realigned.bam.bai'),
		params:
			read_group = rg
		threads: 8
		shell:
			"""
			export TMPDIR=/lscratch/$SLURM_JOB_ID
			echo {params.read_group}
			module load {config[bwa-mem2_version]} {config[samblaster_version]} {config[sambamba_version]}
			bwa-mem2 mem -t {threads} -K 100000000 -M -Y -B 4 -O 6 -E 1 -R {params.read_group} \
				{config[bwa-mem2_ref]} {input} \
			 	| samblaster -M --acceptDupMarks --addMateTags --quiet \
				| sambamba sort -u --compression-level 6 --tmpdir=/lscratch/$SLURM_JOB_ID -t {threads} -o {output.bam} \
					<(sambamba view -S -f bam --compression-level 0 -t $SLURM_CPUS_PER_TASK /dev/stdin)
			"""
elif config['inputFileType'].upper() in ['BAM', 'CRAM']:
	rule realign:
		input:
			lambda wildcards: join('old_bam/', str(SAMPLE_LANEFILE[wildcards.sample][0]))
		output:
			bam = temp('sample_bam/{sample}/{sample}.markDup.bam'),
			bai = temp('sample_bam/{sample}/{sample}.markDup.bai')
		threads: 20
		params:
			read_group = rg
		shell:
			"""
			export TMPDIR=/lscratch/$SLURM_JOB_ID
			echo {params.read_group}
			module load {config[bazam_version]}
			module load {config[bwa-mem2_version]} {config[samblaster_version]} {config[sambamba_version]}
			BAMFILE={input}
			if [ -e {input}.bai ] || [ -e ${{BAMFILE%.bam}}.bai ] || [ -e {input}.crai ] || [ -e ${{BAMFILE%.cram}}.crai ] ; then
				echo "index present"
			else
				sambamba index -t $(({threads}-2)) {input}
			fi
 			case "{input}" in
				*bam)
					java -Xmx12g -jar $BAZAMPATH/bazam.jar -bam {input} \
					| bwa-mem2 mem -t 8 -K 100000000 -M -Y -B 4 -O 6 -E 1 -p -R {params.read_group} {config[bwa-mem2_ref]} - \
			 		| samblaster -M --addMateTags --quiet \
					| sambamba sort -u --tmpdir=/lscratch/$SLURM_JOB_ID -t 8 -o {output.bam} \
						<(sambamba view -S -f bam -l 0 -t 8 /dev/stdin)
					mv {output.bam}.bai {output.bai}
					#sambamba index -t {threads} {output.bam} {output.bai}
					;;
				*cram)
					java -Xmx12g -Dsamjdk.reference_fasta={config[old_cram_ref]} -jar $BAZAMPATH/bazam.jar -bam {input} \
					| bwa-mem2 mem -t 8 -K 100000000 -M -Y -B 4 -O 6 -E 1 -p -R {params.read_group} {config[bwa-mem2_ref]} - \
			 		| samblaster -M --addMateTags --quiet \
					| sambamba sort -u --tmpdir=/lscratch/$SLURM_JOB_ID -t 8 -o {output.bam} \
						<(sambamba view -S -f bam -l 0 -t 6 /dev/stdin)
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
			expand('fastq/{{lane}}{pair}.gz', pair = config['lane_pair_delim'])
		output:
			bam = temp('lane_bam/{lane}.realigned.bam'),
			bai = temp('lane_bam/{lane}.realigned.bam.bai')
		params:
			read_group = rg
		threads: 8
		shell:
			"""
			export TMPDIR=/lscratch/$SLURM_JOB_ID
			echo {params.read_group}
			module load {config[bwa-mem2_version]} {config[samblaster_version]} {config[sambamba_version]}
			bwa-mem2 mem -t {threads} -K 100000000 -M -Y -B 4 -O 6 -E 1 -R {params.read_group} \
				{config[bwa-mem2_ref]} {input} \
			 	| samblaster -M --acceptDupMarks --addMateTags --quiet \
				| sambamba sort -u --tmpdir=/lscratch/$SLURM_JOB_ID -t {threads} -o {output.bam} \
					<(sambamba view -S -f bam -l 0 -t $SLURM_CPUS_PER_TASK /dev/stdin)
			"""
#sambamba default comppression leve is 6. tested on 1/14/21, it generates *.bam.bai automatically.
#bwa-mem2 requires 32g mem with 8 threads.
#for WGS use bwa-mem2
			# module load {config[bwa-mem2_version]}
			# module load {config[samtools_version]};
			# bwa-mem2 mem -t {threads} -K 100000000 -M -B 4 -O 6 -E 1 -R {params.read_group} \

rule merge_lane_bam:
	input:
		bam = lambda wildcards: expand('lane_bam/{lane}.realigned.bam', lane = list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in SAMPLE_LANEFILE[wildcards.sample]]))),
		bai = lambda wildcards: expand('lane_bam/{lane}.realigned.bam.bai', lane = list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in SAMPLE_LANEFILE[wildcards.sample]])))
	output:
		bam = temp('sample_bam/{sample}/{sample}.markDup.bam'),
		bai = temp('sample_bam/{sample}/{sample}.markDup.bai')
	threads: 8
	shell:
		"""
		module load {config[sambamba_version]}
		case "{input.bam}" in
			*\ *)
				sambamba merge -t {threads} /lscratch/$SLURM_JOB_ID/{wildcards.sample}.bam {input.bam}
				sambamba markdup -t {threads} -l 6 \
					--tmpdir=/lscratch/$SLURM_JOB_ID/{wildcards.sample} \
					/lscratch/$SLURM_JOB_ID/{wildcards.sample}.bam \
					{output.bam}
				mv {output.bam}.bai {output.bai}
				;;
			*)
				sambamba markdup -t {threads} -l 6 \
					--tmpdir=/lscratch/$SLURM_JOB_ID/{wildcards.sample} \
					{input.bam} \
					{output.bam}
				mv {output.bam}.bai {output.bai}
				;;
		esac
		"""
#{config[biobambam2_version]}
#			mkdir -p sample_bam/markDupMetrics
# 			bammarkduplicates markthreads={threads} level=6 verbose=0 \
# 					M=sample_bam/markDupMetrics/{wildcards.sample}.markDup.methrics.tsv \
# 					tmpfile=/lscratch/$SLURM_JOB_ID/{wildcards.sample} \
# 					I=/lscratch/$SLURM_JOB_ID/{wildcards.sample}.bam \
# 					O={output.bam} \
# 					index=1 indexfilename={output.bai}
		#	merge:				--SORT_ORDER coordinate \
		#					--CREATE_INDEX true
#no pipe for picard MergeSamFiles; markduplicates requires indexed bam files.Default picard MarkDuplicates compression_level is 5
				# java -Xmx16g -XX:+UseG1GC -XX:ParallelGCThreads={threads} -jar $PICARD_JAR \
				# 	MarkDuplicates \
				# 	--INPUT {output.merged_bam} \
				# 	--OUTPUT {output.bam} \
				# 	--METRICS_FILE {output.metrics} \
				# 	--COMPRESSION_LEVEL 6 \
				# 	--CREATE_INDEX true \
				# 	--ASSUME_SORT_ORDER coordinate

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

localrules: multiqc_fastqc
rule multiqc_fastqc:
	input:
		expand('fastqc/{sample}', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		'fastqc/multiqc.done'
	shell:
		"""
		module load multiqc/1.14
		multiqc -f -o fastqc/multiqc_report_{config[analysis_batch_name]} fastqc/
		touch {output}
		"""
#output: directory('fastqc/multiqc_report')

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

localrules: picard_alignmentQC
rule picard_alignmentQC:
#insert size and alignment metrics
	input:
		bam = 'sample_bam/{sample}/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}/{sample}.markDup.bai'
	output:
		insert_size_metrics = 'picardQC/{sample}.insert_size_metrics.txt',
		insert_size_histogram = 'picardQC/{sample}.insert_size_histogram.pdf',
		alignment_metrics = 'picardQC/{sample}.alignment_metrics.txt'
	threads: 4
	shell:
		"""
		if [[ $(module list 2>&1 | grep "R/" | wc -l) < 1 ]]; then module load {config[R_version]}; fi
		if [[ $(module list 2>&1 | grep "picard" | wc -l) < 1 ]]; then module load {config[picard_version]}; fi
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

#sbatch 8 threads and 32g
localrules: coverage
rule coverage:
	input:
		bam = 'sample_bam/{sample}/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}/{sample}.markDup.bai'
	output:
		thresholds = 'coverage/mosdepth/{sample}.md.thresholds.bed.gz',
		summary = 'coverage/mosdepth/{sample}.md.mosdepth.summary.txt',
		region_summary = temp('coverage/mosdepth/{sample}.region.summary.tsv'),
		xlsx = 'coverage/{sample}.coverage.xlsx'
	threads: 8
	shell:
		"""
		if [[ $(module list 2>&1 | grep "mosdepth" | wc -l) < 1 ]]; then module load {config[mosdepth_version]}; fi
		if [[ $(module list 2>&1 | grep "R/" | wc -l) < 1 ]]; then module load {config[R_version]}; fi
		cd coverage/mosdepth
		mosdepth -t {threads} --no-per-base --by {config[bed]} --mapq 0 --fast-mode --thresholds 10,20,30 \
			{wildcards.sample}.md ../../{input.bam}
		cd ../..
		#mv {wildcards.sample}.md.* coverage/mosdepth/.
		zcat {output.thresholds} \
			 | sed '1 s/^.*$/chr\tstart\tend\tgene\tcoverageTen\tcoverageTwenty\tcoverageThirty/' \
			 > {output.thresholds}.tsv
		echo -e "sample\tlength\tmean" > {output.region_summary}
		tail -n 1 {output.summary} | cut -f 2,4 | sed 's/^/{wildcards.sample}\t/' >> {output.region_summary}
		Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/mosdepth_bed_coverage.R \
			{output.thresholds}.tsv {config[OGL_Dx_research_genes]} {output.region_summary} {output.xlsx}
		rm {output.thresholds}.tsv
		"""
#mv {wildcards.sample}.per-base.bed.gz* coverage/mosdepth/.

localrules: mean_coverage
rule mean_coverage:
	input:
		expand('coverage/mosdepth/{sample}.md.mosdepth.summary.txt', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		'coverage/mean.coverage.done.txt'
	shell:
		"""
		echo -e "sample\tlength\tmean\ton_target_rate" > coverage/{config[analysis_batch_name]}.mean.coverage.summary.tsv
		for file in {input}; do
			filename=$(basename $file)
			sm=$(echo $filename | sed 's/.md.mosdepth.summary.txt//')
			on_target=$(tail -n 2 $file | awk 'NR>1{{print $3/p}} {{p=$3}}')
			tail -n 1 $file | awk -v on_target="$on_target" -v sm="$sm" -F"\t" 'BEGIN{{OFS="\t"}}{{print sm,$2,$4,on_target}}' - >> coverage/{config[analysis_batch_name]}.mean.coverage.summary.tsv
 		done
		touch {output}
		"""

# 30% smaller!
localrules: bam_to_cram
rule bam_to_cram:
	input:
		bam = 'sample_bam/{sample}/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}/{sample}.markDup.bai'
	output:
		cram = 'bam/{sample}.cram',
		crai = 'bam/{sample}.crai'
	threads:
		4
	shell:
		"""
		module load {config[samtools_version]}
		samtools view -T {config[ref_genome]} --threads {threads} --output-fmt cram,store_md=1,store_nm=1 -o {output.cram} {input.bam}
		samtools index -@ {threads} {output.cram} {output.crai}
		"""

# samtools view -O cram,store_md=1,store_nm=1 -o aln.cram aln.bam
# samtools view --input-fmt cram,decode_md=0 -o aln.new.bam aln.cram
#old: samtools sort -O bam -l 0 --threads {threads} -T /lscratch/$SLURM_JOB_ID {input.bam} | \
#		samtools view -T {config[ref_genome]} --threads {threads} -C -o {output.cram} -
# or # samtools view -T {config[ref_genome]} --threads {threads} -C -o {output.cram} {input.bam}


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
				perl /data/OGL/resources/git/CoNVaDING/CoNVaDING.pl -mode StartWithBam \
					-inputDir sample_bam/{wildcards.sample} \
					-outputDir /lscratch/$SLURM_JOB_ID \
					-bed {config[bed]} \
					-useSampleAsControl \
					-controlsDir {config[CoNVaDING_ctr_dir]}_male \
					-rmDup
				mkdir -p CoNVaDING/normalized_coverage_male
				cp /lscratch/$SLURM_JOB_ID/{wildcards.sample}.markDup.aligned.only.normalized.coverage.txt \
					CoNVaDING/normalized_coverage_male/{wildcards.sample}.markDup.aligned.only.normalized.coverage.txt
				chgrp OGL {config[CoNVaDING_ctr_dir]}_male/{wildcards.sample}.markDup.aligned.only.normalized.coverage.txt
				Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/chrRD.R \
					CoNVaDING/normalized_coverage_male/{wildcards.sample}.markDup.aligned.only.normalized.coverage.txt \
					CoNVaDING/normalized_coverage_male/{wildcards.sample}.chrRD.pdf \
					{config[chrRD_highcutoff]} \
					{config[chrRD_lowcutoff]} \
					CoNVaDING/normalized_coverage_male/{wildcards.sample}.abnormalChr.tsv \
					1
				touch {output}
				;;
			"2")
				perl /data/OGL/resources/git/CoNVaDING/CoNVaDING.pl -mode StartWithBam \
					-inputDir sample_bam/{wildcards.sample} \
					-outputDir /lscratch/$SLURM_JOB_ID \
					-bed {config[bed]} \
					-useSampleAsControl \
					-controlsDir {config[CoNVaDING_ctr_dir]}_female \
					-rmDup
				mkdir -p CoNVaDING/normalized_coverage_female
				cp /lscratch/$SLURM_JOB_ID/{wildcards.sample}.markDup.aligned.only.normalized.coverage.txt \
					CoNVaDING/normalized_coverage_female/{wildcards.sample}.markDup.aligned.only.normalized.coverage.txt
				chgrp OGL {config[CoNVaDING_ctr_dir]}_female/{wildcards.sample}.markDup.aligned.only.normalized.coverage.txt
				Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/chrRD.R \
					CoNVaDING/normalized_coverage_female/{wildcards.sample}.markDup.aligned.only.normalized.coverage.txt \
					CoNVaDING/normalized_coverage_female/{wildcards.sample}.chrRD.pdf \
					{config[chrRD_highcutoff]} \
					{config[chrRD_lowcutoff]} \
					CoNVaDING/normalized_coverage_female/{wildcards.sample}.abnormalChr.tsv \
					2
				touch {output}
				;;
			*)
				perl /data/OGL/resources/git/CoNVaDING/CoNVaDING.pl -mode StartWithBam \
					-inputDir sample_bam/{wildcards.sample} \
					-outputDir /lscratch/$SLURM_JOB_ID \
					-bed {config[bed]} \
					-useSampleAsControl \
					-controlsDir {config[CoNVaDING_ctr_dir]} \
					-rmDup
				mkdir -p CoNVaDING/normalized_coverage
				cp /lscratch/$SLURM_JOB_ID/{wildcards.sample}.markDup.aligned.only.normalized.coverage.txt \
					CoNVaDING/normalized_coverage/{wildcards.sample}.markDup.aligned.only.normalized.coverage.txt
				chgrp OGL {config[CoNVaDING_ctr_dir]}/{wildcards.sample}.markDup.aligned.only.normalized.coverage.txt
				Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/chrRD.R \
					CoNVaDING/normalized_coverage/{wildcards.sample}.markDup.aligned.only.normalized.coverage.txt \
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

#localrules: CoNVaDING_2 #10min 1 cpu and <0.5g mem
rule CoNVaDING_2:
	input:
		expand('CoNVaDING/progress1.{sample}', sample=list(SAMPLE_LANEFILE.keys())),
	output:
		'CoNVaDING/progress2.done'
	shell:
		"""
 		filetest0=$((ls CoNVaDING/normalized_coverage/*.markDup.aligned.only.normalized.coverage.txt >> /dev/null 2>&1 && echo TRUE) || echo FALSE)
		if [ $filetest0 == "TRUE" ];
		then
			#cp -a CoNVaDING/normalized_coverage/*.markDup.aligned.only.normalized.coverage.txt {config[CoNVaDING_ctr_dir]}
			perl /data/OGL/resources/git/CoNVaDING/CoNVaDING.pl -mode StartWithMatchScore \
				-inputDir CoNVaDING/normalized_coverage \
				-outputDir  CoNVaDING/MatchScore \
				-controlsDir {config[CoNVaDING_ctr_dir]}
			perl /data/OGL/resources/git/CoNVaDING/CoNVaDING.pl \
  				-mode StartWithBestScore \
  				-inputDir CoNVaDING/MatchScore \
  				-outputDir CoNVaDING/CNV_hiSens \
  				-controlsDir {config[CoNVaDING_ctr_dir]} \
  				-ratioCutOffLow 0.71 \
  				-ratioCutOffHigh 1.35
			# perl /data/OGL/resources/git/CoNVaDING/CoNVaDING.pl \
  			# 	-mode GenerateTargetQcList \
  			# 	-inputDir {config[CoNVaDING_ctr_dir]} \
  			# 	-outputDir CoNVaDING/TargetQcList \
  			# 	-controlsDir {config[CoNVaDING_ctr_dir]} \
  			# 	-ratioCutOffLow 0.71 \
  			# 	-ratioCutOffHigh 1.35
			# perl /data/OGL/resources/git/CoNVaDING/CoNVaDING.pl \
			# 	-mode CreateFinalList \
			# 	-inputDir CoNVaDING/CNV_hiSens \
  			# 	-targetQcList CoNVaDING/TargetQcList \
  			# 	--outputDir CoNVaDING/finalList
		fi
		filetest1=$((ls CoNVaDING/normalized_coverage_male/*.markDup.aligned.only.normalized.coverage.txt >> /dev/null 2>&1 && echo TRUE) || echo FALSE)
		if [ $filetest1 == "TRUE" ];
		then
			#cp -a CoNVaDING/normalized_coverage_male/*.markDup.aligned.only.normalized.coverage.txt {config[CoNVaDING_ctr_dir]}
			perl /data/OGL/resources/git/CoNVaDING/CoNVaDING.pl -mode StartWithMatchScore \
				-inputDir CoNVaDING/normalized_coverage_male \
				-outputDir  CoNVaDING/MatchScore_male \
				-controlsDir {config[CoNVaDING_ctr_dir]}_male \
				-sexChr
			perl /data/OGL/resources/git/CoNVaDING/CoNVaDING.pl \
  				-mode StartWithBestScore \
  				-inputDir CoNVaDING/MatchScore_male \
  				-outputDir CoNVaDING/CNV_hiSens \
  				-controlsDir {config[CoNVaDING_ctr_dir]}_male \
  				-ratioCutOffLow 0.71 \
  				-ratioCutOffHigh 1.35 \
				-sexChr
		fi
		filetest2=$((ls CoNVaDING/normalized_coverage_female/*.markDup.aligned.only.normalized.coverage.txt >> /dev/null 2>&1 && echo TRUE) || echo FALSE)
		if [ $filetest2 == "TRUE" ];
		then
			#cp -a CoNVaDING/normalized_coverage_female/*.markDup.aligned.only.normalized.coverage.txt {config[CoNVaDING_ctr_dir]}
			perl /data/OGL/resources/git/CoNVaDING/CoNVaDING.pl -mode StartWithMatchScore \
				-inputDir CoNVaDING/normalized_coverage_female \
				-outputDir  CoNVaDING/MatchScore_female \
				-controlsDir {config[CoNVaDING_ctr_dir]}_female \
				-sexChr
			perl /data/OGL/resources/git/CoNVaDING/CoNVaDING.pl \
  				-mode StartWithBestScore \
  				-inputDir CoNVaDING/MatchScore_female \
  				-outputDir CoNVaDING/CNV_hiSens \
  				-controlsDir {config[CoNVaDING_ctr_dir]}_female \
  				-ratioCutOffLow 0.71 \
  				-ratioCutOffHigh 1.35 \
				-sexChr
		fi
		for i in CoNVaDING/CNV_hiSens/*.shortlist.txt; do awk -F "\t" '{{print FILENAME"\t"$0}}' $i >> CoNVaDING/shortlist.temp; done
		awk -F"\t" 'BEGIN{{OFS="\t"}} {{sub(/CoNVaDING\/CNV_hiSens\//,""); sub(/.markDup.aligned.only.best.score.shortlist.txt/,""); print }}' CoNVaDING/shortlist.temp \
			| grep -v -P 'CHR\tSTART' - > CoNVaDING/SHORTlist.txt && \
			echo -e "SAMPLE\tCHR\tSTART\tSTOP\tGENE\tNUMBER_OF_TARGETS\tNUMBER_OF_TARGETS_PASS_SHAPIRO-WILK_TEST\tABBERATION" \
			| cat - CoNVaDING/SHORTlist.txt > CoNVaDING/tmpout && mv CoNVaDING/tmpout CoNVaDING/{config[analysis_batch_name]}.SHORTlist.txt
		rm CoNVaDING/shortlist.temp
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
		gvcf = 'deepvariant/gvcf/{sample}.dv.g.vcf.gz',
		filteredvcf = temp('deepvariant/vcf/{sample}.dv.filtered.vcf.gz'),
		filteretbi = temp('deepvariant/vcf/{sample}.dv.filtered.vcf.gz.tbi'),
		phasedvcf = 'deepvariant/vcf/{sample}.dv.phased.vcf.gz',
		phasedtbi = 'deepvariant/vcf/{sample}.dv.phased.vcf.gz.tbi'
	threads: 16
	shell:
		"""
		module load {config[deepvariant_version]}
		PROJECT_WD="$PWD"
		N_SHARDS="4"
		mkdir -p /lscratch/$SLURM_JOB_ID/{wildcards.sample}
		WORK_DIR=/lscratch/$SLURM_JOB_ID/{wildcards.sample}
		cp {input} $WORK_DIR
		cd $WORK_DIR
		run_deepvariant --model_type WES --num_shards $N_SHARDS \
			--ref {config[ref_genome]} \
			--regions {config[padded_bed]} \
			--reads $WORK_DIR/$(basename {input.bam})  \
			--output_vcf $WORK_DIR/$(basename {output.vcf}) \
			--output_gvcf $WORK_DIR/$(basename {output.gvcf}) \
			--sample_name {wildcards.sample} \
			--intermediate_results_dir $WORK_DIR
		cd $PROJECT_WD
		cp $WORK_DIR/$(basename {output.vcf})* deepvariant/vcf
		cp $WORK_DIR/$(basename {output.gvcf})* deepvariant/gvcf
		module unload {config[deepvariant_version]}
		module load {config[samtools_version]}
		bcftools norm --multiallelics -any --output-type u {output.vcf} \
			| bcftools norm -d exact --output-type u - \
			| bcftools filter --include 'FILTER="PASS" & FORMAT/AD[0:1]>2' --output-type z --output {output.filteredvcf}
		sleep 2
		tabix -f -p vcf {output.filteredvcf}
		module load {config[whatshap_version]}
		whatshap phase --reference {config[ref_genome]} --indels {output.filteredvcf} {input.bam} | bgzip -f > {output.phasedvcf}
		tabix -f -p vcf {output.phasedvcf}
		"""

#v1.1.0 added model.bin model.xml model.mapping to the working directory, causing problems for multiple samples
#instead of seperating to three steps, can also go to lscratch work directory, but add $PWD to files below.
#v1.0.0 added model files to temp folder.
# ( time seq 0 $((N_SHARDS-1)) | parallel -q --halt 2 --line-buffer \
# 	make_examples --mode calling \
# 	--ref {config[ref_genome]} \
# 	--reads {input.bam} \
# 	--examples $WORK_DIR/make_examples.tfrecord@$N_SHARDS.gz \
# 	--gvcf $WORK_DIR/gvcf.tfrecord@$N_SHARDS.gz \
# 	--regions {config[padded_bed]} \
# 	--sample_name {wildcards.sample} --task {{}} )
# cd $WORK_DIR
# ( time call_variants --outfile $WORK_DIR/call_variants_output.tfrecord.gz \
# 	--examples $WORK_DIR/make_examples.tfrecord@$N_SHARDS.gz \
# 	--checkpoint "/opt/models/wes/model.ckpt" --use_openvino )
# cd $PROJECT_WD
# ( time postprocess_variants \
# 	--ref {config[ref_genome]} \
# 	--infile $WORK_DIR/call_variants_output.tfrecord.gz \
# 	--outfile {output.vcf} \
# 	--nonvariant_site_tfrecord_path $WORK_DIR/gvcf.tfrecord@$N_SHARDS.gz \
# 	--gvcf_outfile {output.gvcf} \
# 	--sample_name {wildcards.sample} )
#--sample_name {wildcards.sample}; If default (not specified) then name from RG
# QUAL>0 replaced with Pass filter 1/29/2021


localrules: merge_deepvariant_vcf
rule merge_deepvariant_vcf:
	input:
		vcf = expand('deepvariant/vcf/{sample}.dv.phased.vcf.gz', sample=list(SAMPLE_LANEFILE.keys())),
		tbi = expand('deepvariant/vcf/{sample}.dv.phased.vcf.gz.tbi', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		'deepvariant/deepvariantVcf.merge.done.txt'
	threads: 8
	shell:
		"""
		if [[ $(module list 2>&1 | grep "samtools" | wc -l) < 1 ]]; then module load {config[samtools_version]}; fi
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

rule glnexus:
	input:
		vcf = expand('deepvariant/gvcf/{sample}.dv.g.vcf.gz', sample=list(SAMPLE_LANEFILE.keys())),
		bam = expand('sample_bam/{sample}/{sample}.markDup.bam', sample=list(SAMPLE_LANEFILE.keys())),
		bai = expand('sample_bam/{sample}/{sample}.markDup.bai', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		'deepvariant/deepvariant.gvcf.merge.done.txt'
	threads: 24
	shell:
		"""
		module load {config[glnexus_version]} {config[samtools_version]} {config[whatshap_version]} parallel
		WORK_DIR="/lscratch/${{SLURM_JOB_ID}}"
		glnexus_cli --dir /lscratch/$SLURM_JOB_ID/glnexus --config DeepVariantWES --bed {config[padded_bed]} \
			--threads {threads} \
			{input.vcf} \
			| bcftools norm --multiallelics -any --output-type u --no-version \
			| bcftools norm --check-ref s --fasta-ref {config[ref_genome]} --output-type u --no-version - \
			| bcftools +fill-tags - -Ou -- -t AC,AC_Hom,AC_Het,AN,AF \
			| bcftools annotate --threads {threads} --set-id 'dvg_%CHROM\:%POS%REF\>%ALT' --no-version - -Oz -o $WORK_DIR/glnexus.vcf.gz
		tabix -f -p vcf $WORK_DIR/glnexus.vcf.gz
		head -n 19 /data/OGL/resources/whatshap/vcf.contig.filename.{config[genomeBuild]}.txt > $WORK_DIR/contig.txt
		CONTIGFILE="$WORK_DIR/contig.txt"
		mkdir -p /lscratch/$SLURM_JOB_ID/chr
		mkdir -p /lscratch/$SLURM_JOB_ID/phased
		cat $CONTIGFILE | parallel -C "\t" -j 19 "bcftools filter -r {{1}} --output-type z $WORK_DIR/glnexus.vcf.gz -o $WORK_DIR/chr/{{2}}.vcf.gz"
		cat $CONTIGFILE | parallel -C "\t" -j 19 "tabix -f -p vcf $WORK_DIR/chr/{{2}}.vcf.gz"
		cat $CONTIGFILE | parallel -C "\t" -j 19 --tmpdir $WORK_DIR --eta --halt 2 --line-buffer \
		 	--tag "whatshap phase --reference {config[ref_genome]} \
			--indels $WORK_DIR/chr/{{2}}.vcf.gz {input.bam} \
			| bgzip -f > $WORK_DIR/phased/{{2}}.phased.vcf.gz"
		cat $CONTIGFILE | parallel -C "\t" -j 19 "tabix -f -p vcf $WORK_DIR/phased/{{2}}.phased.vcf.gz"
		PHASEDCHRFILE=""
		cut -f 2 $CONTIGFILE > $WORK_DIR/temp.chr.txt
		while read line; do PHASEDCHRFILE+=" /lscratch/${{SLURM_JOB_ID}}/phased/$line.phased.vcf.gz"; done < $WORK_DIR/temp.chr.txt
		echo "chr files are $PHASEDCHRFILE"
		bcftools concat --threads {threads} --output-type z $PHASEDCHRFILE > deepvariant/{config[analysis_batch_name]}.glnexus.phased.vcf.gz
		tabix -f -p vcf deepvariant/{config[analysis_batch_name]}.glnexus.phased.vcf.gz
		touch {output}
		"""

#AC_Hemi removed because samtools 1.13 set it as 0 for chrX 12/23/21
#try phasing gvcf files after providing ped files
#get all gvcf files together and use glnexus to combine samples.
#configuration of DeepVariantWES removed 139 orf15 variant.

#freebayes avoids indel_realign, base-quality_recalibration, better to run on all of the bam files
#freebayes can identify MNP and complex indels.
#need to have unique RG ID fields for each sample. Use the RG made by the wrapper.
#Using CleanSam and FixMateInformation, freebayes was able to call one additional ins variant using the panel NA12878 data.

rule freebayes_phasing:
	input:
		bam = 'sample_bam/{sample}/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}/{sample}.markDup.bai'
	output:
		vcf = 'freebayes/vcf/{sample}.vcf.gz',
		filteredvcf = temp('freebayes/vcf/{sample}.filtered.vcf.gz'),
		tbi = temp('freebayes/vcf/{sample}.filtered.vcf.gz.tbi'),
		phasedvcf = 'freebayes/vcf/{sample}.phased.vcf.gz',
		phasedvcf_tbi = 'freebayes/vcf/{sample}.phased.vcf.gz.tbi'
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
		bcftools filter --regions-file {config[padded_bed]} --output-type v {output.vcf} \
			| bcftools norm --multiallelics -any --output-type v - \
			| vt decompose_blocksub -p -m -d 2 - \
			| bcftools norm --check-ref s --fasta-ref {config[ref_genome]} --output-type v - \
			| bcftools norm -d exact --output-type v - \
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
#added --regions-file 1/29/2021
# bcftools norm -d exact means remove duplicates if they are identical (keep the first instance)
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


localrules: merge_freebayes
rule merge_freebayes:
	input:
		vcf = expand('freebayes/vcf/{sample}.phased.vcf.gz', sample=list(SAMPLE_LANEFILE.keys())),
		tbi = expand('freebayes/vcf/{sample}.phased.vcf.gz.tbi', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		'freebayes/freebayes.merge.done.txt'
	threads: 8
	shell:
		"""
		if [[ $(module list 2>&1 | grep "samtools" | wc -l) < 1 ]]; then module load {config[samtools_version]}; fi
		case "{input.vcf}" in
			*\ *)
				bcftools merge --merge none --missing-to-ref --output-type z --threads {threads} {input.vcf} \
				> freebayes/{config[analysis_batch_name]}.freebayes.vcf.gz
				sleep 2
				tabix -f -p vcf freebayes/{config[analysis_batch_name]}.freebayes.vcf.gz
				;;
			*)
				cp -p -l {input.vcf} freebayes/{config[analysis_batch_name]}.freebayes.vcf.gz
				cp -p -l {input.tbi} freebayes/{config[analysis_batch_name]}.freebayes.vcf.gz.tbi
				;;
		esac
		touch {output}
		"""
#bcftools merge: freebayes INFO/QA minimum value among all samples is retained in the INFO.

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

localrules: merge_dv_fb_vcfs
rule merge_dv_fb_vcfs:
	input:
		'deepvariant/deepvariantVcf.merge.done.txt',
		'deepvariant/deepvariant.gvcf.merge.done.txt',
		'freebayes/freebayes.merge.done.txt'
	output:
		'prioritization/dv_fb.merge.done.txt'
	threads: 8
	shell:
		"""
		if [[ $(module list 2>&1 | grep "samtools" | wc -l) < 1 ]]; then module load {config[samtools_version]}; fi
		WORK_DIR=/lscratch/$SLURM_JOB_ID
		bcftools isec -p $WORK_DIR/dv -w 2 --collapse none --output-type u --threads {threads} \
			deepvariant/{config[analysis_batch_name]}.glnexus.phased.vcf.gz \
			deepvariant/{config[analysis_batch_name]}.dv.phased.vcf.gz
		bcftools +fill-tags $WORK_DIR/dv/0001.bcf -Ov -- -t AC,AC_Hom,AC_Het,AN,AF \
			| sed 's#0/0:.:.:.#0/0:10:10:10,0#g' - \
			| bcftools annotate --threads {threads} --set-id 'dv_%CHROM\:%POS%REF\>%ALT' --no-version - -Oz -o $WORK_DIR/dv/dv.hf.vcf.gz
		tabix -f -p vcf $WORK_DIR/dv/dv.hf.vcf.gz
		bcftools concat --threads {threads} -a --rm-dups none --no-version \
			deepvariant/{config[analysis_batch_name]}.glnexus.phased.vcf.gz $WORK_DIR/dv/dv.hf.vcf.gz -Oz \
			-o $WORK_DIR/dv.glnexus.hf.vcf.gz
		tabix -f -p vcf $WORK_DIR/dv.glnexus.hf.vcf.gz
		rm -r -f $WORK_DIR/dv
		bcftools isec --threads {threads} -p $WORK_DIR --collapse none -Oz \
			$WORK_DIR/dv.glnexus.hf.vcf.gz \
			freebayes/{config[analysis_batch_name]}.freebayes.vcf.gz
		rm $WORK_DIR/0003.vcf.gz &
		#bcftools annotate --threads {threads} --set-id 'dv_%CHROM\:%POS%REF\>%ALT' \
		#	--no-version $WORK_DIR/0000.vcf -Oz -o $WORK_DIR/dv.vcf.gz
		#rm $WORK_DIR/0000.vcf &
		bcftools annotate --threads {threads} --set-id 'fb_%CHROM\:%POS%REF\>%ALT' -x ^INFO/QA,FORMAT/RO,FORMAT/QR,FORMAT/AO,FORMAT/QA,FORMAT/GL \
			--no-version $WORK_DIR/0001.vcf.gz -Ov - \
			| sed 's#0/0:.:.:.#0/0:10:10:10,0#g' - \
			| bcftools +fill-tags - -Oz -o $WORK_DIR/fb.vcf.gz -- -t AC,AC_Hom,AC_Het,AN,AF
		rm $WORK_DIR/0001.vcf.gz &
		zcat $WORK_DIR/0002.vcf.gz | sed 's/\tdv/\tfbDv/' | bgzip -@ {threads} > $WORK_DIR/dvFb.vcf.gz
		#bcftools annotate --threads {threads} --set-id 'dvFb_%CHROM\:%POS%REF\>%ALT' \
		#	--no-version $WORK_DIR/0002.vcf -Oz -o $WORK_DIR/dvFb.vcf.gz
		rm $WORK_DIR/0002.vcf.gz &
		tabix -f -p vcf $WORK_DIR/0000.vcf.gz
		tabix -f -p vcf $WORK_DIR/fb.vcf.gz
		tabix -f -p vcf $WORK_DIR/dvFb.vcf.gz
		bcftools concat --threads {threads} -a --rm-dups none --no-version \
			$WORK_DIR/dvFb.vcf.gz $WORK_DIR/0000.vcf.gz $WORK_DIR/fb.vcf.gz -Oz \
			-o prioritization/{config[analysis_batch_name]}.vcf.gz
		tabix -f -p vcf prioritization/{config[analysis_batch_name]}.vcf.gz
		touch {output}
		"""
###crossmap requires chrM MT, which needs to be super careful when using sed because it may change the INFO in the header.
# if [[ {config[genomeBuild]} == "GRCh38" ]]; then
# 	module load {config[crossmap_version]}
# 	hg19ref=/data/OGL/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta
# 	crossmap vcf /data/OGL/resources/ucsc/hg38ToHg19.over.chain.gz \
# 		prioritization/{config[analysis_batch_name]}.vcf.gz \
# 		$hg19ref \
# 		$WORK_DIR/GRCh37.vcf
# 	sed -e 's/^chr//' -e 's/<ID=chr/<ID=/' $WORK_DIR/GRCh37.vcf \
# 		| bcftools norm --check-ref s --fasta-ref $hg19ref --output-type u - \
# 		| bcftools sort -m 26G -T $WORK_DIR/samtools -Ou - \
# 		| bcftools norm --threads $(({threads}-4)) -d exact --output-type z - -o prioritization/{config[analysis_batch_name]}.GRCh37.vcf.gz
# 	tabix -f -p vcf prioritization/{config[analysis_batch_name]}.GRCh37.vcf.gz
# fi

#glnexus 0/0 genotypes could have DP GQ values less than 10 12/24/21
#bcftools +fill-tags AC_Hemi does not work in samtools/1.13
# bcftools concat --threads {threads} -a --rm-dups none --no-version \
# 	deepvariant/{config[analysis_batch_name]}.dv.phased.vcf.gz \
# 	freebayes/{config[analysis_batch_name]}.freebayes.vcf.gz \
# 	-Oz -o prioritization/{config[analysis_batch_name]}.vcf.gz

#localrules: scramble #>1 hour with deletion calling
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
		annovarR = temp('scramble_anno/{sample}.forR.txt'),
		anno = 'scramble_anno/{sample}.scramble.tsv',
		anno_xlsx = 'scramble_anno/{sample}.scramble.xlsx',
		del_anno = 'scramble_anno/{sample}.scramble.del.tsv'
	resources: res=1
	shell:
		"""
		if [[ $(module list 2>&1 | grep "R/" | wc -l) < 1 ]]; then module load {config[R_version]}; fi
		if [[ $(module list 2>&1 | grep "annovar" | wc -l) < 1 ]]; then module load {config[annovar_version]}; fi
		if [[ {config[genomeBuild]} == "GRCh38" ]]; then
			ver=hg38
		else
			ver=hg19
		fi
		if [[ $(wc -l {input.mei} | cut -d " " -f 1) == 1 ]]
		then
			touch {output.avinput}
			touch {output.annovarR}
			touch {output.anno}
			touch {output.anno_xlsx}
		else
			cut -f 1 {input.mei} | awk -F ":" 'BEGIN{{OFS="\t"}} NR>1 {{print $1,$2,$2,"0","-"}}' > {output.avinput}
			table_annovar.pl {output.avinput} \
				$ANNOVAR_DATA/$ver \
				-buildver $ver \
				-remove \
				-out scramble_anno/{wildcards.sample} \
				--protocol refGene \
				-operation g \
				--argument '-splicing 100 -hgvs' \
				--polish -nastring . \
				--thread 1
			awk -F"\t" 'BEGIN{{OFS="\t"}} NR==1 {{print "Func_refGene","Gene","Intronic","AA"}} NR>1 {{print $6,$7,$8,$10}}' scramble_anno/{wildcards.sample}."$ver"_multianno.txt | paste {input.mei} - > {output.annovarR}
			rm scramble_anno/{wildcards.sample}."$ver"_multianno.txt
			Rscript /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/scramble_anno.R {output.annovarR} {config[SCRAMBLEdb]} {config[OGL_Dx_research_genes]} {config[HGMDtranscript]} {wildcards.sample} {output.anno} {output.anno_xlsx}
		fi
		if [[ $(wc -l {input.deletion} | cut -d " " -f 1) == 1 ]]
		then
			touch {output.del_anno}
		else
			if [[ $(module list 2>&1 | grep "annotsv" | wc -l) < 1 ]]; then module load {config[annotsv_version]}; fi
			tail -n +2 {input.deletion} | awk -F"\t" 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,"DEL"}}' > {input.deletion}.bed
			AnnotSV -genomeBuild {config[genomeBuild]} -SVinputFile {input.deletion}.bed -SVinputInfo 0 -svtBEDcol 4 -outputFile {output.del_anno}.temp
			Rscript /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/scramble_del_edit.R {output.del_anno}.temp.tsv {config[scrambleDELdb]} {output.del_anno}
			rm {output.del_anno}.temp.tsv
		fi
		"""
#--intronhgvs 100
# if [[ {config[genomeBuild]} == "GRCh38" ]]; then
# 	crossmap region /data/OGL/resources/ucsc/hg38ToHg19.over.chain.gz {input.deletion}.bed {input.deletion}.GRCh37.bed
# 	mv {input.deletion}.GRCh37.bed {input.deletion}.bed
# fi

#localrules: manta
rule manta:
	input:
		bam = 'sample_bam/{sample}/{sample}.markDup.bam',
 		bai = 'sample_bam/{sample}/{sample}.markDup.bai'
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
			-SVinputInfo 1 -genomeBuild {config[genomeBuild]} \
			-outputDir $RUNDIR \
			-outputFile $RUNDIR/manta.{wildcards.sample}.annotated.tsv
		if [ -e $RUNDIR/manta.{wildcards.sample}.annotated.tsv ];
		then
			cp $RUNDIR/manta.{wildcards.sample}.annotated.tsv {output}
		else
			touch {output}
		fi
		"""

#AnnotSV can pick sample from a multi-sample vcf file by using ?? test to find out the sample operation
#if [[ {config[genomeBuild]} == "GRCh38" ]]; then
# 	module load {config[crossmap_version]}
# 	crossmap vcf /data/OGL/resources/ucsc/hg38ToHg19.over.chain.gz \
# 		$RUNDIR/results/variants/diploidSV.vcf.gz \
# 		/data/OGL/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta \
# 		$RUNDIR/results/variants/diploidSV.GRCh37.vcf.gz
# 	mv $RUNDIR/results/variants/diploidSV.GRCh37.vcf.gz $RUNDIR/results/variants/diploidSV.vcf.gz
# 	tabix -f -p vcf $RUNDIR/results/variants/diploidSV.vcf.gz
# fi

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

localrules: bcm_locus
rule bcm_locus:
	input:
		bam = 'sample_bam/{sample}/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}/{sample}.markDup.bai'
	output:
		vcf = temp('bcmlocus/{sample}.vcf'),
		avinput = temp('bcmlocus/{sample}.avinput'),
		bcm_out = 'bcmlocus/{sample}.bcmlocus.tsv'
	shell:
		"""
		if [[ $(module list 2>&1 | grep "mosdepth" | wc -l) < 1 ]]; then module load {config[mosdepth_version]}; fi
		if [[ $(module list 2>&1 | grep "R/" | wc -l) < 1 ]]; then module load {config[R_version]}; fi
		mkdir -p bcmlocus/mosdepth
		cd bcmlocus/mosdepth
		mosdepth -t {threads} --no-per-base --by {config[bcmlocus_bed]} --use-median --mapq 0 --fast-mode \
			{wildcards.sample}.md ../../{input.bam}
		cd ../..
		if [[ $(module list 2>&1 | grep "samtools" | wc -l) < 1 ]]; then module load {config[samtools_version]}; fi
		if [[ $(module list 2>&1 | grep "freebayes" | wc -l) < 1 ]]; then module load {config[freebayes_version]}; fi
		if [[ $(module list 2>&1 | grep "annovar" | wc -l) < 1 ]]; then module load {config[annovar_version]}; fi
		freebayes -f {config[ref_genome]} --max-complex-gap 80 -p 6 -C 3 -F 0.05 \
			--genotype-qualities --strict-vcf --use-mapping-quality \
			--targets /data/OGL/resources/bed/OPN1LWe2e5.bed \
			{input.bam} \
			| bcftools norm --multiallelics -any --output-type u - \
			| bcftools annotate --set-id '%CHROM\:%POS%REF\>%ALT' -x ^INFO/AF --output-type u --no-version \
			| bcftools norm --check-ref s --fasta-ref {config[ref_genome]} --output-type u --no-version - \
			| bcftools norm -d exact --output-type v - \
			> {output.vcf}
		convert2annovar.pl -format vcf4old {output.vcf} -includeinfo --outfile {output.avinput}
		if [[ {config[genomeBuild]} == "GRCh38" ]]; then
			ver=hg38
		else
			ver=hg19
		fi
		if [[ -s {output.avinput} ]]; then
			table_annovar.pl {output.avinput} \
				$ANNOVAR_DATA/$ver \
				-buildver $ver \
				-remove \
				-out {output.avinput} \
				--protocol refGeneWithVer \
				-operation g \
				--argument '-hgvs' \
				--polish -nastring . \
				--thread 1 \
				--otherinfo
			sed -i "1 s/Otherinfo1\tOtherinfo2\tOtherinfo3\tOtherinfo4\tOtherinfo5\tOtherinfo6\tOtherinfo7\tOtherinfo8\tOtherinfo9\tOtherinfo10/CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGT_FIELDS/" bcmlocus/{wildcards.sample}.avinput."$ver"_multianno.txt
			if [[ $(module list 2>&1 | grep "R/" | wc -l) < 1 ]]; then module load {config[R_version]}; fi
			Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/bcmlocus.R \
				/data/OGL/resources/bcmlocus.xlsx \
				{wildcards.sample} bcmlocus/{wildcards.sample}.avinput."$ver"_multianno.txt {output.bcm_out}
			rm bcmlocus/{wildcards.sample}.avinput."$ver"_multianno.txt
		else
			touch {output.bcm_out}
		fi
		"""
## Next step: get CN information

localrules: combine_bcmlocus
rule combine_bcmlocus:
	input:
		expand('bcmlocus/{sample}.bcmlocus.tsv', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		'bcmlocus/combine.bcmlocus.done.txt'
	shell:
		"""
		echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFunc.refGeneWithVer\tGene.refGeneWithVer\tGeneDetail.refGeneWithVer\tExonicFunc.refGeneWithVer\tAAChange.refGeneWithVer\tHGVSp\tAnnotation\tFunction\tACMG_Class\tNote\tSample\tINFO\tFORMAT\tGT_FIELDS" > bcmlocus/{config[analysis_batch_name]}.bcmlocus.all.tsv
		for file in {input}; do
			tail -n +2 $file >> bcmlocus/{config[analysis_batch_name]}.bcmlocus.all.tsv
 		done
		touch {output}
		"""

localrules: automap_roh
rule automap_roh:
	input:
		vcf = 'freebayes/vcf/{sample}.vcf.gz'
	output:
		tsv = temp('AutoMap/{sample}/{sample}.HomRegions.tsv'),
		bed = temp('AutoMap/{sample}/{sample}.HomRegions.bed'),
		annotated = 'AutoMap/{sample}/{sample}.HomRegions.annot.tsv'
	shell:
		"""
		if [[ $(module list 2>&1 | grep "samtools" | wc -l) < 1 ]]; then module load {config[samtools_version]}; fi
		if [[ $(module list 2>&1 | grep "bedtools" | wc -l) < 1 ]]; then module load {config[bedtools_version]}; fi
		if [[ $(module list 2>&1 | grep "R/" | wc -l) < 1 ]]; then module load {config[R_version]}; fi
		if [[ {config[genomeBuild]} == "GRCh38" ]]; then
			ver=hg38
		else
			ver=hg19
		fi
		mkdir -p /lscratch/$SLURM_JOB_ID/AutoMap
		zcat {input.vcf} > /lscratch/$SLURM_JOB_ID/AutoMap/{wildcards.sample}.vcf
		bash /data/OGL/resources/git/AutoMap/AutoMap_v1.2.sh \
			--vcf /lscratch/$SLURM_JOB_ID/AutoMap/{wildcards.sample}.vcf \
			--out AutoMap --genome $ver --chrX \
			--minsize 0.5 --minvar 15
		echo "AutoMap1.2 done"
		if [[ $(grep -v ^# {output.tsv} | wc -l) == 0 ]]; then
			touch {output.bed}
			touch {output.annotated}
			echo "no ROH region detected."
		else
			grep -v ^# {output.tsv} | cut -f 1-3 > {output.bed}
			module load {config[annotsv_version]}
			AnnotSV -SVinputFile {output.bed} \
				-SVinputInfo 1 -genomeBuild {config[genomeBuild]} \
				-outputDir AutoMap/{wildcards.sample} \
				-outputFile AutoMap/{wildcards.sample}/{wildcards.sample}.annotated.tsv
			Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/automap.R {output.tsv} AutoMap/{wildcards.sample}/{wildcards.sample}.annotated.tsv {config[OGL_Dx_research_genes]} {output.annotated}
			rm AutoMap/{wildcards.sample}/{wildcards.sample}.annotated.tsv
		fi
		"""

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
