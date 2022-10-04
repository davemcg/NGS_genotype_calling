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
metadata.close()

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
		rg_out = str(LANEFILE_READGROUP[lane_file + config['lane_pair_delim'][0] + '.gz'][0])
		return(rg_out)

# import CREST hg19 regions, chr1 to chrY
# REGIONS_file = config['regions']
# if '/home/$USER' in REGIONS_file:
# 	REGIONS_file = os.environ['HOME'] + REGIONS_file.split('$USER')[-1]
# REGIONS = open(REGIONS_file).readlines()
# REGIONS = [r.strip() for r in REGIONS]

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
	chr = '|'.join(CHRS),
	lane = '|'.join(list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in [y for sub in list(SAMPLE_LANEFILE.values()) for y in sub]])))
#	region = '|'.join(REGIONS)

rule all:
	input:
		# expand('CRESTanno/{sample}.predSV.xlsx', sample=list(SAMPLE_LANEFILE.keys())) if config['CREST'] == 'TRUE' else 'dummy.txt',
		expand('gvcfs/{sample}.g.vcf.gz', sample=list(SAMPLE_LANEFILE.keys())) if config['GATKgvcf'] == 'TRUE' else 'dummy.txt',
		# expand('recal_bam/{sample}.recal.bam', sample=list(SAMPLE_LANEFILE.keys())) if config['recal_bam'] == 'TRUE' else 'dummy.txt',
		expand('bam/{sample}.cram', sample=list(SAMPLE_LANEFILE.keys())) if config['cram'] == 'TRUE' else expand('bam/{sample}.bam', sample=list(SAMPLE_LANEFILE.keys())),
		# 'GATK_metrics/multiqc_report' if config['multiqc'] == 'TRUE' else 'dummy.txt',
		'fastqc/multiqc_report' if config['multiqc'] == 'TRUE' else 'dummy.txt',
		# expand('picardQC/{sample}.insert_size_metrics.txt', sample=list(SAMPLE_LANEFILE.keys())) if config['picardQC'] == 'TRUE' else 'dummy.txt',
		#'deepvariant/deepvariantVcf.merge.done.txt' if config['deepvariant'] == 'TRUE' else 'dummy.txt',
		'prioritization/dv_fb.merge.done.txt' if config['freebayes_phasing'] == 'TRUE' else 'dummy.txt',
		'coverage/mean.coverage.done.txt' if config['coverage'] == 'TRUE' else 'dummy.txt',
		expand('manta/manta.{sample}.annotated.tsv', sample=list(SAMPLE_LANEFILE.keys())),
		expand('scramble_anno/{sample}.scramble.tsv', sample=list(SAMPLE_LANEFILE.keys())) if config['SCRAMble'] == 'TRUE' else 'dummy.txt',
		expand('AutoMap/{sample}/{sample}.HomRegions.annot.tsv', sample=list(SAMPLE_LANEFILE.keys())),
		'bcmlocus/combine.bcmlocus.done.txt'


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
			expand('fastq/{{lane}}{pair}.gz', pair = config['lane_pair_delim'])
		output:
			bam = temp('lane_bam/{lane}.bam'),
			bai = temp('lane_bam/{lane}.bam.bai')
		params:
			read_group = rg
		threads: 32
		shell:
			"""
			export TMPDIR=/lscratch/$SLURM_JOB_ID
			echo {params.read_group}
			module load {config[bwa-mem2_version]} {config[samblaster_version]} {config[sambamba_version]}
			bwa-mem2 mem -t $(({threads}-2)) -K 100000000 -M -Y -B 4 -O 6 -E 1 -R {params.read_group} \
				{config[bwa-mem2_ref]} {input} \
			 	| samblaster -M --addMateTags --quiet \
				| sambamba sort -u --tmpdir=/lscratch/$SLURM_JOB_ID -t $(({threads}-2)) -o {output.bam} \
					<(sambamba view -S -f bam -l 0 -t $(({threads}-2)) /dev/stdin)
			"""
	localrules: cp_lane_bam
	rule cp_lane_bam:
		input:
			bam = lambda wildcards: expand('lane_bam/{lane}.bam', lane = list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in SAMPLE_LANEFILE[wildcards.sample]]))),
			bai = lambda wildcards: expand('lane_bam/{lane}.bam.bai', lane = list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in SAMPLE_LANEFILE[wildcards.sample]])))
		output:
			bam = 'sample_bam/{sample}.markDup.bam',
			bai = 'sample_bam/{sample}.markDup.bai'
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
			bam = 'sample_bam/{sample}.markDup.bam',
			bai = 'sample_bam/{sample}.markDup.bai'
		threads: 32
		params:
			read_group = rg
		shell:
			"""
			export TMPDIR=/lscratch/$SLURM_JOB_ID
			module load {config[bazam_version]} {config[bwa-mem2_version]} {config[samblaster_version]} {config[sambamba_version]}
			BAMFILE={input}
			if [ -e {input}.bai ] || [ -e ${{BAMFILE%.bam}}.bai ] || [ -e {input}.crai ] || [ -e ${{BAMFILE%.cram}}.crai ] ; then
				echo "index present"
			else
				samtools index -@ $(({threads}-2)) {input}
			fi
			case "{input}" in
				*bam)
					java -Xmx16g -jar $BAZAMPATH/bazam.jar -bam {input} \
					| bwa-mem2 mem -t $(({threads}/2)) -K 100000000 -M -Y -B 4 -O 6 -E 1 -p -R {params.read_group} {config[bwa-mem2_ref]} - \
			 		| samblaster -M --addMateTags --quiet \
					| sambamba sort -u --tmpdir=/lscratch/$SLURM_JOB_ID -t $(({threads}/2)) -o {output.bam} \
						<(sambamba view -S -f bam -l 0 -t $(({threads}/2)) /dev/stdin)
					mv {output.bam}.bai {output.bai}
					;;
				*cram)
					java -Xmx16g -Dsamjdk.reference_fasta={config[old_cram_ref]} -jar $BAZAMPATH/bazam.jar -bam {input} \
					| bwa-mem2 mem -t $(({threads}/2)) -K 100000000 -M -Y -B 4 -O 6 -E 1 -p -R {params.read_group} {config[bwa-mem2_ref]} - \
			 		| samblaster -M --addMateTags --quiet \
					| sambamba sort -u --tmpdir=/lscratch/$SLURM_JOB_ID -t $(({threads}/2)) -o {output.bam} \
						<(sambamba view -S -f bam -l 0 -t $(({threads}/2)) /dev/stdin)
					mv {output.bam}.bai {output.bai}
					;;
			esac
			"""
else:
	rule align: # < 10 min, 20g mem, 32 threads for bp exome.
		input:
			# config['lane_pair_delim'] is the string differentiating
			# the forward from reverse
			# e.g. ['_R1_001', '_R2_001'] if the file names are
			# sample17_R1_001.fastq.gz and sample17_R2_001.fastq.gz
			# for a set of paired end fastq
			# if you don't have a paired fastq set, give as ['']
			expand('fastq/{{lane}}{pair}.gz', pair = config['lane_pair_delim'])
		output:
			bam = temp('lane_bam/{lane}.bam'),
			bai = temp('lane_bam/{lane}.bam.bai')
		params:
			read_group = rg
		threads: 32
		shell:
			"""
			export TMPDIR=/lscratch/$SLURM_JOB_ID
			echo {params.read_group}
			module load {config[bwa-mem2_version]} {config[samblaster_version]} {config[sambamba_version]}
			bwa-mem2 mem -t $(({threads}-2)) -K 100000000 -M -Y -B 4 -O 6 -E 1 -R {params.read_group} \
				{config[bwa-mem2_ref]} {input} \
			 	| samblaster -M --acceptDupMarks --addMateTags --quiet \
				| sambamba sort -u --tmpdir=/lscratch/$SLURM_JOB_ID -t $(({threads}-2)) -o {output.bam} \
					<(sambamba view -S -f bam -l 0 -t $(({threads}-2)) /dev/stdin)
			"""
	rule merge_lane_bam:
		input:
			bam = lambda wildcards: expand('lane_bam/{lane}.bam', lane = list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in SAMPLE_LANEFILE[wildcards.sample]]))),
			bai = lambda wildcards: expand('lane_bam/{lane}.bam.bai', lane = list(set([re.split(r'|'.join(config['lane_pair_delim']),x.split('/')[-1])[0] for x in SAMPLE_LANEFILE[wildcards.sample]])))
		output:
			bam = 'sample_bam/{sample}.markDup.bam',
			bai = 'sample_bam/{sample}.markDup.bai'
		threads: 16
		shell:
			"""
			module load {config[sambamba_version]} {config[biobambam2_version]}
			mkdir -p sample_bam/markDupMetrics
			case "{input.bam}" in
				*\ *)
					sambamba merge -t {threads} /lscratch/$SLURM_JOB_ID/{wildcards.sample}.bam {input.bam}
					bammarkduplicates markthreads={threads} level=6 verbose=0 \
						M=sample_bam/markDupMetrics/{wildcards.sample}.markDup.methrics.tsv \
						tmpfile=/lscratch/$SLURM_JOB_ID/{wildcards.sample} \
						I=/lscratch/$SLURM_JOB_ID/{wildcards.sample}.bam \
						O={output.bam} \
						index=1 indexfilename={output.bai}
					;;
				*)
					bammarkduplicates markthreads={threads} level=6 verbose=0 \
						M=sample_bam/markDupMetrics/{wildcards.sample}.markDup.methrics.tsv \
						tmpfile=/lscratch/$SLURM_JOB_ID/{wildcards.sample} \
						I={input.bam} \
						O={output.bam} \
						index=1 indexfilename={output.bai}
					;;
			esac
			"""
#2/28/21 added (then removed again) --ignoreUnmated to samblaster because of error Can't find first and/or second of pair in sam block of length 1 for id: C7F3HANXX:5:1112:14712:22259
#5/17/2022: add --remove-duplicates to sambamba markdup so that manta will not produce error? --Used manta/1.6.0-fork-jykr.
#5/18/2022: biobambam2 is much faster than sambamba markdup (faster than 5x)
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
	threads: 4
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
		bam = 'sample_bam/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}.markDup.bai'
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
		mosdepth -t {threads} --no-per-base --by {config[bed]} --use-median --mapq 0 --fast-mode --thresholds 10,20,30 \
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

localrules: mean_coverage
rule mean_coverage:
	input:
		expand('coverage/mosdepth/{sample}.md.mosdepth.summary.txt', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		'coverage/mean.coverage.done.txt'
	shell:
		"""
		echo -e "sample\tlength\tmean" > coverage/{config[analysis_batch_name]}.mean.coverage.summary.tsv
		for file in {input}; do
			filename=$(basename $file)
			sm=$(echo $filename | sed 's/.md.mosdepth.summary.txt//')
			tail -n 1 $file | cut -f 2,4 | sed 's/^/'"$sm"'\t/' >> coverage/{config[analysis_batch_name]}.mean.coverage.summary.tsv
 		done
		touch {output}
		"""


# 30% smaller!
rule bam_to_cram:
	input:
		bam = 'sample_bam/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}.markDup.bai'
	output:
		cram = 'bam/{sample}.cram',
		crai = 'bam/{sample}.crai'
	threads:
		8
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
	shell:
		"""
		if [[ $(module list 2>&1 | grep "annovar" | wc -l) < 1 ]]; then module load {config[annovar_version]}; fi
		if [[ $(module list 2>&1 | grep "R/" | wc -l) < 1 ]]; then module load {config[R_version]}; fi
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
			module load {config[annotsv_version]} {config[crossmap_version]}
			tail -n +2 {input.deletion} | awk -F"\t" 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,"DEL"}}' > {input.deletion}.bed
			AnnotSV -genomeBuild {config[genomeBuild]} -SVinputFile {input.deletion}.bed -SVinputInfo 0 -svtBEDcol 4 -outputFile {output.del_anno}.temp
			Rscript /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/scramble_del_edit.R {output.del_anno}.temp {config[scrambleDELdb]} {output.del_anno}.temp
			rm {output.del_anno}.temp
		fi
		"""
#Try this strategy,if not working well, consider split by chr as in the GATK.
#freebayes_phasing finished in 1-2 hours for exome, if taking more than 4 hours, there's likely node problem
rule freebayes_phasing:
	input:
		bam = 'sample_bam/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}.markDup.bai'
	output:
		vcf = 'freebayes/vcf/{sample}.vcf.gz',
		filteredvcf = temp('freebayes/vcf/{sample}.filtered.vcf.gz'),
		tbi = temp('freebayes/vcf/{sample}.filtered.vcf.gz.tbi'),
		phasedvcf = 'freebayes/vcf/{sample}.phased.vcf.gz',
		phasedvcf_tbi = 'freebayes/vcf/{sample}.phased.vcf.gz.tbi'
	threads: 18
	shell:
		"""
		module load {config[freebayes_version]}
		module load {config[vcflib_version]}
		module load {config[samtools_version]}
		module load {config[vt_version]}
		freebayes-parallel {config[freebayes_exome_region]} $(({threads}-2)) -f {config[ref_genome]} \
			--limit-coverage 1000 {input.bam} --min-alternate-fraction 0.05 \
			--min-mapping-quality 1 --genotype-qualities --strict-vcf --use-mapping-quality \
			| bgzip -f > {output.vcf}
		sleep 2
		tabix -f -p vcf {output.vcf}
		if [[ {config[real_bed]} == "TRUE" ]]; then
			bcftools filter --regions-file {config[padded_bed]} --output-type v {output.vcf} \
				| bcftools norm --multiallelics -any --output-type v - \
				| vt decompose_blocksub -p -m -d 2 - \
				| bcftools norm --check-ref s --fasta-ref {config[ref_genome]} --output-type v - \
				| bcftools norm -d exact --output-type v - \
				| vcffilter -f "( QUAL > 15 & QA / AO > 15 & SAF > 0 & SAR > 0 & RPR > 0 & RPL > 0 & AO > 2 & DP > 3 ) | ( QUAL > 30 & QA / AO > 25 & ( SAF = 0 | SAR = 0 | RPR = 0 | RPL = 0 ) & AO > 2 & DP > 3 )" \
				| bgzip -f > {output.filteredvcf}
		else
			bcftools norm --multiallelics -any --output-type v {output.vcf} \
				| vt decompose_blocksub -p -m -d 2 - \
				| bcftools norm --check-ref s --fasta-ref {config[ref_genome]} --output-type v - \
				| bcftools norm -d exact --output-type v - \
				| vcffilter -f "( QUAL > 15 & QA / AO > 15 & SAF > 0 & SAR > 0 & RPR > 0 & RPL > 0 & AO > 2 & DP > 3 ) | ( QUAL > 30 & QA / AO > 25 & ( SAF = 0 | SAR = 0 | RPR = 0 | RPL = 0 ) & AO > 2 & DP > 3 )" \
				| bgzip -f > {output.filteredvcf}
		fi
		sleep 2
		tabix -f -p vcf {output.filteredvcf}
		module unload {config[freebayes_version]}
		module unload {config[vcflib_version]}
		module unload {config[vt_version]}
		module load {config[whatshap_version]}
		whatshap phase --reference {config[ref_genome]} --indels {output.filteredvcf} {input.bam} | bgzip -f > {output.phasedvcf}
		tabix -f -p vcf {output.phasedvcf}
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
			--minsize 1 --minvar 20
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

#Automap requires at least $10k variants passed QC filters, thus panel data fail frequently. Redirect by ">/dev/null 2>&1 || touch {output.tsv}" does not show error in sinteractive, but error in snakemake run.
		# [ -s file.txt ]: check whether a file is empty
rule merge_freebayes:
	input:
		vcf = expand('freebayes/vcf/{sample}.phased.vcf.gz', sample=list(SAMPLE_LANEFILE.keys())),
		tbi = expand('freebayes/vcf/{sample}.phased.vcf.gz.tbi', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		'freebayes/freebayes.merge.done.txt'
	threads: 8
	shell:
		"""
		module load {config[samtools_version]}
		case "{input.vcf}" in
			*\ *)
				bcftools merge --merge none --missing-to-ref --output-type z --threads {threads} {input.vcf} \
				> freebayes/{config[analysis_batch_name]}.freebayes.vcf.gz
				tabix -f -p vcf freebayes/{config[analysis_batch_name]}.freebayes.vcf.gz
				;;
			*)
				cp -p -l {input.vcf} freebayes/{config[analysis_batch_name]}.freebayes.vcf.gz
				cp -p -l {input.tbi} freebayes/{config[analysis_batch_name]}.freebayes.vcf.gz.tbi
				;;
		esac
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
		N_SHARDS="32"
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
			--intermediate_results_dir $WORK_DIR \
			--call_variants_extra_args="use_openvino=true"
		cd $PROJECT_WD
		cp $WORK_DIR/$(basename {output.vcf})* deepvariant/vcf
		cp $WORK_DIR/$(basename {output.gvcf})* deepvariant/gvcf
		module unload {config[deepvariant_version]}
		module load {config[samtools_version]}
		bcftools norm --multiallelics -any --output-type u {output.vcf} \
			| bcftools norm -d exact --output-type u - \
			| bcftools filter --threads $(({threads}-4)) --include 'FILTER="PASS" & FORMAT/AD[0:1]>2' --output-type z --output {output.filteredvcf}
		sleep 2
		tabix -f -p vcf {output.filteredvcf}
		module load {config[whatshap_version]}
		whatshap phase --reference {config[ref_genome]} --indels {output.filteredvcf} {input.bam} | bgzip -f > {output.phasedvcf}
		tabix -f -p vcf {output.phasedvcf}
		"""
#deepvariant PASS filter requires Alt AD > 1. I used > 2 for more stringent filtering.

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

rule glnexus:
	input:
		vcf = expand('deepvariant/gvcf/{sample}.dv.g.vcf.gz', sample=list(SAMPLE_LANEFILE.keys())),
		#bam = expand('sample_bam/{sample}.markDup.bam', sample=list(SAMPLE_LANEFILE.keys())),
		#bai = expand('sample_bam/{sample}.markDup.bai', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		'deepvariant/deepvariant.gvcf.merge.done.txt'
	threads: 24
	shell:
		"""
		module load {config[glnexus_version]} {config[samtools_version]}
		WORK_DIR="/lscratch/${{SLURM_JOB_ID}}"
		glnexus_cli --dir /lscratch/$SLURM_JOB_ID/glnexus --config DeepVariantWES --bed {config[padded_bed]} \
			--threads {threads}  --mem-gbytes 96 \
			{input.vcf} \
			| bcftools norm --multiallelics -any --output-type u --no-version \
			| bcftools norm --check-ref s --fasta-ref {config[ref_genome]} --output-type u --no-version - \
			| bcftools +fill-tags - -Ou -- -t AC,AC_Hom,AC_Het,AN,AF \
			| bcftools annotate --threads {threads} --set-id 'dvg_%CHROM\:%POS%REF\>%ALT' --no-version - -Oz -o deepvariant/{config[analysis_batch_name]}.glnexus.vcf.gz
		tabix -f -p vcf deepvariant/{config[analysis_batch_name]}.glnexus.vcf.gz
		touch {output}
		"""
#glnexus arranges samples in the multi-sample vcf files based on the inputs in the command line; bcftools merge sorts samples.
#dv_whatshap done in 20-60min for BP exome.
rule dv_whatshap:
	input:
		glnexus = 'deepvariant/deepvariant.gvcf.merge.done.txt',
		bam = 'sample_bam/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}.markDup.bai'
	output:
		phasedvcf = 'deepvariant/vcf/{sample}.dv.glnexus.phased.vcf.gz',
		phasedtbi = 'deepvariant/vcf/{sample}.dv.glnexus.phased.vcf.gz.tbi'
	threads: 4
	shell:
		"""
		module load {config[samtools_version]} {config[whatshap_version]}
		cp deepvariant/{config[analysis_batch_name]}.glnexus.vcf.gz* /lscratch/$SLURM_JOB_ID
		bcftools view --threads {threads} -Oz --samples {wildcards.sample} /lscratch/$SLURM_JOB_ID/{config[analysis_batch_name]}.glnexus.vcf.gz \
			-o /lscratch/$SLURM_JOB_ID/{wildcards.sample}.vcf.gz
		tabix -f -p vcf /lscratch/$SLURM_JOB_ID/{wildcards.sample}.vcf.gz
		whatshap phase --reference {config[ref_genome]} --indels /lscratch/$SLURM_JOB_ID/{wildcards.sample}.vcf.gz {input.bam} | bgzip -f -@ {threads} > {output.phasedvcf}
		tabix -f -p vcf {output.phasedvcf}
		"""

localrules: merge_glnexus_phased_vcf
rule merge_glnexus_phased_vcf:
	input:
		vcf = expand('deepvariant/vcf/{sample}.dv.glnexus.phased.vcf.gz', sample=list(SAMPLE_LANEFILE.keys())),
		tbi = expand('deepvariant/vcf/{sample}.dv.glnexus.phased.vcf.gz.tbi', sample=list(SAMPLE_LANEFILE.keys()))
	output:
		'deepvariant/deepvariant.glnexus.phased.merge.done.txt'
	threads: 8
	shell:
		"""
		if [[ $(module list 2>&1 | grep "samtools" | wc -l) < 1 ]]; then module load {config[samtools_version]}; fi
		case "{input.vcf}" in
			*\ *)
				bcftools merge --merge none --output-type z --threads {threads} {input.vcf} \
					> deepvariant/{config[analysis_batch_name]}.dv.glnexus.phased.vcf.gz
				sleep 2
				tabix -f -p vcf deepvariant/{config[analysis_batch_name]}.dv.glnexus.phased.vcf.gz
				;;
			*)
				cp -p -l {input.vcf} deepvariant/{config[analysis_batch_name]}.dv.glnexus.phased.vcf.gz
				cp -p -l {input.tbi} deepvariant/{config[analysis_batch_name]}.dv.glnexus.phased.vcf.gz.tbi
				;;
		esac
		touch {output}
		"""

localrules: merge_dv_fb_vcfs
rule merge_dv_fb_vcfs:
	input:
		'deepvariant/deepvariantVcf.merge.done.txt',
		'deepvariant/deepvariant.glnexus.phased.merge.done.txt',
		'freebayes/freebayes.merge.done.txt'
	output:
		'prioritization/dv_fb.merge.done.txt'
	threads: 8
	shell:
		"""
		if [[ $(module list 2>&1 | grep "samtools" | wc -l) < 1 ]]; then module load {config[samtools_version]}; fi
		WORK_DIR=/lscratch/$SLURM_JOB_ID
		bcftools isec -p $WORK_DIR/dv -w 2 --collapse none --output-type u --threads {threads} \
			deepvariant/{config[analysis_batch_name]}.dv.glnexus.phased.vcf.gz \
			deepvariant/{config[analysis_batch_name]}.dv.phased.vcf.gz
		bcftools +fill-tags $WORK_DIR/dv/0001.bcf -Ov -- -t AC,AC_Hom,AC_Het,AN,AF \
			| sed 's#0/0:.:.:.#0/0:10:10:10,0#g' - \
			| bcftools annotate --threads {threads} --set-id 'dv_%CHROM\:%POS%REF\>%ALT' --no-version - -Oz -o $WORK_DIR/dv/dv.hf.vcf.gz
		tabix -f -p vcf $WORK_DIR/dv/dv.hf.vcf.gz
		bcftools concat --threads {threads} -a --rm-dups none --no-version \
			deepvariant/{config[analysis_batch_name]}.dv.glnexus.phased.vcf.gz $WORK_DIR/dv/dv.hf.vcf.gz -Oz \
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
# if [[ $(module list 2>&1 | grep "samtools" | wc -l) < 1 ]]; then module load {config[samtools_version]}; fi
# WORK_DIR=/lscratch/$SLURM_JOB_ID
# bcftools isec --threads {threads} -p $WORK_DIR --collapse none --no-version -Oz \
# 	deepvariant/{config[analysis_batch_name]}.dv.glnexus.phased.vcf.gz \
# 	freebayes/{config[analysis_batch_name]}.freebayes.vcf.gz
# rm $WORK_DIR/0003.vcf* &
# bcftools annotate --threads {threads} --set-id 'dv_%CHROM\:%POS%REF\>%ALT' \
# 	--no-version $WORK_DIR/0000.vcf.gz -Oz -o $WORK_DIR/dv.vcf.gz && rm $WORK_DIR/0000.vcf* &
# bcftools annotate --threads {threads} --set-id 'fb_%CHROM\:%POS%REF\>%ALT' -x ^INFO/QA,FORMAT/RO,FORMAT/QR,FORMAT/AO,FORMAT/QA,FORMAT/GL \
# 	--no-version $WORK_DIR/0001.vcf.gz -Ou - \
# 	| bcftools +fill-tags - -Ov -- -t AC,AC_Hom,AC_Het,AN,AF \
# 	| sed 's#0/0:.:.:.#0/0:10:10:10,0#g' - \
# 	| bgzip -f > $WORK_DIR/fb.vcf.gz
# rm $WORK_DIR/0001.vcf* &
# bcftools annotate --threads {threads} --set-id 'dvFb_%CHROM\:%POS%REF\>%ALT' \
# 	--no-version $WORK_DIR/0002.vcf.gz -Oz -o $WORK_DIR/dvFb.vcf.gz
# rm $WORK_DIR/0002.vcf* &
# tabix -f -p vcf $WORK_DIR/dv.vcf.gz
# tabix -f -p vcf $WORK_DIR/fb.vcf.gz
# tabix -f -p vcf $WORK_DIR/dvFb.vcf.gz
# bcftools concat --threads {threads} -a --rm-dups none --no-version \
# 	$WORK_DIR/dvFb.vcf.gz $WORK_DIR/dv.vcf.gz $WORK_DIR/fb.vcf.gz -Oz \
# 	-o prioritization/{config[analysis_batch_name]}.vcf.gz
# tabix -f -p vcf prioritization/{config[analysis_batch_name]}.vcf.gz
# if [[ {config[genomeBuild]} == "GRCh38" ]]; then
# 	module load {config[crossmap_version]}
# 	hg19ref=/data/OGL/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta
# 	crossmap vcf /data/OGL/resources/ucsc/hg38ToHg19.over.chain.gz \
# 		prioritization/{config[analysis_batch_name]}.vcf.gz \
# 		$hg19ref \
# 		$WORK_DIR/GRCh37.vcf
# 	sed -e 's/^chrM/MT/' -e 's/<ID=chrM/<ID=MT/' $WORK_DIR/GRCh37.vcf \
# 		| sed -e 's/^chr//' -e 's/<ID=chr/<ID=/' - \
# 		| bcftools norm --check-ref s --fasta-ref $hg19ref --output-type u - \
# 		| bcftools sort -m 20G -T $WORK_DIR/ -Ou - \
# 		| bcftools norm --threads $(({threads}-4)) -d exact --output-type z - -o prioritization/{config[analysis_batch_name]}.GRCh37.vcf.gz
# 	tabix -f -p vcf prioritization/{config[analysis_batch_name]}.GRCh37.vcf.gz
# fi

#exome some MT sequences were covered based on mosdepth. 5/21/2022.
#bcftools 0: private to the first vcf; 1: private to the 2nd vcf; 2: records shared by both from 1st vcf; 3: records shared by both from 2nd vcf.
#<30min, 3gb max mem for wgs
rule manta:
	input:
 		bam = 'sample_bam/{sample}.markDup.bam',
 		bai = 'sample_bam/{sample}.markDup.bai'
	output:
		'manta/manta.{sample}.annotated.tsv'
	threads: 32
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
		cp $RUNDIR/manta.{wildcards.sample}.annotated.tsv {output}
		"""

rule bcm_locus:
	input:
		bam = 'sample_bam/{sample}.markDup.bam',
		bai = 'sample_bam/{sample}.markDup.bai'
	output:
		bam = 'bam/bcmlocus/{sample}.bcm.bam',
		bai = 'bam/bcmlocus/{sample}.bcm.bai',
		vcf = temp('bcmlocus/{sample}.vcf'),
		avinput = temp('bcmlocus/{sample}.avinput'),
		bcm_out = 'bcmlocus/{sample}.bcmlocus.tsv'
	threads: 8
	shell:
		"""
		export TMPDIR=/lscratch/$SLURM_JOB_ID
		module load {config[samtools_version]} {config[bazam_version]} {config[bwa-mem2_version]} {config[samblaster_version]} {config[sambamba_version]}
		RG=$(samtools view -H {input.bam} | grep "^@RG" | head -n 1 | sed 's/\t/\\\\t/g')
		java -Xmx16g -jar $BAZAMPATH/bazam.jar -bam {input.bam} --regions chrX:153929000-154373500 \
		| bwa-mem2 mem -t $(({threads}/2)) -K 100000000 -M -Y -B 4 -O 6 -E 1 -p -R $RG {config[GRCh38Decoy2]} - \
		| samblaster -M --addMateTags --quiet \
		| sambamba sort -u --tmpdir=/lscratch/$SLURM_JOB_ID -t $(({threads}/2)) -o {output.bam} \
			<(sambamba view -S -f bam -l 0 -t $(({threads}/2)) /dev/stdin)
		mv {output.bam}.bai {output.bai}
		if [[ $(module list 2>&1 | grep "mosdepth" | wc -l) < 1 ]]; then module load {config[mosdepth_version]}; fi
		if [[ $(module list 2>&1 | grep "R/" | wc -l) < 1 ]]; then module load {config[R_version]}; fi
		mkdir -p bcmlocus/mosdepth
		cd bcmlocus/mosdepth
		mosdepth -t {threads} --no-per-base --by {config[bcmlocus_bed]}  --use-median --mapq 0 --fast-mode \
			{wildcards.sample}.md ../../{output.bam}
		cd ../..
		if [[ $(module list 2>&1 | grep "samtools" | wc -l) < 1 ]]; then module load {config[samtools_version]}; fi
		if [[ $(module list 2>&1 | grep "freebayes" | wc -l) < 1 ]]; then module load {config[freebayes_version]}; fi
		if [[ $(module list 2>&1 | grep "annovar" | wc -l) < 1 ]]; then module load {config[annovar_version]}; fi
		freebayes -f {config[GRCh38Decoy2]} --max-complex-gap 90 -p 6 -C 3 -F 0.05 \
			--genotype-qualities --strict-vcf --use-mapping-quality \
			--targets /data/OGL/resources/bed/OPN1LWe2e5.bed \
			{output.bam} \
			| bcftools norm --multiallelics -any --output-type u - \
			| bcftools annotate --set-id '%CHROM\:%POS%REF\>%ALT' -x ^INFO/AF --output-type u --no-version \
			| bcftools norm --check-ref s --fasta-ref {config[GRCh38Decoy2]} --output-type u --no-version - \
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
		if [[ $(module list 2>&1 | grep "picard" | wc -l) < 1 ]]; then module load {config[picard_version]}; fi
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
