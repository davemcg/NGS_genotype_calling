#!/bin/bash
#SBATCH --gres=lscratch:100
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --partition=quick
#SBATCH --time=1:0:0

config=$1
sample=$2

#selecting mapped reads using subset_mapped_reads.sh reduced sample processing time to 44 hours.
#selecting reads using subset_chr1-M_reads.sh: 40 hours (without P flag. Adding P flag makes the run > 48h).

set -e
module load $(grep "^samtools_version:" $config | head -n 1 | cut -d"'" -f 2)

export TMPDIR=/lscratch/$SLURM_JOB_ID
module load bazam/1.0.1 bwa/0.7.17 samblaster/0.1.26 sambamba/1.0.0
mkdir -p sample_bam/subset/
#samtools view --threads $SLURM_CPUS_PER_TASK -b sample_bam/$sample.markDup.bam --output sample_bam/subset/$sample.markDup.chrX.bam chrX
#samtools index -@ $SLURM_CPUS_PER_TASK sample_bam/subset/$sample.markDup.chrX.bam

RG=$(samtools view -H sample_bam/"$sample".markDup.bam | grep "^@RG" | head -n 1 | sed 's/\t/\\t/g')
echo $RG
echo "Start realignment"
date
java -Xmx32g -jar $BAZAMPATH/bazam.jar -bam sample_bam/$sample.markDup.bam --regions chrX:153929000-154373500 | bwa mem -t 16 -M -Y -B 4 -O 6 -E 1 -p -R $RG /data/OGL/resources/genomes/GRCh38/GRCh38Decoy2.fa - | samblaster -M --addMateTags --quiet | sambamba sort -u --tmpdir=/lscratch/$SLURM_JOB_ID -t 16 -o bcmlocus/bam/$sample.bcm.bam <(sambamba view -S -f bam -l 0 -t 16 /dev/stdin)
date
echo "Realignment completed"

if [[ $(module list 2>&1 | grep "mosdepth" | wc -l) -lt 1 ]]; then module load mosdepth/0.3.3; fi
mkdir -p bcmlocus/mosdepth
cd bcmlocus/mosdepth
mosdepth -t 20 --no-per-base --by /data/OGL/resources/bcm/xgen-v2-chrX154mb_v1.bed --use-median --mapq 0 --fast-mode $sample.md ../../bcmlocus/bam/$sample.bcm.bam
cd ../..
if [[ $(module list 2>&1 | grep "freebayes" | wc -l) -lt 1 ]]; then module load freebayes/1.3.5; fi
if [[ $(module list 2>&1 | grep "annovar" | wc -l) -lt 1 ]]; then module load annovar/2020-06-08; fi
if [[ $(module list 2>&1 | grep "vcflib" | wc -l) -lt 1 ]]; then module load vcflib/1.0.3; fi
if [[ $(module list 2>&1 | grep "vt/" | wc -l) -lt 1 ]]; then module load vt/0.57721; fi
freebayes -f /data/OGL/resources/genomes/GRCh38/bwa-mem2/GRCh38Decoy2.fa --max-complex-gap 80 -p 10 -C 3 -F 0.05 --genotype-qualities --strict-vcf --use-mapping-quality --targets /data/OGL/resources/bed/OPN1LWe2e5.bed bcmlocus/bam/$sample.bcm.bam | vcffilter -f "QUAL > 10" | bcftools +fill-tags - -Ou -- -t VAF | bcftools norm --multiallelics -any --output-type u - | bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/GRCh38/bwa-mem2/GRCh38Decoy2.fa  --output-type u --no-version - | bcftools annotate --set-id 'haplo_%CHROM\:%POS%REF\>%ALT' -x ^INFO/AO,^FORMAT/GT,FORMAT/DP,FORMAT/VAF --output-type u --no-version | bcftools norm -d exact --output-type z -o bcmlocus/vcf/$sample.haplo.vcf.gz
tabix -p vcf bcmlocus/vcf/$sample.haplo.vcf.gz
freebayes -f /data/OGL/resources/genomes/GRCh38/bwa-mem2/GRCh38Decoy2.fa --limit-coverage 1000 --min-alternate-fraction 0.05 -C 3 -p 10 --min-mapping-quality 0 --genotype-qualities --strict-vcf --use-mapping-quality --targets /data/OGL/resources/bed/OPN1LWe1e6_MWe1.bed bcmlocus/bam/$sample.bcm.bam | vcffilter -f "QUAL > 10 & SAF > 0 & SAR > 0 & RPR > 0 & RPL > 0 & AO > 2 & DP > 5" | bcftools +fill-tags - -Ou -- -t VAF | bcftools norm --multiallelics -any --output-type u --force - | vt decompose_blocksub -p -m -d 2 - | bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/GRCh38/bwa-mem2/GRCh38Decoy2.fa --output-type u - | bcftools annotate --set-id 'sml_%CHROM\:%POS%REF\>%ALT' -x ^INFO/AO,^FORMAT/GT,FORMAT/DP,FORMAT/VAF --output-type v --no-version | bcftools norm -d exact --output-type z -o bcmlocus/vcf/$sample.small.vcf.gz
tabix -p vcf bcmlocus/vcf/$sample.small.vcf.gz
bcftools concat -a --rm-dups none --no-version bcmlocus/vcf/$sample.small.vcf.gz bcmlocus/vcf/$sample.haplo.vcf.gz -Oz -o bcmlocus/vcf/$sample.vcf.gz
tabix -p vcf bcmlocus/vcf/$sample.vcf.gz
