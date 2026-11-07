#!/bin/bash
#SBATCH --gres=lscratch:50
#SBATCH --cpus-per-task=4
#SBATCH --mem=8g
#SBATCH --partition=norm,quick
#SBATCH --time=4:0:0

sample=$1
module load nextflow/25.10.0 samtools
export NXF_SINGULARITY_CACHEDIR=/data/OGL/resources/git/mtdna-server-2/work/singularity
WORK_DIR=$PWD
LWORK_DIR=/lscratch/$SLURM_JOB_ID
MUTSERVE_HOME="/data/OGL/resources/git/mutserve"
samtools view -@ 4 -b $WORK_DIR/sample_bam/$sample.markDup.bam chrM -o $LWORK_DIR/$sample.chrM.bam
samtools index -@ 4 $LWORK_DIR/$sample.chrM.bam
cd /lscratch/$SLURM_JOB_ID
#sed -e "s#bamfile#$WORK_DIR/sample_bam/$sample.markDup.bam#" /data/OGL/resources/git/mtdna-server-2.config > mtdna-server-2.config
sed -e "s#bamfile#$sample.chrM.bam#" -e "s#mtDNA-server-2#mtDNA-server-2/$sample#" /data/OGL/resources/git/mtdna-server-2.config > mtdna-server-2.$sample.config

nextflow run /data/OGL/resources/git/mtdna-server-2/main.nf -c mtdna-server-2.$sample.config \
  --output_reports mtDNA-server-2/$sample \
  -profile singularity
rm -rf mtDNA-server-2/$sample/mutserve_raw_files $sample.chrM.bam*
rm -rf $WORK_DIR/mito/$sample/*
mkdir -p $WORK_DIR/mito/$sample $WORK_DIR/mito/mitomap/

mv mtDNA-server-2/$sample/* $WORK_DIR/mito/$sample
cd $WORK_DIR

awk -F'\t' -v OFS='\t' '
NR==FNR {
    if (FNR > 1) {
        gsub(/^"|"$/, "", $2)
        gsub(/^"|"$/, "", $3)
        gsub(/^"|"$/, "", $4)

        hap = hap (hap ? ";" : "") $2
        rank = rank (rank ? ";" : "") $3
        qual = qual (qual ? ";" : "") $4
    }
    next
}

FNR == 1 {
    print $0, "\"Haplogroup\"", "\"Rank\"", "\"Quality\""
    next
}

{
    print $0, "\"" hap "\"", "\"" rank "\"", "\"" qual "\""
}
' mito/$sample/haplogroups.txt mito/$sample/haplocheck.txt \
 | sed 's/\.markDup\.bam//' > mito/$sample.haplocheck.txt
 
module load vcfanno R/4.3.0

awk -F'\t' 'BEGIN {
    OFS="\t"
    print "##fileformat=VCFv4.2"
    print "##contig=<ID=chrM,length=16569>"
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">"
    print "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"
    print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","S1"
}
NR > 1 {
    id="chrM-" $3 "-" $4 "-" $5
    print "chrM",$3,id,$4,$5,"30","PASS",".","GT:GQ:DP","0/1:50:100"
}' mito/$sample/variants.annotated.txt \
  | bcftools norm --multiallelics -any --output-type u - \
  | bcftools norm --check-ref w --fasta-ref $MUTSERVE_HOME/rCRS.fasta --output-type v - \
  | bgzip -f > $LWORK_DIR/$sample.mt.norm.vcf.gz
tabix -p vcf -f $LWORK_DIR/$sample.mt.norm.vcf.gz
vcfanno ~/git/NGS_genotype_calling/NGS_generic_OGL/mitomap.conf $LWORK_DIR/$sample.mt.norm.vcf.gz \
 > $LWORK_DIR/$sample.mt.mm.vcf
bcftools query -H -f '%ID\t%INFO/Pubmed\t%INFO/mitomapDisease\t%INFO/mitomapDiseaseStatus\n' $LWORK_DIR/$sample.mt.mm.vcf \
 | sed -e '1 s/^#//' -e '1 s/\[[0-9]*\]//g' \
 > $LWORK_DIR/$sample.mt.mm.tsv
Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/mtDNA-server2.mitomap.R \
  mito/$sample/variants.annotated.txt \
  $LWORK_DIR/$sample.mt.mm.tsv \
  mito/mitomap/$sample.AnnotatedVariants.txt
  
#WARN: Access to undefined parameter `output_auxiliary` -- Initialise it to a default value eg. `params.output_auxiliary = some_value`

#34040344 submitted 11:40 10/27/2026