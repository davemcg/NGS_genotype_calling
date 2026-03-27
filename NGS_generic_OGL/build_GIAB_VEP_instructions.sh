##GIAB: https://www.nist.gov/programs-projects/genome-bottle


##Download the high-confidence truth set into folder NA12878/:
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/supplementaryFiles/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.bed
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/supplementaryFiles/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.vcf.gz
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh37/supplementaryFiles/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.vcf.gz.tbi

##make high-confidence bed files
module load bedtools
intersectBed -a NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.bed -b /data/OGVFB/OGL_NGS/bed/OGLv1.sorted.bed > OGLv1_GIAB_highconf.bed

##prepare rtg format
module load rtg
rtg format -o human_g1k_v37decoy.sdf /data/OGVFB/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta

sinteractive --mem=64g ##(16g is not sufficient for vcfeval)

rtg vcfeval -b NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.vcf.gz -c /data/guanb/miseq_analysis/0804_12samples/freebayes/14_NA12878.freebayes.filtered.vcf.gz -t human_g1k_v37decoy.sdf -e OGLv1_GIAB_highconf.bed --vcf-score-field QUAL -o OGLv1_NA12878_freebayes_test2 

rtg rocplot OGLv1_NA12878_freebayes_test2/snp_roc.tsv.gz --png=OGLv1_NA12878_freebayes_test2/snp_roc_freebayes.png
rtg rocplot OGLv1_NA12878_freebayes_test2/non_snp_roc.tsv.gz --png=OGLv1_NA12878_freebayes_test2/non_snp_roc_freebayes.png

##Did not produce ROC curves? Why? Becuase did not use -f argument.

rtg GATK file -


#GATK file, gvcf-to-vcf

#select variants:
GATK -m 16g SelectVariants \
 -R /data/guanb/resource/GATK_resource_bundle/b37-2.8/human_g1k_v37_decoy.fasta \
 -V GATK.filterSNP-INDEL.vcf.gz \
 -o NA12878.GATK.filterSNP-INDEL.vcf.gz \
 -sn 14_NA12878 \
 --excludeFiltered

rtg vcfeval -b NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.vcf.gz -c /data/guanb/miseq_analysis/0804_12samples/prioritization/NA12878.GATK.filterSNP-INDEL.vcf.gz -t human_g1k_v37decoy.sdf -e OGLv1_GIAB_highconf.bed -o GATK_NA12878

##default vcf-score-field was GQ above.

rtg vcfeval -b NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.vcf.gz -c /data/guanb/miseq_analysis/0804_12samples/prioritization/NA12878.GATK.filterSNP-INDEL.vcf.gz -t human_g1k_v37decoy.sdf -e OGLv1_GIAB_highconf.bed --vcf-score-field QUAL -o GATK_NA12878_test2


##Test 10/21/19 using vt
module load  vt/0.577
module load samtools/1.9
zcat /data/guanb/miseq_analysis/0804_12samples/freebayes/14_NA12878.freebayes.filtered.vcf.gz | vt decompose -s - | vt normalize -r /data/OGVFB/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta - | bgzip -c > 14_NA12878.freebayes.filtered.vt.vcf.gz

zcat NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.vcf.gz | vt decompose -s - | vt normalize -r /data/OGVFB/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta - | bgzip -c > NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.vt.vcf.gz
 
zcat /data/guanb/miseq_analysis/0804_12samples/prioritization/NA12878.GATK.filterSNP-INDEL.vcf.gz | vt decompose -s - | vt normalize -r /data/OGVFB/resources/1000G_phase2_GRCh37/human_g1k_v37_decoy.fasta - | bgzip -c > NA12878.GATK.filterSNP-INDEL.vt.vcf.gz


rtg vcfeval -b NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.vt.vcf.gz -c 14_NA12878.freebayes.filtered.vt.vcf.gz -t human_g1k_v37decoy.sdf -e OGLv1_GIAB_highconf.bed --vcf-score-field QUAL -o freebayes_1021

rtg vcfeval -b NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.vt.vcf.gz -c NA12878.GATK.filterSNP-INDEL.vt.vcf.gz -t human_g1k_v37decoy.sdf -e OGLv1_GIAB_highconf.bed --vcf-score-field QUAL -o GATK_1021

