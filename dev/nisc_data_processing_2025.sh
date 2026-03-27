#!/bin/bash
#SBATCH --gres=lscratch:10
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --time=2:00:00

#Create a project folder, then globus transfer data to the project folder such as eyeGENE_ORF15
#run outside the project folder

ProjectDir=${1%/}

WORK_DIR=$PWD
if (( $(echo $WORK_DIR | grep "/data/OGL" | wc -l) > 0 )); then
 find $ProjectDir -type f ! -group OGL -print -exec chgrp OGL -- {} +
 chmod --recursive g+rw $ProjectDir
fi

TIMESTAMP=$(date "+%Y%m%d-%H%M%S")
mkdir -p $ProjectDir-$TIMESTAMP/old_bam $ProjectDir-$TIMESTAMP/vcf_nisc/gvcf $ProjectDir-$TIMESTAMP/vcf_nisc/vcf $ProjectDir-$TIMESTAMP/summary

module load samtools parallel

find $ProjectDir -type f -name "*.bam*" -exec mv {} $ProjectDir-$TIMESTAMP/old_bam \;
find $ProjectDir -type f -name "*.g.vcf.gz*" -exec mv {} $ProjectDir-$TIMESTAMP/vcf_nisc/gvcf \;
find $ProjectDir -type f -name "*.vcf*" -exec mv {} $ProjectDir-$TIMESTAMP/vcf_nisc/vcf \;
find $ProjectDir -type f -name "*.txt" -exec mv {} $ProjectDir-$TIMESTAMP/summary \; 
find $ProjectDir -type f -name "*.xls" -exec mv {} $ProjectDir-$TIMESTAMP/summary \;

#for vcf in $ProjectDir-$TIMESTAMP/vcf_nisc/vcf/*.vcf; do bgzip -@ 8 $vcf; tabix -@ 8 $vcf.gz; done
find "$ProjectDir-$TIMESTAMP/vcf_nisc/vcf" -type f -name '*.vcf' -print0 \
| parallel -0 --no-run-if-empty -j 8 '
    bgzip -@1 "{}" &&
    tabix -@1 -p vcf "{}.gz"
'

mv $ProjectDir/*.* $ProjectDir-$TIMESTAMP/
rm -rf $ProjectDir
mv $ProjectDir-$TIMESTAMP $ProjectDir
cd $ProjectDir
bash ~/git/NGS_genotype_calling/make_metadata_scripts/make_bam_metadatafile.sh


#to be continued
