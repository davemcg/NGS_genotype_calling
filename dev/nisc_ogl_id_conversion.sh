cd /data/OGL/genome/
mkdir -p nisc23-1/old_bam
find nisc_2023_02_transfer/ -name "*D*.bam*" -exec mv {} nisc23-1/old_bam/ \;
mv nisc2022Oct/old_bam/*.bam* nisc23-1/old_bam/
rm -r nisc2022Oct
cd nisc23-1/
bash ~/git/NGS_genotype_calling/make_bam_metadatafile.sh
#use R to make the file matching niscID with OGL_ID using the R script in the folder O:\Research-Projects\Research_NGS_Seq\2022-10-NISC
module load R/3.6.3
Rscript ~/git/NGS_genotype_calling/dev/nisc_ogl_id_conversion.R IDconversion.tsv metadata_file.csv metadata_file.1.csv
awk -F"," 'BEGIN{OFS=","} {sub($1, $4, $3); print $4,$2,$3 }' metadata_file.1.csv > metadata_file.csv
rm metadata_file.1.csv
