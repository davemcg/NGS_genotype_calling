#!/bin/bash
set -e
sample=$1
#ml parallel
#parallel -j 8 bash ~/git/NGS_genotype_calling/dev/scramble_anno.sh {} < sampleList.txt
#head -n 1 G4V9_1_BP93806.scramble.tsv > all.exome.scramble.txt && for file in *.scramble.tsv; do tail -n +2 $file >> all.exome.scramble.txt; done

if [[ $(module list 2>&1 | grep "annovar" | wc -l) < 1 ]]; then module load annovar/2019-10-24; fi
if [[ $(module list 2>&1 | grep "R/" | wc -l) < 1 ]]; then module load R/3.6.3; fi
if [[ GRCh38 == "GRCh38" ]]; then
	ver=hg38
else
	ver=hg19
fi
if [[ $(wc -l scramble/"$sample"_MEIs.txt | cut -d " " -f 1) == 1 ]]
then
	touch scramble_anno/"$sample".avinput
	touch scramble_anno/"$sample".forR.txt
	touch scramble_anno/"$sample".scramble.tsv
else
	cut -f 1 scramble/"$sample"_MEIs.txt | awk -F ":" 'BEGIN{OFS="\t"} NR>1 {print $1,$2,$2,"0","-"}' > scramble_anno/"$sample".avinput
	table_annovar.pl scramble_anno/"$sample".avinput $ANNOVAR_DATA/$ver -buildver $ver -remove -out scramble_anno/"$sample" --protocol refGene -operation g --argument '-splicing 100 -hgvs' --polish -nastring . --thread 1
	awk -F"\t" 'BEGIN{OFS="\t"} NR==1 {print "Func_refGene","Gene","Intronic","AA"} NR>1 {print $6,$7,$8,$10}' scramble_anno/"$sample"."$ver"_multianno.txt | paste scramble/"$sample"_MEIs.txt - > scramble_anno/"$sample".forR.txt
	rm scramble_anno/"$sample"."$ver"_multianno.txt
	Rscript /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/scramble_anno.R scramble_anno/"$sample".forR.txt /data/OGL/resources/scramble/SCRAMBLEvariantClassification.GRCh38.2022-09.xlsx /data/OGL/resources/OGLpanelGeneDxORcandidate.xlsx /data/OGL/resources/HGMD/HGMDtranscript.txt "$sample" scramble_anno/"$sample".scramble.tsv scramble_anno/"$sample".scramble.xlsx
fi
rm scramble_anno/"$sample".avinput
rm scramble_anno/"$sample".forR.txt