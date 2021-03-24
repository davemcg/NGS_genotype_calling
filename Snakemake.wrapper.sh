#!/bin/bash
#SBATCH --gres=lscratch:100
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g

# to run snakemake as batch job
# run in the data folder for this project, fastq files must be in the folder fastq.
# sbatch --cpus-per-task=8 --mem=32g --time=12:0:0 ~/git/NGS_genotype_calling/Snakemake.wrapper.sh config_panel.yaml OGL733v1 panel
# $2: Libary Information, such as OGLv1
# $3: If typing in "panel", the panel.Snakefile and panel.cluster.json will be run. If "exome", then exome snakefile/exome.json will be used. Anything else including empty will run the WGS pipeline. The config file is the $1.
# $4: --notemp --dryrun --unlock or empty
#When there is two or more *metadata_file*.csv present in the folder, then -e *metadata_file.csv will produce "binary operator expected". Thus changed to only single file.

cp /data/OGL/resources/NGS_genotype_calling.git.log .
mkdir -p 00log
module load snakemake/5.24.1 || exit 1
#current version 5.7.4
sbcmd="sbatch --cpus-per-task={threads} \
--mem={cluster.mem} \
--time={cluster.time} \
--partition={cluster.partition} \
--output={cluster.output} \
--error={cluster.error} \
{cluster.extra}"

lib=$2
ngstype=$3

if [ -e metadata_file.csv ];
then
	echo "metadata_file provided"
else
	for fastq1 in fastq/*.gz; do
	filename=$(basename $fastq1)
	header=$(zcat $fastq1 | head -1)
	id=$(echo $header | cut -d: -f 3,4 | sed 's/\:/\./g')
	sm=$(echo $filename | cut -d_ -f 1 | sed 's/\-/\_/g')
	echo "$sm,$filename,@RG\\\tID:$id"_"$sm\\\tSM:$sm\\\tLB:$lib"_"$sm\\\tPL:ILLUMINA" >> metadata_file.csv
	done
fi
#if sample name has underscore, then sm=$(echo $filename | cut -d_ -f 1,2)
#https://software.broadinstitute.org/gatk/documentation/article.php?id=6472

#id=$(echo $header | cut -d: -f 3,4,10 | sed 's/\:/\./g') edited 8/27/19
#echo "$sm,$filename,@RG\\\tID:$id\\\tSM:$sm\\\tLB:$lib"_"$sm\\\tPL:ILLUMINA" edited 8/27/19
#removed R1_001 from for fastq1 in fastq/*R1_001.fastq.gz; 7/9/19
#RG information: https://software.broadinstitute.org/gatk/documentation/article?id=6472
#https://support.sentieon.com/appnotes/read_groups/
#id added sample no field which is in position 10 in the fastq file; When working with another Instrument, check and see whether the id field will be unique.
#pu is removed
#pu=$(echo $header | cut -d: -f 3,4,10 | sed 's/\:/\./g');

# restart-times changed to 1 because CREST had job failure.
case "${ngstype^^}" in
	"PANEL")
		snakemake -s /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/panel.Snakefile \
		-pr --local-cores 2 --jobs 1999 \
		--cluster-config /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/panel.cluster.json \
		--cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
		-k --restart-times 1 \
		--resources parallel=4 \
		--configfile $1 $4
		;;
	"EXOME"|"WES")
		snakemake -s /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/exome.Snakefile \
		-pr --local-cores 2 --jobs 1999 \
		--cluster-config /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/exome.cluster.json \
		--cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
		-k --restart-times 1 \
		--resources parallel=4 \
		--configfile $1 $4
		;;
	*)
		snakemake -s /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/Snakefile \
		-pr --local-cores 2 --jobs 1999 \
		--cluster-config /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/cluster.json \
		--cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
		-k --restart-times 1 \
		--resources parallel=4 \
		--configfile $1 $4
		;;
esac

# --notemp Ignore temp() declaration;
# --dryrun
# --unlock
# --dag


# if [[ "${ngstype^^}" =~ ^(PANEL|EXOME)$ ]];
# then
	# snakemake -s /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/panel.Snakefile \
	# -pr --local-cores 2 --jobs 1999 \
	# --cluster-config /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/panel.cluster.json \
	# --cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
	# -k --restart-times 4 \
	# --resources parallel=4 \
	# --configfile $1
# else
	# snakemake -s /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/Snakefile \
	# -pr --local-cores 2 --jobs 1999 \
	# --cluster-config /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/cluster.json \
	# --cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
	# -k --restart-times 4 \
	# --resources parallel=4 \
	# --configfile $1
# fi
