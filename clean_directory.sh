#!/bin/bash
#SBATCH --gres=lscratch:10
#SBATCH --cpus-per-task=8
#SBATCH --mem=8g
#SBATCH --time=2:00:00

calling_configure_file=$1
variantPrioritization_configure_file=$2

set -e

TIMESTAMP=$(date "+%Y%m%d-%H%M%S")


#module load R/4.3.0 
module load parallel
analysis_batch_name=$(grep "^analysis_batch_name:" $calling_configure_file | head -n 1 | cut -d"'" -f 2)
ngstype=$(grep "^datatype:" $variantPrioritization_configure_file | head -n 1 | cut -d"'" -f 2)
ped_file=$(grep "^ped:" $variantPrioritization_configure_file | head -n 1 | cut -d"'" -f 2)

export ngstype

#copy sample names to OGL.ped adding datatype and analysis_batch_name
#Compile this file manually.
# awk -F"\t" -v batch="$analysis_batch_name" -v dt="$ngstype" 'BEGIN{OFS="\t"} NR==1 {print "NGS_type", "batch", $0} NR>1 {print dt, batch, $0}' \
 # prioritization/$ped_file | cat /data/OGL/sample_info/master.v3.ped - > /lscratch/$SLURM_JOB_ID/temp.ped
# awk '!a[$0]++' /lscratch/$SLURM_JOB_ID/temp.ped > /data/OGL/sample_info/master.v3.ped

#subsequent batch
# awk -F"\t" -v batch="$analysis_batch_name" -v dt="$ngstype" 'BEGIN{OFS="\t"} NR>1 {print $0, batch, dt}' prioritization/$ped_file >> /data/OGL/resources/OGLsample/OGL.ped

#chgrp OGL /data/OGL/sample_info/master.ped #what happens for another user?

#copy relevant files to resources

find prioritization/gemini_tsv_filtered/ -name "*.tsv" | parallel -j 8 'gzip -c {} \
 > /data/OGL/resources/GeneSearch/$ngstype/gemini_tsv_filtered/$(echo {/} | cut -d. -f 1).tsv.gz && chgrp OGL /data/OGL/resources/GeneSearch/$ngstype/gemini_tsv_filtered/$(echo {/} | cut -d. -f 1).tsv.gz' 


#chgrp --recursive OGL /data/OGL/resources/GeneSearch/genome/gemini_tsv_filtered/temp$TIMESTAMP
#mv /data/OGL/resources/GeneSearch/genome/gemini_tsv_filtered/temp$TIMESTAMP/* /data/OGL/resources/GeneSearch/genome/gemini_tsv_filtered/v3
#rm -r /data/OGL/resources/GeneSearch/genome/gemini_tsv_filtered/temp$TIMESTAMP



# ( cat $CONTIGFILE | parallel -C "\t" -j 21 --tmpdir $WORK_DIR --eta --halt 2 --line-buffer \
		 	# --tag "whatshap phase --reference $WORK_DIR/$(basename {config[ref_genome]}) \
			# --indels --ignore-read-groups $WORK_DIR/filtered/{{2}}.filtered.vcf.gz $WORK_DIR/{wildcards.sample}.markDup.bam \
			# | bgzip -f --threads $(({threads}-10)) > $WORK_DIR/phased/{{2}}.phased.vcf.gz" \
		# ) && echo "whatshap on chr completed" || exit 7

#vcf, SV files to resources

case "${ngstype^^}" in
 "PANEL")
  find prioritization/gemini_tsv_filtered/ -name "*.tsv.gz" | parallel -j 8 'cp {} /data/OGL/resources/GeneSearch/panel/gemini_tsv_filtered/$(echo {/} | cut -d. -f 1).tsv.gz && chgrp OGL /data/OGL/resources/GeneSearch/panel/gemini_tsv_filtered/$(echo {/} | cut -d. -f 1).tsv.gz' 
  ;;
 "AMPLICON"|"AMP")
  echo "Amplicon, do not retain files to resources"
  ;;
 "EXOME"|"WES"|"ES")
  ngstype="exome"
  find prioritization/gemini_tsv_filtered/ -type f -name "*.tsv.gz" | parallel -j 8 'destination_file="/data/OGL/resources/GeneSearch/$ngstype/gemini_tsv_filtered/$(echo {/} | cut -d. -f 1).tsv.gz"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  find prioritization/ -type f -name "*.gt3.anno3.dvg.vcf.gz*" | parallel -j 8 'destination_file="/data/OGL/resources/OGLsample/annotatedVCF/$ngstype/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  find deepvariant/gvcf/ -type f -name "*.vcf.gz*" | parallel -j 8 'destination_file="/data/OGL/resources/OGLsample/${ngstype}_dv_gvcf/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  find deepvariant/vcf/ -type f -name "*.dv.phased.vcf.gz*" | parallel -j 8 'destination_file="/data/OGL/resources/OGLsample/${ngstype}_dv_vcf/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  ' # These are hard-filtered and phased.
  find clair3/gvcf -type f -name "*.gvcf.gz*" | parallel -j 8 'destination_file="/data/OGL/resources/OGLsample/${ngstype}_clair3_gvcf/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  find clair3/vcf -type f -name "*.filtered.vcf.gz*" | parallel -j 8 'destination_file="/data/OGL/resources/OGLsample/${ngstype}_clair3_vcf/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  find freebayes/vcf -type f -name "*.phased.vcf.gz*" | parallel -j 8 'destination_file="/data/OGL/resources/OGLsample/${ngstype}_freebayes_vcf/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  find manta -type f -name "*.*" | parallel -j 8 'destination_file="/data/OGL/resources/manta/${ngstype}/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  find scramble_anno -type f -name "*.tsv" | parallel -j 8 'destination_file="/data/OGL/resources/scramble/exome/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  ;;
 *)
  ngstype="genome"
  find prioritization/gemini_tsv_filtered/ -type f -name "*.tsv.gz" | parallel -j 8 'destination_file="/data/OGL/resources/GeneSearch/$ngstype/gemini_tsv_filtered/$(echo {/} | cut -d. -f 1).tsv.gz"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  find prioritization/ -type f -name "*.gt3.anno3.dvg.vcf.gz*" | parallel -j 8 'destination_file="/data/OGL/resources/OGLsample/annotatedVCF/$ngstype/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  find clinSV/ -type f \( -iname '*.clinsv.sv-cnv.pass.vcf.gz*' \
  -o -iname '*.clinsv.sv-cnv.rare_pass_gene.vcf.gz*' \
  -o -iname '*.clinsv.rare_pass_gene.annotated.tsv.gz' \
  -o -iname '*.clinSV.RARE_PASS_GENE.eG.tsv' \) \
  | parallel -j 8 'destination_file="/data/OGL/resources/clinSV/genome/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  find deepvariant/gvcf/ -type f -name "*.vcf.gz*" | parallel -j 8 'destination_file="/data/OGL/resources/OGLsample/${ngstype}_dv_gvcf/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  find deepvariant/vcf/ -type f -name "*.dv.filtered.vcf.gz*" | parallel -j 8 'destination_file="/data/OGL/resources/OGLsample/${ngstype}_dv_vcf/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  ' # these are hard-filtered un-phased.
  find clair3/gvcf -type f -name "*.gvcf.gz*" | parallel -j 8 'destination_file="/data/OGL/resources/OGLsample/${ngstype}_clair3_gvcf/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  find clair3/vcf -type f -name "*.filtered.vcf.gz*" | parallel -j 8 'destination_file="/data/OGL/resources/OGLsample/${ngstype}_clair3_vcf/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  find freebayes/vcf -type f -name "*.phased.vcf.gz*" | parallel -j 8 'destination_file="/data/OGL/resources/OGLsample/${ngstype}_freebayes_vcf/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  find manta -type f -name "*.*" | parallel -j 8 'destination_file="/data/OGL/resources/manta/${ngstype}/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  find jax-cnv -type f -name "*.*" | parallel -j 8 'destination_file="/data/OGL/resources/jaxCNV/genome/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  find orf15/vcf -type f -name "*.*" | parallel -j 8 'destination_file="/data/OGL/resources/OGLsample/orf15_clair3/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  find scramble_anno -type f -name "*.tsv" | parallel -j 8 'destination_file="/data/OGL/resources/scramble/genome_mei/{/}"
  cp {} "$destination_file"
  grp=$(stat -c %G -- "$destination_file" 2>/dev/null)
  [ "$grp" = "OGL" ] || chgrp OGL -- "$destination_file"
  '
  ;;
esac

rm -rf coverage/mean.coverage.done
rm -rf bcmlocus/combine.bcmlocus.done
rm -rf orf15/combine.orf15.done
rm -rf mutserve/haplocheck.done
rm -rf freebayes/freebayes.merge.done* 
rm -rf prioritization/dv_fb.merge.done
rm -rf fastqc/multiqc.done
rm -rf .snakemake
rm -rf 00log
rm -rf slurm*.out
rm -rf sample_bam
rm -rf deepvariant/deepvariant_phased_glnexusVcf.merge.done
rm -rf deepvariant/deepvariant.glnexus.phased.merge.done
rm -rf deepvariant/deepvariant.gvcf.merge.done
rm -rf deepvariant/deepvariantVcf.merge.done
rm -rf clair3/clair3.merge.done
find fastqc -name ".snakemake_timestamp" -exec rm {} \;
rm -rf freebayes/chr_split
rm -rf prioritization/slurm*.out
rm -rf prioritization/.snakemake
rm -rf prioritization/00log
rm -rf prioritization/madeline/madeline.done
rm -rf prioritization/temp
rm -rf prioritization/*.gt3.anno3.dvg.vcf.gz*
rm -rf prioritization/*.gemini.db

#rm -rf old_bam bam cram
echo "File deletion task done"

#copy index case files to resources and change group ownership to OGL
#make indexSample files

# Rscript ~/git/variant_prioritization/dev/pick_affected_index_sample_from_ped_columnName.R prioritization/$ped_file /data/OGL/resources/OGLsample/OGL.indexSample.$TIMESTAMP.ped $datatype.$analysis_batch_name.indexSample.tsv

# tail -n 2+ /data/OGL/resources/OGLsample/OGL.indexSample.$TIMESTAMP.ped >> /data/OGL/resources/OGLsample/OGL.indexSample.ped

#some scripts are in Z:\resources\OGLsample
#Steps: remove duplicated samples and others "special cases", copy vcf files and SV tsv files to resources, merge samples in another scripts etc. These special cases can be in a file.
#test later.
