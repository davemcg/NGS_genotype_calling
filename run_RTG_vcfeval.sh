#!/bin/bash

vcf_to_compare=$1
output_dir=$2

/home/mcgaugheyd/rtg/rtg-core-non-commercial-3.6.2/./rtg RTG_MEM=16G \
vcfeval --threads=$SLURM_CPUS_PER_TASK --ref-overlap \
--baseline=/data/mcgaugheyd/projects/nei/mcgaughey/NA12878/truth/NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.18_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.reg.chr.vcf.gz \
--bed-regions=/data/mcgaugheyd/projects/nei/mcgaughey/NA12878/truth/v2.19_intersect_TruSeq_exome_targeted_regions.hg19.bed \
-c $1 -o $2 \
-t /home/mcgaugheyd/rtg/rtg-core-non-commercial-3.6.2/hg19.sdf
