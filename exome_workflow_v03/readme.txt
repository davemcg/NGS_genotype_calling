This workflow is similar to exome_workflow_v01, but instead of using hard filters on the variants, VariantRecalibrator is used to calculate VQSLOD to more precisely remove bad variants.

This requires >30 human exomes. As of 2015/03/15, I have 20 ddl exomes and 27 ccgo exomes. Plus a few other random groups (W23) and (fam5/34). So, plenty to use VariantRecalibrators.

The first stage (bwa-gatk_workflow_v1.sh) is identical.

The second stage starts out the same: use GenotypeGVCFs to merge the GVCFs together. Then it differs.

 
