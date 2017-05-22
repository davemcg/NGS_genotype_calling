# Workflows to process raw sequence into VCF
- Forked over from https://github.com/davemcg/biowulf2-bin

# For NISC exomes/WGS
- Refer to readme.md in NISC_workflow/

# Overview
1. (Re)Align with bwa-mem against 1000G GRCh37 with decoy
- NGS/generic/workflow/run_bwa-mem_hg37d5_fromBam.sh (if starting from a bam)
- NGS/generic/workflow/run_bwa-mem_hg37d5.sh (if starting from fastq)
2. Call GVCF with GATK (v3.5-0 right now)
- src/process_and_callGVCF.sh
3. Filter GVCF to VCF
- src/GVCF_to_hardFilteredVCF.sh (first choice)
- src/GVCF_to_VQSRfilteredVCF.sh (if WGS or >30 exomes)
4. Annotate variants
- New repository: https://github.com/davemcg/variant_prioritization
