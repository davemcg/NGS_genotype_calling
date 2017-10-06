# Purpose
Creating a new VCF from ega-archive.org WGS/WES/NGS data

# Background
Specifically for EGAD00001002656. But this should generally work for all EGA NGS data. They are delivered as [CRAM](http://www.internationalgenome.org/faq/what-are-cram-files/) files. One per person/sample. 

# Outline
1. CRAM -> lane specific BAM
2. Realign to GRCh37 or GRCh38 with bwa-mem
3. Combine lane specific BAM to sample specific BAM
4. GATK to call GVCF
5. Cohort GVCF to VCF
6. Annotate

# Notes
4, 5, and 6 are already taken care of by the various scripts in this repository.  
